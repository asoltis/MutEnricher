from __future__ import division
import sys, os
import argparse
import numpy as np
from cyvcf2 import VCF
import multiprocessing as mp
import cPickle 
import time
from collections import Counter
import random
from scipy.stats.mstats import gmean
from scipy.stats import chi2 as chi2
from math_funcs import WAP
from scipy.special import betainc
from region_covariate_clustering import covariate_cluster as covc
from pysam import TabixFile
import gzip

'''
MutEnricher non-coding analysis module code.
Created by Anthony R. Soltis (anthony.soltis.ctr@usuhs.edu)
'''

#######
# RUN #
#######
def run(parser, args, version):
    
    print '------------------------MUTENRICHER NON-CODING------------------------'
    print 'MutEnricher version: %s'%(version)
  
    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################
    
    vargs = vars(args)

    # Parse required
    regions_fn = vargs['regions']
    if not os.path.isfile(regions_fn): parser.error('Regions file does not exist!')
    list_fn = vargs['vcfs']
    if not os.path.isfile(list_fn): parser.error('VCF list file does not exist!')
    for line in open(list_fn):
        fn,sn = tuple(line.strip().split('\t'))
        if not os.path.isfile(fn) or not os.path.isfile(fn+'.tbi'): parser.error('VCF or index file missing for sample %s.'%(sn))

    # Set VCF data # 
    VCFs,sample_names = [],[]
    for line in open(list_fn).readlines():
        v,n = line.strip().split('\t')
        VCFs.append(v)
        sample_names.append(n)
    ns = len(VCFs) # Number of samples
    if len(set(sample_names)) < ns:
        print 'WARNING: non-unique sample IDs provided in VCF list file! Exiting.'
        sys.exit()

    # OPTIONS
    print '\n----------------------------INITIALIZATION----------------------------'
    outdir = vargs['outdir']
    os.system('mkdir -p %s'%(outdir))
    print 'Output directory for results: %s'%(outdir)
    if not outdir.endswith('/'): outdir+='/'
    prefix = vargs['prefix'] + '_'
    print 'Analysis prefix: %s'%(prefix)
    nprocessors = min(mp.cpu_count(),vargs['nprocessors'])
    snps_only = vargs['snps_only']
    # Mappability information
    mapr = vargs['map_regions']
    if mapr != None:
        if not os.path.isfile(mapr) or not os.path.isfile(mapr+'.tbi'): parser.error('Mappability file and/or its index not found!')
        print 'Loaded mappable regions from input BED file.'
    # covariates data
    use_covars = False
    cov_fn = vargs['cov_fn']
    if cov_fn != None:
        if not os.path.isfile(cov_fn): parser.error('Supplied covariates file does not exist!')
        else: use_covars = True 
    weights_fn = vargs['weights_fn']
    if weights_fn != None and not os.path.isfile(weights_fn): parser.error('Supplied covariates weight file does not exist!')
    if cov_fn != None and weights_fn != None:
        cov_fn_covs = open(cov_fn).readlines()[0].strip().split('\t')[1:]
        cv_count = 0 
        for l in open(weights_fn): 
            if l.strip().split('\t')[0] not in cov_fn_covs:
                parser.error('%s in covariates weight file not in covariate file!'%(l.strip().split('\t')[0]))
            cv_count += 1
        if cv_count != len(cov_fn_covs):
            parser.error('Mismatch between covariate numbers in covariates file and weights file!')
    clust_precompute = False
    if vargs['cov_precomp_dir'] != None:
        clust_precompute = True
        use_covars = True
        if not vargs['cov_precomp_dir'].endswith('/'): vargs['cov_precomp_dir']+='/'
    min_rclust_size = vargs['min_rclust_size']
    # check if local background requested
    use_local = False
    if vargs['use_local'] == True and use_covars == True:
        print '  --use-local selected but covariates provided. Skipping this option and using covariates for backgrounds.'
        use_local = False
    elif vargs['use_local'] == True and use_covars == False:
        use_local = True
    # hotspot options
    max_hs_dist = vargs['max_hs_dist']
    min_clust_vars = vargs['min_clust_vars']
    no_wap = vargs['no_wap']
    # blacklist
    bl_fn = vargs['blacklist_fn']
    if bl_fn != None:
        if not os.path.isfile(bl_fn): parser.error('Blacklist file does not exist!')
    # Affinity propagation options
    ap_iters = vargs['ap_iters']
    ap_convits=vargs['ap_convits']
    if ap_convits > ap_iters: parser.error('AP total iterations must be greater than convergence iterations! Exiting')
    ap_alg = vargs['ap_alg']
    if ap_alg not in ['fast','slow']: 
        parser.error("AP algorithm type must be one of either 'fast' or 'slow'.")
    
    # create log file for parameters
    oflog = open(outdir + prefix[:-1] + '.log','w')
    oflog.writelines('%s: %s\n'%('MutEnricher version',version))
    oflog.writelines('%s: %s\n'%('Regions.bed',vargs['regions']))
    oflog.writelines('%s: %s\n'%('VCFs',vargs['vcfs']))
    if use_covars: oflog.writelines('background_type: covariate\n')
    elif use_local: oflog.writelines('background_type: local\n')
    else: oflog.writelines('background_type: global\n')
    for param in sorted(vargs.keys()):
        if param in ['regions','vcfs']: continue
        oflog.writelines('%s: %s\n'%(param,str(vargs[param])))
    oflog.writelines('\n')

    # Timing
    time0_master = time.time()

    # initialize parallel pool
    pool = mp.Pool(processes=nprocessors)
    print 'Set pool with %d processors'%(nprocessors)

    #############
    # EXECUTION #
    #############

    # Load in regions
    print '\n----------------------------LOADING REGIONS---------------------------'
    print 'Loading regions...'
    regions = load_regions(regions_fn,mapr)
    # Create region2index dictionary
    reg2index = {}
    for r in regions:
        ind = r.index
        rstr = r.region_string
        reg2index[rstr] = ind
    print 'Loaded %d regions from input BED file.'%(len(regions))
    oflog.writelines('Loaded %d regions from input BED file.\n'%(len(regions)))

    # Load blacklist if provided
    if bl_fn != None:
        print '\nLoading blacklist variants file...'
        blacklist = load_blacklist(bl_fn)
        print 'Blacklist variants loaded.'
    else: blacklist=None

    # Regional clustering 
    if use_covars:
        print '\n----------------------REGION COVARIATE CLUSTERING----------------------'
        time0_cov = time.time()
        region_clusters = {}
        if not clust_precompute:
            cpath = outdir+'apcluster_regions/'
            os.system('mkdir -p %s'%(cpath))
            region_clusters = covc(cov_fn,weights_fn,cpath,ap_iters,ap_convits,ap_alg,pool)
        else:
            print 'Using pre-computed clusters.'
            cpath = vargs['cov_precomp_dir']
            #contigs = os.listdir(cpath)
            contigs = [x for x in os.listdir(cpath) if os.path.isdir(cpath+x)]
            for c in contigs:
                region_clusters[c] = []
                for line in open(cpath+c+'/clusters.txt').readlines():
                    l = line.strip().split('\t')
                    num,exemplar,members = l[0],l[1],l[2]
                    region_clusters[c].append([num,exemplar,members.split(';')])
        
        # log cluster info
        tot_r_from_covars, tot_clusters = 0, 0
        for c in region_clusters:
            for cl in region_clusters[c]:
                tot_clusters += 1
                tot_r_from_covars += len(cl[2])
        oflog.writelines('Clustered %d regions with covariates into %d clusters.\n'%(tot_r_from_covars,tot_clusters))

        # Map to regions
        for contig in region_clusters:
            for clust in region_clusters[contig]:
                members = clust[2]
                indexes = []
                for m in members: 
                    indexes.append(reg2index[m])
                for ind in indexes:
                    r = regions[ind]
                    for ind2 in indexes:
                        if ind2 == ind: continue
                        r.covariate_clusters.append(ind2)
        
        # Check for regions with too small clusters
        small_clust_regions = set()
        for r in regions:
            if len(r.covariate_clusters) < min_rclust_size:
                small_clust_regions.add(reg2index[r.region_string])
        
        # Calculate local background rate for small cluster regions
        print '\nCalculating local background rates for %d regions with fewer than %d'%(len(small_clust_regions),min_rclust_size)
        print 'cluster members.'
        chunk_size_local = 50
        num_chunks_local = int(np.ceil(len(small_clust_regions)/chunk_size_local))
        chunks_local = []
        scr = sorted(list(small_clust_regions))
        for i in range(0,num_chunks_local):
            rr = range(i*chunk_size_local,min((i+1)*chunk_size_local,len(scr)))
            chunks_local.append([regions[index] for index in scr[rr[0]:rr[-1]+1]]) 
        res = [pool.apply_async(get_local_background,args=(chunk,VCFs,sample_names,snps_only,blacklist,True,mapr)) for chunk in chunks_local]
        for r in res:
            rget = r.get()
            for rr in rget:
                assert regions[rr.index].index == rr.index
                regions[rr.index] = rr
        print '\nRegion covariate clustering analysis complete (%0.2f min.).'%((time.time() - time0_cov)/60)

    pool.close()
    pool = mp.Pool(processes=nprocessors)

    # Extract mutation info from VCFs
    print '\n---------------------------MUTATION COUNTING---------------------------'
    print 'Extracting mutations from %d VCF files...'%(len(VCFs))
    time0_mutc = time.time()
    dones,was_done = [],0
    chunk_size = 1000
    num_chunks = int(np.ceil(len(regions)/chunk_size))
    chunks = []
    for i in range(0,num_chunks):
        rr = range(i*chunk_size,min((i+1)*chunk_size,len(regions)))
        chunks.append(regions[rr[0]:rr[-1]+1])
    print '  Divided %d regions into %d region chunks.'%(len(regions),num_chunks)
    res = [pool.apply_async(count_mutations_from_vcfs,args=(VCFs,sample_names,c,snps_only,blacklist,mapr),callback=dones.append) for c in chunks]
    while len(dones) != num_chunks:
        if len(dones)%10==0 and was_done != len(dones):
            was_done = len(dones)
            print '    %d of %d region chunks complete.'%(len(dones),num_chunks)
    regions = []
    for r in res:
        rget = r.get()
        for rr in rget: regions.append(rr)
    regions = sorted(regions,key=lambda x:x.index) # re-sort new regions by index
    print 'Done extracting mutations in VCF files (%0.2f min.).'%((time.time() - time0_mutc)/60)
    
    # Log mutations
    total_muts, tot_r_with_mutation = 0, 0
    for r in regions:
        if r.num_mutations > 0:
            total_muts += r.num_mutations
            tot_r_with_mutation += 1
    oflog.writelines('%d total somatic mutations identified in %d regions.\n'%(total_muts, tot_r_with_mutation))

    # close and re-open pool
    pool.close()
    pool = mp.Pool(processes=nprocessors)
    
    # Get global per-sample background rates
    print '\n------------------------BACKGROUND CALCULATIONS------------------------'
    time0_gbg = time.time()
  
    # Get cluster background
    if use_covars:
        print 'Calculating covariate background rates...'
        for r in regions: 
            if r.num_mutations == 0: continue
            get_cluster_bg_rates(r,regions,scr)
    elif use_local:
        print 'Calculating local region background rates...'
        chunk_size = 100
        regs_for_local = [r for r in regions if r.num_mutations > 0] 
        num_chunks = int(np.ceil(len(regs_for_local)/chunk_size))
        chunks = []
        for i in range(0,num_chunks):
            rr = range(i*chunk_size,min((i+1)*chunk_size,len(regs_for_local)))
            chunks.append(regs_for_local[rr[0]:rr[-1]+1])
        res = [pool.apply_async(get_local_background,args=(chunk,VCFs,sample_names,snps_only,blacklist,False,mapr)) for chunk in chunks]
        for r in res:
            rget = r.get()
            for rr in rget:
                assert regions[rr.index].index == rr.index
                regions[rr.index] = rr
        print 'Local backgrounds obtained.'
    else:
        print 'Calculating global per-sample region background rates...'
        global_bg_rates = get_global_bg_rates(regions,sample_names)
    print 'Background rates calculated (%0.2f min.).'%((time.time()-time0_gbg)/60)
   
    print '\n-----------------------CANDIDATE HOTSPOT FINDING-----------------------'
    print 'Finding candidate hotspots in regions...'
    time0_hs = time.time()
    tot_r_with_hs, tot_hotspots = 0, 0
    for r in regions:
        if r.num_mutations < min_clust_vars: continue
        r.find_clusters(dist=max_hs_dist,min_clust_vars=min_clust_vars)
        if len(r.clusters) > 0:
            tot_r_with_hs += 1
            tot_hotspots += len(r.clusters)
    oflog.writelines('Identified %d candidate hotspots in %d regions for testing.\n'%(tot_hotspots, tot_r_with_hs))
    print 'Candidate hotspot regions obtained (%0.2f min.).'%((time.time()-time0_hs)/60)
   
    # Calculate region and hotspot enrichment p-values
    print '\n--------------------------ENRICHMENT ANALYSES--------------------------'
    print 'Calculating region and hotspot enrichment p-values with negative binomial tests...'
    print '\nPerforming regional enrichment analysis...'
    time0_enr = time.time()
    if use_covars:
        print '  Using clustered covariate background regions.'
        res = [pool.apply_async(get_region_enrichments_covar,args=(r,ns)) for r in regions if r.num_mutations>0]
        rget = [r.get() for r in res]
        #for r in rget:
        #    assert regions[r.index].index == r.index
        #    regions[r.index] = r
        for r in rget:
            index,pvf = r
            regions[index].region_pval = pvf
    elif use_local:
        print '  Using local background rates.'
        res = [pool.apply_async(get_region_enrichments_local,args=(r,ns)) for r in regions if r.num_mutations>0]
        rget = [r.get() for r in res]
        #for r in rget:
        #    assert regions[r.index].index == r.index
        #    regions[r.index] = r
        for r in rget:
            index,bgprob,bgtype,pvf = r
            regions[index].background_prob = bgprob
            regions[index].enrichment_bg_type = bgtype
            regions[index].region_pval = pvf
    else:
        print '  Using global background rates.'
        res = [pool.apply_async(get_region_enrichments_global,args=(r,global_bg_rates,ns)) for r in regions if r.num_mutations>0]
        rget = [r.get() for r in res]
        #for r in rget:
        #    assert regions[r.index].index == r.index
        #    regions[r.index] = r
        for r in rget:
            index,bgprob,bgtype,pvf = r
            regions[index].background_prob = bgprob
            regions[index].enrichment_bg_type = bgtype
            regions[index].region_pval = pvf
    print '\nRegion enrichment p-values obtained (%0.2f min.).'%((time.time() - time0_enr)/60)
    
    # get hotspot enrichments
    print '\nCalculating hotspot enrichment p-values with negative binomial tests...'
    time0_hse = time.time()
    if use_covars:
        print '  Using clustered covariate background regions.'
        clust_output = []
        for r in regions:
            if len(r.clusters) > 0:
                oce = get_hotspot_enrichments_covar(r,regions,ns,scr)
                clust_output.append(oce) 
    elif use_local:
        print '  Using local background rates.'
        clust_results = [pool.apply_async(get_hotspot_enrichments_local,args=(r,ns)) for r in regions if len(r.clusters)>0]
        clust_output = [p.get() for p in clust_results]
    else:
        print '  Using global background rates.'
        clust_results = [pool.apply_async(get_hotspot_enrichments_global,args=(r,global_bg_rates,ns)) for r in regions if len(r.clusters)>0]
        clust_output = [p.get() for p in clust_results]
    clust_out_adjust = []
    for co in clust_output: 
        if co == None: continue
        for c in co: clust_out_adjust.append(c)
    cluster_enrichments = sorted(clust_out_adjust,key=lambda x:x[-5])
    # FDR correct p-values
    pvals = []
    for cr in cluster_enrichments: pvals.append(cr[-5])
    fdrs = fdr_BH(np.array(pvals))
    for i,cr in enumerate(cluster_enrichments): cr.insert(9,fdrs[i])
    print '\nHotspot enrichment p-values obtained (%0.2f min.).'%((time.time() - time0_hse)/60)
     
    # Assign hotspots to their originating regions
    for i,cr in enumerate(cluster_enrichments):
        rindex = reg2index[cr[1]]
        region = regions[rindex]
        assert region.index == rindex
        region.cluster_enrichments.append(cr)
  
    # perform WAP hotspot analysis
    if no_wap:
        print '\nSkipping weighted average proximity (WAP) procedure.'
    else:
        print '\nPerforming weighted average proximity (WAP) hotspot enrichments...'
        time0_wap = time.time()
        dones,was_done = [],0
        res = [pool.apply_async(hotspot_wap_pval,args=(r,),callback=dones.append) for r in regions if r.num_mutations>=min_clust_vars]
        nreg = len([r for r in regions if r.num_mutations>=min_clust_vars])
        while len(dones) < nreg:
            if len(dones)%5000 == 0 and was_done != len(dones):
                was_done = len(dones)
                print '  %d of %d regions complete.'%(was_done,nreg)
        rget = [r.get() for r in res]
        n_wap_tests = 0
        for r in rget:
            assert regions[r.index].index == r.index
            regions[r.index] = r
            n_wap_tests += 1
        print '\nWAP analysis complete (%0.2f min.).'%((time.time() - time0_wap)/60)
        oflog.writelines('Performed WAP analysis on %d regions with >= %d mutations.\n'%(n_wap_tests,min_clust_vars))

    # Compute Fisher combined p-values 
    print '\nCombining region and WAP p-values with Fisher method...'
    for r in regions:
        if r.num_mutations < 1: continue
        elif r.WAP_pval == None: 
            r.WAP = np.nan
            r.WAP_pval = np.nan
            r.fisher_pval = r.region_pval
        else: r.compute_fisher()

    # Write region enrichment output
    print '\nCorrecting p-values and writing regional and WAP enrichment analyses output...'
    ofr = open(outdir+prefix+'region_WAP_enrichments.txt','w') 
    ofr.writelines('\t'.join(['Region','region_name','num_mutations','length','effective_length','bg_type','bg_prob','region_pval',
                              'WAP','WAP_pval','Fisher_pval','FDR_BH','num_samples','position_counts','mutation_counts','samples'])+'\n')
    r_pvals = []
    reg_enr = [r for r in regions if r.num_mutations>0]
    reg_enr.sort(key=lambda x:x.fisher_pval)
    for r in reg_enr:
        r_pvals.append(r.fisher_pval)
    r_fdrs = fdr_BH(np.array(r_pvals))
    tot_r_output = 0
    for i,r in enumerate(reg_enr):
        r.fisher_qval = r_fdrs[i]
        bgp = r.background_prob
        if bgp == None: bgp = -1
        
        rs,name,num,rlen,eff_len,bg_type,bgp = r.region_string,r.name,r.num_mutations,r.length,r.length*ns,r.enrichment_bg_type,r.background_prob
        rpv = r.region_pval
        w0,wap_pv,fish_pv,fdr = r.WAP,r.WAP_pval,r.fisher_pval,r.fisher_qval
        samples = r.mutations_by_sample.keys()
        nsamps,sampstr = len(samples),';'.join(sorted(samples))
        pos_counter = Counter(r.positions)
        pos_counts = []
        for pos in sorted(pos_counter.keys()):
            pos_counts.append('%d_%d'%(pos,pos_counter[pos]))
        pos_counts = ';'.join(pos_counts)
        mut_counts = []
        for mut in sorted(r.samples_by_mutations.keys()):
            count = len(r.samples_by_mutations[mut])
            mut_counts.append('%s_%d'%(mut,count)) 
        mut_counts = ';'.join(mut_counts)

        # write output line 
        ol = '%s\t%s\t%d\t%d\t%d\t%s\t%0.2g\t%0.3g\t%0.2g\t%0.3g\t%0.3g\t%0.3g\t%d\t%s\t%s\t%s\n'%(rs,name,num,rlen,eff_len,
                                                                                                   bg_type,bgp,rpv,w0,wap_pv,fish_pv,fdr,
                                                                                                   nsamps,pos_counts,mut_counts,sampstr)
        ofr.writelines(ol)
        tot_r_output += 1
    ofr.close()
    oflog.writelines('\n%d region enrichment results reported.\n'%(tot_r_output))

    # Write hotspot enrichment output 
    print 'Writing hotspot enrichment analysis output...'
    of = open(outdir+prefix+'hotspot.txt','w')
    of.writelines('\t'.join(['Hotpsot','region','region_name','num_mutations','hotspot_length','effective_length',
                             'bg_type','bg_prob','pval','BH_qval','num_samples','position_counts','mutation_counts','samples'])+'\n')
    for i,cr in enumerate(cluster_enrichments):
        ostr = '%s\t%s\t%s\t%d\t%d\t%d\t%s\t%0.3g\t%0.3g\t%0.3g\t%d\t%s\t%s\t%s\n'%(tuple(cr))
        of.writelines(ostr)
    of.close()
    oflog.writelines('%d region hotspot enrichment results reported.\n'%(len(cluster_enrichments)))

    # Save regions analysis as python pickle object
    cPickle.dump(regions,open(outdir+prefix+'region_data.pkl','w'))

    # close pool, log file
    pool.close()
    oflog.close()

    # Finish statement
    print '\nAnalysis finished in %0.2f minutes.\n'%((time.time() - time0_master)/60)

############
# END MAIN #
############

#####################
# General functions #
#####################
def load_regions(regions_fn,mapr):
    '''
    Load in regions BED file. 
    Returns list of Region class variables.
    '''
    regions = []
    if mapr != None: mapr = TabixFile(mapr)
    
    if regions_fn.endswith('.gz'): FH = gzip.open(regions_fn,'rb')
    else: FH = open(regions_fn)
    for i,r in enumerate(FH): # testing
        region = Region(r,i,mapr)
        regions.append(region)
    return regions

def load_blacklist(bl_fn):
    '''
    Load blacklist regions and return dictionary.
    Requires first four columns to be: contig (chromosome), position, ref base, alt base
    '''
    blacklist = {}
    for line in open(bl_fn).readlines():
        l = line.strip().split('\t')
        chrom,pos,ref,alt = tuple(l[0:4])

        if chrom not in blacklist: blacklist[chrom] = set()
        var_info = '%s_%s_%s'%(pos,ref,alt)
        blacklist[chrom].add(var_info)

    return blacklist

def count_mutations_from_vcfs(VCFs,names,regions,snps_only,blacklist=None,mapr=None):
    '''
    Read in somatic VCF files and extract mutations in regions.
    '''
    for j,vcf_f in enumerate(VCFs):
        vcf = VCF(vcf_f,gts012=True,lazy=False)
        name = names[j]
        tot_added = 0
        for i,r in enumerate(regions):
            if mapr == None:
                r_strings = [r.region_string]
            else: r_strings = r.mappable_regions
            chrom = r.chrom
            tot = 0
            for rstr in r_strings:
                for v in vcf(rstr):
                    # variant filtering
                    filt = v.FILTER
                    if filt != None: continue # None is PASS with cyvcf2
                    
                    # Get variant info
                    pos = v.POS
                    ref = v.REF.encode('ascii')
                    for a in v.ALT:
                        alt = a.encode('ascii')
                        var_info = '%d_%s_%s'%(pos,ref,alt)
                    
                        # Check for duplicate lines
                        # strange examples where same exact variant is reported >1 time
                        if name in r.mutations_by_sample:
                            if var_info in r.mutations_by_sample[name]: continue

                        # check if black-listed site
                        if blacklist != None:
                            if chrom in blacklist:
                                if var_info in blacklist[chrom]: continue

                        if name not in r.mutations_by_sample: r.mutations_by_sample[name] = []
                        r.mutations_by_sample[name].append(var_info)
                        r.positions.append(pos)
                        if pos not in r.samples_by_positions:
                            r.samples_by_positions[pos] = [name]
                        else: r.samples_by_positions[pos].append(name)
                        if var_info not in r.samples_by_mutations:
                            r.samples_by_mutations[var_info] = [name]
                        else: r.samples_by_mutations[var_info].append(name)
                        tot+=1
            r.num_mutations += tot
            tot_added += tot
   
    return regions

def get_local_background(regions,VCFs,names,snps_only,blacklist=None,for_covar=False,mapr=None):
    # Windows
    wins = [100e3,500e3,1e6]

    # load mappable regions file if supplied
    if mapr != None:
        mappable = TabixFile(mapr)
        
    # Loop over VCFs
    for j,vcf_f in enumerate(VCFs):
        vcf = VCF(vcf_f,gts012=True,lazy=False)
        name = names[j]

        for r in regions:
            if for_covar == False: # need all computations for covariate analysis; otherwise just need for mutated samples
                if name not in r.mutations_by_sample: continue # skip if no mutations in region for current sample
            lb = {}
            bgs_per_win = []
            for win in wins:
                half_win = win/2
                midp = int(round((r.start+r.stop)/2))
                start = max(1,midp-half_win)
                stop = midp + half_win
                region_string = '%s:%d-%d'%(r.chrom,start,stop)
                
                # Get mappable regions in window (if supplied)
                r_strings = []
                win_l = 0
                if mapr == None: 
                    r_strings.append(region_string)
                    win_l += (stop-start)+1
                else:
                    tbq = [t.split('\t') for t in mappable.fetch(region_string)]
                    for t in tbq:
                        qstart,qstop = int(t[1])+1,int(t[2])
                        if qstart < start: ostart = start
                        else: ostart = qstart
                        if qstop > stop: ostop = stop
                        else: ostop = qstop
                        win_l += ostop - ostart + 1
                        MR = '%s:%d-%d'%(r.chrom,ostart,ostop)
                        r_strings.append(MR)
                
                if win_l == 0:
                    bgs_per_win.append(0)
                    continue

                # Count mutations
                mut_count = 0
                for rstr in r_strings:
                    for v in vcf(rstr):
                        # variant filtering
                        filt = v.FILTER
                        if filt != None: continue # None is PASS with cyvcf2
                         
                        # Get variant info
                        for a in v.ALT:
                            # check if black-listed site
                            if blacklist != None:
                                if r.chrom in blacklist:
                                    var_info = '%d_%s_%s'%(v.POS,v.REF.encode('ASCII'),a.encode('ascii'))
                                    if var_info in blacklist[r.chrom]: continue
                            mut_count += 1
                
                bgs_per_win.append(mut_count / win_l)

            # Get max rate
            maxi,maxr = 0,bgs_per_win[0]
            #print bgs_per_win
            for i,bgr in enumerate(bgs_per_win):
                if bgr > maxr:
                    maxi,maxr = i,bgr

            # Set current region's local background for current VCF sample
            r.local_backgrounds[name] = {}
            r.local_backgrounds[name]['win_size'] = wins[maxi]
            r.local_backgrounds[name]['bg_rate'] = maxr

    return regions

def get_global_bg_rates(regions,names):
    '''
    Determine per-sample global background rates by counting mutations in each sample at each region.
    '''
    # Initialize output dictionary
    bg_rates = {}
    for n in names: bg_rates[n] = [0,0] # Initialize with two-element list for counts and total length

    # Get background rates
    for r in regions:
        rl = r.length
        for s in names:
            try:
                nm = len(r.mutations_by_sample[s])
            except KeyError: nm = 0
            bg_rates[s][0] += nm
            bg_rates[s][1] += rl

    # Calculate rates and return
    for s in bg_rates:
        rate = bg_rates[s][0] / bg_rates[s][1]
        bg_rates[s] = rate
    return bg_rates

def get_cluster_bg_rates(r,regions,scr):
    '''
    Determine per-region background rates from clustered regions. 
    If region is a member of small cluster, the local background rate is used.
    '''
    # Get appropriate background
    csf = [s for s in r.mutations_by_sample]
    bgpf_l = []
    if r.index in scr:
        r.enrichment_bg_type = 'local'
        for s in csf: bgpf_l.append(r.local_backgrounds[s]['bg_rate'])
    else:
        r.enrichment_bg_type = 'clustered_regions'
        members = [r.index] + r.covariate_clusters
        regs = [regions[i] for i in members]
        bgf_rates = get_global_bg_rates(regs,csf)
        for s in csf: bgpf_l.append(bgf_rates[s])
    bgpf = gmean(bgpf_l) # geometric mean
    r.background_prob = bgpf

def get_region_enrichments_global(r,bg_rates,ns):
    '''
    Function for determining full region enrichments using global background rates.
    Return tuple data.
    '''
    # calculate p-value for full region
    xf = r.length * ns
    kf = r.num_mutations
    rsamps = r.mutations_by_sample.keys()
    bg_l = []
    for s in rsamps: bg_l.append(bg_rates[s])
    bg = gmean(bg_l)
    
    # Calculate full region p-value
    try: 
        pvf = betainc(kf,xf-kf+1,bg)
    except:
        print '  error at region %s with length: %d, num mutations: %d, and bg: %f'%(r.region_string,r.length,kf,bg)
        pvf = 1

    # Return tuple
    return (r.index,bg,'global',pvf)

def get_region_enrichments_local(r,ns):
    '''
    Function for determining full region enrichments using local background rates.
    Return tuple data.
    '''
    # calculate p-value for full region
    xf = r.length * ns
    kf = r.num_mutations
    rsamps = r.mutations_by_sample.keys()
    bg_l = []
    for s in rsamps: bg_l.append(r.local_backgrounds[s]['bg_rate'])
    bg = gmean(bg_l)
    
    # Calculate full region p-value
    try: 
        pvf = betainc(kf,xf-kf+1,bg)
    except:
        print '  error at region %s with length: %d, num mutations: %d, and bg: %f'%(r.region_string,r.length,kf,bg)
        pvf = 1
    
    # Return tuple
    return (r.index,bg,'local',pvf)

def get_region_enrichments_covar(r,ns):
    '''
    Function for determining full region enrichments using clustered covariate background rates.
    Return tuple data.
    '''
    # calculate p-value for full region
    xf = r.length * ns
    kf = r.num_mutations #- 1
    bg = r.background_prob

    # Calculate full region p-value
    try: 
        pvf = betainc(kf,xf-kf+1,bg)
    except:
        print '  error at region %s with length: %d, num mutations: %d, and bg: %f'%(r.region_string,r.length,kf,bg)
        pvf = 1

    # Return tuple
    return (r.index,pvf)

def get_hotspot_enrichments_covar(r,regions,ns,scr):
    '''
    Compute hotspot enrichment p-values using negative binomial test with clustered-region background rates.
    '''
    # initialize output
    enrich = []
    # Find clusters and compute p-values
    if len(r.clusters) > 0:
        for clust in r.clusters:       
            # Cluster info
            clust_name = r.clusters[clust]['name']

            # data for NB test 
            k = r.clusters[clust]['count']
            c_len = r.clusters[clust]['length']
            x = c_len * ns # cluster length times number of samples
            cs = r.clusters[clust]['samples']
            
            # Get background            
            bgp_l,bg_type = [],''
            if r.index in scr:
                bg_type = 'local'
                for s in cs: bgp_l.append(r.local_backgrounds[s]['bg_rate'])
            else:
                bg_type = 'clustered_regions'
                members = [r.index] + r.covariate_clusters
                regs = [regions[i] for i in members]
                bg_rates = get_global_bg_rates(regs,cs)  
                for s in cs: bgp_l.append(bg_rates[s]) 
            bgp = gmean(bgp_l) # geometric mean

            # Calculate p-values
            pv = betainc(k,x-k+1,bgp)

            # Get counts for positions and mutations
            counter = Counter(r.clusters[clust]['positions'])
            count_s = []
            for pos in sorted(counter.keys()):
                count_s.append('%d_%d'%(pos,counter[pos]))
            count_s = ';'.join(count_s)
            muts = []
            for mut in sorted(r.samples_by_mutations.keys()):
                pos = int(mut.split('_')[0])
                if pos in counter:
                    for s in r.samples_by_mutations[mut]:
                        if s in cs: muts.append(mut)
            mut_counter = Counter(muts)
            count_m = []
            for mut in sorted(mut_counter.keys()):
                count_m.append('%s_%d'%(mut,mut_counter[mut]))
            count_m = ';'.join(count_m)

            # Output info
            ol = [clust_name,r.region_string,r.name,k,c_len,x,bg_type,bgp,pv,len(cs),count_s,count_m,';'.join(sorted(list(cs)))]
            enrich.append(ol)
        
        # return
        if len(enrich) > 0:
            return enrich

def get_hotspot_enrichments_local(r,ns):
    '''
    Compute hotspot enrichment p-values using negative binomial test with local background rates.
    '''
    # initialize output
    enrich = []
    # Find clusters and compute p-values
    if len(r.clusters) > 0:
        for clust in r.clusters: 
            # Cluster info
            clust_name = r.clusters[clust]['name']

            # data for NB test 
            k = r.clusters[clust]['count']
            c_len = r.clusters[clust]['length']
            x = c_len * ns # cluster length times number of samples
            cs = r.clusters[clust]['samples']
            
            # Get background            
            bgp_l = []
            bg_type = 'local'
            for s in cs: bgp_l.append(r.local_backgrounds[s]['bg_rate'])
            bgp = gmean(bgp_l) # geometric mean

            # Calculate p-values
            pv = betainc(k,x-k+1,bgp)

            # Get counts for positions and mutations
            counter = Counter(r.clusters[clust]['positions'])
            count_s = []
            for pos in sorted(counter.keys()):
                count_s.append('%d_%d'%(pos,counter[pos]))
            count_s = ';'.join(count_s)
            muts = []
            for mut in sorted(r.samples_by_mutations.keys()):
                pos = int(mut.split('_')[0])
                if pos in counter:
                    for s in r.samples_by_mutations[mut]:
                        if s in cs: muts.append(mut)
            mut_counter = Counter(muts)
            count_m = []
            for mut in sorted(mut_counter.keys()):
                count_m.append('%s_%d'%(mut,mut_counter[mut]))
            count_m = ';'.join(count_m)

            # Output info
            ol = [clust_name,r.region_string,r.name,k,c_len,x,bg_type,bgp,pv,len(cs),count_s,count_m,';'.join(sorted(list(cs)))]
            enrich.append(ol)
        
        # return
        if len(enrich) > 0:
            return enrich

def get_hotspot_enrichments_global(r,global_bg_rates,ns):
    '''
    Compute hotspot enrichment p-values using negative binomial test with global background rates.
    '''
    # initialize output
    enrich = []
    bg_type = 'global'
    # Find clusters and compute p-values
    if len(r.clusters) > 0:
        for clust in r.clusters:    
            # Cluster info
            clust_name = r.clusters[clust]['name']

            # data for NB test 
            k = r.clusters[clust]['count']
            c_len = r.clusters[clust]['length']
            x = c_len * ns # cluster length times number of samples
            cs = r.clusters[clust]['samples']
            bgp_l = []
            for s in cs: bgp_l.append(global_bg_rates[s])
            bgp = gmean(bgp_l) # geometric mean
            pv = betainc(k,x-k+1,bgp)

            # Get counts for positions and mutations
            counter = Counter(r.clusters[clust]['positions'])
            count_s = []
            for pos in sorted(counter.keys()):
                count_s.append('%d_%d'%(pos,counter[pos]))
            count_s = ';'.join(count_s)
            muts = []
            for mut in sorted(r.samples_by_mutations.keys()):
                pos = int(mut.split('_')[0])
                if pos in counter:
                    for s in r.samples_by_mutations[mut]:
                        if s in cs: muts.append(mut)
            mut_counter = Counter(muts)
            count_m = []
            for mut in sorted(mut_counter.keys()):
                count_m.append('%s_%d'%(mut,mut_counter[mut]))
            count_m = ';'.join(count_m)
  
            # Output info
            ol = [clust_name,r.region_string,r.name,k,c_len,x,bg_type,bgp,pv,len(cs),count_s,count_m,';'.join(sorted(list(cs)))]
            enrich.append(ol)
        
        # return
        if len(enrich) > 0:
            return enrich

#####################
# Class definitions #
#####################
class Region:
    '''
    Region class.
    '''
    def __init__(self,region,index,mappable=None):
        '''
        Initialize region. Parse BED file line for region.
        '''
        # Data from BED file
        r = region.strip('\n').split('\t')
        self.index = index
        self.chrom = r[0]
        self.start = int(r[1])+1 # Add one to change 0-indexing in BED
        self.stop = int(r[2])
        self.region_string = '%s:%d-%d'%(self.chrom,self.start,self.stop)
        if len(r) > 3:
            self.name = r[3]
        else: self.name = 'region_%d'%(index)

        # Get length according to mappable regions (if supplied)
        self.mappable_regions = []
        if mappable == None: 
            self.length = self.stop - self.start + 1
        else:
            length = 0
            tbq = [t.split('\t') for t in mappable.fetch(self.region_string)]
            for t in tbq:
                qstart,qstop = int(t[1])+1,int(t[2])
                if qstart < self.start: ostart = self.start
                else: ostart = qstart
                if qstop > self.stop: ostop = self.stop
                else: ostop = qstop
                length += ostop - ostart + 1
                MR = '%s:%d-%d'%(self.chrom,ostart,ostop)
                self.mappable_regions.append(MR)
            self.length = length

        # set up other values
        self.region_pval = 1
        self.region_qval = 1
        self.background_prob = None
        self.enrichment_bg_type = None 
        self.WAP = None
        self.WAP_pval = None
        self.fisher_pval = 1
        self.fisher_qval = 1
        self.clusters = {} # Set up dictionary for clusters 
        self.cluster_enrichments = [] # List of cluster enrichment outputs for region
        self.mutations_by_sample = {} # Dictionary of samples and the mutations from them
        self.samples_by_mutations = {} # Dictionary of mutations and the samples containing them
        self.samples_by_positions = {} # Dictionary of mutation positions and the samples containing them
        self.num_mutations = 0 # counter for total mutations in region from all samples
        self.positions = [] # List of coordinates within region for mutations
        self.covariate_clusters = [] # Initialize as empty list for values
        self.local_backgrounds = {}
    
    def find_clusters(self,dist=50,min_clust_vars=3):
        '''
        Find candidate mutation clusters by merging mutations within 'dist' of each other.
        '''
        self.positions.sort() # Sort positions
        #pos_to_test = sorted(list(set(self.positions))) # Considers only unique positions
        pos_to_test = self.positions # considering all mutations
        num_pos = len(pos_to_test)
        clust_num,counts = 0,0 

        # Loop over positions and merge if within dist
        cluster,c_samples = [],set()
        for i in range(0,num_pos-1): 
            curr = pos_to_test[i]
            next_p = pos_to_test[i+1]

            # Add current position if at first index
            if i == 0:
                cluster.append(curr)
                counts += 1
                for s in self.samples_by_positions[curr]: c_samples.add(s)

            # Test if neighboring position is within distance
            # If true, add neighbor to current cluster
            if (next_p - curr) <= dist: 
                cluster.append(next_p)
                counts += 1
                for s in self.samples_by_positions[next_p]: c_samples.add(s)
            
                # Case of overlap at last index - write out if of appropriate size
                if i == num_pos-2 and counts >= min_clust_vars:
                    self.clusters[clust_num] = {}
                    self.clusters[clust_num]['positions'] = cluster
                    self.clusters[clust_num]['count'] = counts
                    self.clusters[clust_num]['samples'] = c_samples
                    self.clusters[clust_num]['length'] = cluster[-1] - cluster[0] + 1
                    self.clusters[clust_num]['name'] = '%s:%d-%d'%(self.chrom,cluster[0],cluster[-1])

            # If not true, write out appropriate clusters
            else: 
                # Write out cluster if of appropriate size
                if i < (num_pos - 2) and counts >= min_clust_vars:
                    self.clusters[clust_num] = {}
                    self.clusters[clust_num]['positions'] = cluster
                    self.clusters[clust_num]['count'] = counts
                    self.clusters[clust_num]['samples'] = c_samples
                    self.clusters[clust_num]['length'] = cluster[-1] - cluster[0] + 1
                    self.clusters[clust_num]['name'] = '%s:%d-%d'%(self.chrom,cluster[0],cluster[-1])
               
                    # Re-set values
                    clust_num += 1
                    cluster = [next_p] # set to current next position
                    c_samples = set()
                    for s in self.samples_by_positions[next_p]: c_samples.add(s)
                    counts = 1

                # In case at last position
                elif i == num_pos-2 and counts >= min_clust_vars:
                    self.clusters[clust_num] = {}
                    self.clusters[clust_num]['positions'] = cluster
                    self.clusters[clust_num]['count'] = counts
                    self.clusters[clust_num]['samples'] = c_samples
                    self.clusters[clust_num]['length'] = cluster[-1] - cluster[0] + 1
                    self.clusters[clust_num]['name'] = '%s:%d-%d'%(self.chrom,cluster[0],cluster[-1])
                    
                    if min_clust_vars <= 1:
                        clust_num+=1
                        self.clusters[clust_num] = {}
                        self.clusters[clust_num]['positions'] = next_p
                        self.clusters[clust_num]['count'] = 1
                        self.clusters[clust_num]['length'] = 1
                        self.clusters[clust_num]['name'] = '%s:%d-%d'%(self.chrom,next_p,next_p)
                        c_samples = set()
                        for s in self.samples_by_positions[next_p]: c_samples.add(s)
                        self.clusters[clust_num]['samples'] = c_samples
                    
                # Eliminate if not of sufficient size
                else:
                    cluster = [next_p]
                    c_samples = set()
                    for s in self.samples_by_positions[next_p]: c_samples.add(s)
                    counts = 1
    
    def compute_fisher(self):
        ''' Combine region and WAP p-values with Fisher's method.'''
        pvals = [self.region_pval, self.WAP_pval]
        fpv = fishers_method(pvals)
        self.fisher_pval = fpv

##################
# Math functions #
##################
def hotspot_wap_pval(region,tau=6,theta=0.25,nperm=1000000):
    
    random.seed(1)

    # true WAP
    w0 = WAP(region.positions,tau)

    # Get number of mutations per position
    n_per_p_C = Counter(region.positions)
    n_per_p = []
    for p in n_per_p_C: n_per_p.append(n_per_p_C[p])

    # permute
    rlen = region.length
    k = 0 # counter for number of permuted WAPs >= w0
    pv = 1
    for n in range(0,nperm):
        pos = []
        for i in range(0,len(n_per_p)):

            rp = int(round(rlen * random.random()))
            pos += n_per_p[i] * [rp]

        w = WAP(pos,tau)
        if w >= w0: k+=1
        pv = (k+1) / (n+1+1)
        
        crit = 2 * 1.96 * np.sqrt((n+1-k+1) / ((n+1+3) * (k+1))) # test value for cut-off
        
        if (n+1) > 1000 and crit < theta:
            break

    region.WAP = w0
    region.WAP_pval = pv
    return region

def fishers_method(pvals):
    '''
    Fisher's combined p-value method applied to list of p-values.
    '''
    for i,p in enumerate(pvals):
        if p == 0: pvals[i] = 1e-322 #pseudo
    chi = -2 * sum([np.log(x) for x in pvals])
    pv = 1-chi2.cdf(x=chi,df=2*len(pvals)) # df = 2 x num p-values
   
    return pv

def fdr_BH(pvals):
    '''
    Compute Benjamini-Hochberg FDR q-values on sorted array of p-values.
    '''
    total = pvals.size
    fdrs0 = total*pvals/range(1,total+1)
    fdrs = []
    # Preserve monotonicity
    for i in range(0,total):
        fdrs.append(min(fdrs0[i:]))
    return fdrs


