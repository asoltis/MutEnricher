from __future__ import division
import sys, os
import argparse
from numpy import array, ceil
from cyvcf2 import VCF
import multiprocessing as mp
import cPickle 
import time
from collections import Counter
import random
from scipy.stats.mstats import gmean
from scipy.special import betainc
from gene_covariate_clustering import covariate_cluster as covc
from pysam import TabixFile

'''
MutEnricher coding analysis module code.
Created by Anthony R. Soltis (anthony.soltis.ctr@usuhs.edu)
'''

#######
# RUN #
#######
def run(parser,args,version):

    print '--------------------------MUTENRICHER CODING--------------------------'
    print 'MutEnricher version: %s'%(version)
  
    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################

    vargs = vars(args)
    
    # Parse required
    gtf_fn = vargs['GTF']
    if not os.path.isfile(gtf_fn): parser.error('GTF file does not exist!')
    list_fn = vargs['vcfs']
    genefield = vargs['genefield']
    MAF = vargs['maf'] 
    use_maf = False
    if MAF != None:
        if not os.path.isfile(MAF): parser.error('Supplied MAF file does not exist!')
        else:
            use_maf = True
    ns = None # initialize number of samples variable here, set later
    sample_names = None # initialize sample names variable
    if use_maf == False:
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
    # background variants type
    bg_vtype = vargs['bg_vars_type']
    if bg_vtype not in ['all','silent']: parser.error("--bg-vars-type must be one of either 'all' or 'silent'.")
    if bg_vtype == 'all': print 'Considering all variants in background rate calculations.'
    elif bg_vtype == 'silent': print 'Considering only silent variants in background rate calculations.'
    # annotation type
    tType = vargs['tType']
    if use_maf: tType = 'maf' # over-ride if MAF being used
    elif tType not in ['annovar','illumina']:
        if os.path.isfile(tType):
            annrows = set([x.strip().split('\t')[0] for x in open(tType).readlines()])
            if len(annrows) != 2: parser.error('Custom annotations file must contain Gene and Effect as row names')
            elif 'Effect' not in annrows and 'Gene' not in annrows:
                parser.error('Custom annotations file must only contain Gene and Effect as row names')
            else: print 'Using non-silent annotation terms from input file %s'%(tType)
        else: parser.error('Invalid annotation type!')
    else: print 'Annotation type: %s'%(tType)
    snps_only = vargs['snps_only']
    if snps_only: print 'Performing SNPs-only analysis'
    else: print 'Considering both SNPs and indels in analysis.'
    exome_only = vargs['exome_only']
    if exome_only: print 'Considering only exonic gene coordinates.'
    gene_list = vargs['gene_list']
    genes_to_use = []
    if gene_list != None:
        if not os.path.isfile(gene_list): parser.error('Supplied gene list file does not exist!')
        else:
            genes_to_use = [x.strip() for x in open(gene_list).readlines()]
            genes_to_use = set(genes_to_use)
            print 'Only considering %d genes from input gene list in analysis.'%(len(genes_to_use))
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
    min_clust_size = vargs['min_clust_size']
    by_contig = vargs['by_contig']
    # check if local background requested
    use_local = False
    if vargs['use_local'] == True and use_covars == True:
        print '  --use-local selected but covariates provided. Skipping this option and using covariates for backgrounds.'
        use_local = False
    elif vargs['use_local'] == True and use_covars == False:
        use_local = True
    if use_local == True and use_maf == True:
        print '  --use-local option not available when reading from MAF. Using global bacgkround rates instead.'
        use_local = False
    # hotspot options
    max_hs_dist = vargs['max_hs_dist']
    min_clust_vars = vargs['min_clust_vars']
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

    # create logging file
    oflog = open(outdir + prefix[:-1] + '.log','w')
    oflog.writelines('%s: %s\n'%('MutEnricher version',version))
    oflog.writelines('%s: %s\n'%('GTF',vargs['GTF']))
    oflog.writelines('%s: %s\n'%('VCFs',vargs['vcfs']))
    if use_covars: oflog.writelines('background_type: covariate\n')
    elif use_local: oflog.writelines('background_type: local\n')
    else: oflog.writelines('background_type: global\n')
    for param in sorted(vargs.keys()):
        if param in ['GTF','vcfs']: continue
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
    
    # Load in genes
    print '\n-----------------------------LOADING GENES----------------------------'
    print 'Loading GTF...'
    GTF = load_gtf(gtf_fn,genefield,gene_list,genes_to_use)
    print 'GTF loaded.'
    print 'Loading genes...'
    genes = []
    g_index_counter = 0
    
    if use_maf:
        print ' loading genes from MAF...'
        genes_in_maf = set([x.strip().split('\t')[0] for x in open(MAF).readlines()[1:]]) # Unique genes in MAF
        genes_not_recognized,gnr_count = [],0
        if mapr != None: mappable = TabixFile(mapr)
        else: mappable = None
        for g in genes_in_maf:
            if g not in GTF:
                genes_not_recognized.append(g)
                gnr_count += 1
                continue
            else:
                gene = Gene(GTF[g],g,g_index_counter,mappable,exome_only)
                if gene.coding_length < 50: continue
                else:
                    genes.append(gene)
                    g_index_counter += 1
        print ' %d genes in MAF were not found in GTF and are being skipped.'%(gnr_count)
        print '  writing skipped genes to file: %s%s.'%(outdir,'unrecognized_genes.txt')
        ofgnr = open('%s%s'%(outdir,'unrecognized_genes.txt'),'w')
        for g in genes_not_recognized: ofgnr.writelines(g+'\n')
        ofgnr.close()
        oflog.writelines('Skipped %d genes from input MAF not present in GTF file.\n'%(gnr_count))
    else:
        if mapr != None: mappable = TabixFile(mapr)
        else: mappable = None
        for g in sorted(GTF.keys()):
            gene = Gene(GTF[g],g,g_index_counter,mappable,exome_only)
            if gene.coding_length < 50: continue
            else:
                genes.append(gene)
                g_index_counter += 1
    
    # Create gene2index dictionary
    gene2index = {}
    for g in genes:
        ind = g.index
        name = g.name
        gene2index[name] = ind
    print 'Loaded %d genes from input GTF file.'%(len(genes))
    oflog.writelines('Loaded %d genes from input GTF file.\n'%(len(genes)))

    # load nonsilent terms
    terms = load_nonsilent_terms(tType)
    
    # Load blacklist if provided
    if bl_fn != None:
        print '\nLoading blacklist variants file...'
        blacklist = load_blacklist(bl_fn)
        print 'Blacklist variants loaded.'
    else: blacklist=None

    # Gene clustering 
    if use_covars:
        print '\n-----------------------GENE COVARIATE CLUSTERING-----------------------'
        time0_cov = time.time()
        gene_clusters = {}
        if not clust_precompute:
            cpath = outdir+'apcluster_genes/'
            os.system('mkdir -p %s'%(cpath))
            contigs = None
            if by_contig:
                contigs = {}
                for g in genes: contigs[g.name] = g.chrom
            gene_clusters = covc(cov_fn,weights_fn,cpath,ap_iters,ap_convits,ap_alg,pool,contigs)
        else:
            print 'Using pre-computed clusters.'
            cpath = vargs['cov_precomp_dir']
            contigs = [x for x in os.listdir(cpath) if os.path.isdir(cpath+x)]
            for c in contigs:
                gene_clusters[c] = []
                for line in open(cpath+c+'/clusters.txt').readlines():
                    l = line.strip().split('\t')
                    num,exemplar,members = l[0],l[1],l[2]
                    gene_clusters[c].append([num,exemplar,members.split(';')])
    
        # log cluster info
        tot_g_from_covars, tot_clusters = 0, 0
        for c in gene_clusters:
            for cl in gene_clusters[c]:
                tot_clusters += 1
                tot_g_from_covars += len(cl[2])
        oflog.writelines('Clustered %d genes with covariates into %d clusters.\n'%(tot_g_from_covars,tot_clusters))

        # Map to genes
        for contig in gene_clusters:
            for clust in gene_clusters[contig]:
                members = clust[2]
                indexes = []
                for m in members: 
                    if m not in gene2index: continue
                    indexes.append(gene2index[m])
                for ind in indexes:
                    g = genes[ind]
                    for ind2 in indexes:
                        if ind2 == ind: continue
                        g.covariate_clusters.append(ind2)
        
        # Check for geness with too small clusters 
        small_clust_genes = set()
        for g in genes:
            if len(g.covariate_clusters) < min_clust_size:
                small_clust_genes.add(gene2index[g.name])
        if not use_maf:    
            # Calculate local background rate for small cluster genes
            print '\nCalculating local background rates for %d genes with fewer than %d'%(len(small_clust_genes),min_clust_size)
            print 'cluster members.'
            chunk_size_local = 50
            num_chunks_local = int(ceil(len(small_clust_genes)/chunk_size_local))
            chunks_l = []
            scr = sorted(list(small_clust_genes))
            for i in range(0,num_chunks_local):
                rr = range(i*chunk_size_local,min((i+1)*chunk_size_local,len(scr)))
                chunks_l.append([genes[index] for index in scr[rr[0]:rr[-1]+1]])
            res = [pool.apply_async(get_local_background,args=(chunk,VCFs,sample_names,terms,tType,snps_only,blacklist,True,mapr,bg_vtype)) for chunk in chunks_l]
            for r in res:
                rget = r.get()
                for rr in rget:
                    assert genes[rr.index].index == rr.index
                    genes[rr.index] = rr
            del(res)
        else:
            scr = sorted(list(small_clust_genes))
        print '\nGene covariate clustering analysis complete (%0.2f min.).'%((time.time() - time0_cov)/60)
    
    pool.close()
    pool = mp.Pool(processes=nprocessors)

    # Extract mutation info from VCFs
    print '\n---------------------------MUTATION COUNTING---------------------------'
    time0_mutc = time.time()    
    if use_maf:
        print 'Extracting mutations from MAF input file...'
        genes,sample_names = count_mutations_from_maf(MAF,genes,gene2index,terms,tType,snps_only,blacklist)
        ns = len(sample_names) # set number of samples
        oflog.writelines('%d samples detected in input MAF.\n'%(ns))
        
        if use_covars: # Get local bacgkround rates for small cluster size genes after counting
            print '\nCalculating local background rates for %d genes with fewer than %d'%(len(scr),min_clust_size)
            print 'cluster members.'
        
            min_bg = 1e-8
            for i in scr:
                bg_len = genes[i].total_length
                samples = genes[i].mutations_by_sample.keys()
                
                for s in samples:
                    if bg_vtype == 'all':
                        bg_muts = genes[i].mutations_by_sample[s]['bg'] + len(genes[i].mutations_by_sample[s]['nonsilent'])
                    elif bg_vtype == 'silent':
                        bg_muts = genes[i].mutations_by_sample[s]['bg']
                    if bg_muts == 0: bgr = min_bg
                    else: bgr = bg_muts / bg_len
                    genes[i].local_backgrounds[s] = {}
                    genes[i].local_backgrounds[s]['bg_rate'] = bgr 
    else:
        print 'Extracting mutations from %d VCF files...'%(len(VCFs))
        dones,was_done = [],0
        chunk_size = 1000
        num_chunks = int(ceil(len(genes)/chunk_size))
        chunks = []
        for i in range(0,num_chunks):
            rr = range(i*chunk_size,min((i+1)*chunk_size,len(genes)))
            chunks.append(genes[rr[0]:rr[-1]+1])
        print '  Divided %d genes into %d gene chunks.'%(len(genes),num_chunks)

        res = [pool.apply_async(count_mutations_from_vcfs,args=(VCFs,sample_names,c,terms,tType,snps_only,blacklist,mapr),callback=dones.append) for c in chunks]
        while len(dones) != num_chunks:
            if len(dones)%10==0 and was_done != len(dones):
                was_done = len(dones)
                print '    %d of %d gene chunks complete.'%(len(dones),num_chunks)
        genes = []
        for r in res:
            rget = r.get() 
            for rr in rget: genes.append(rr)
        del(res)
        genes = sorted(genes,key=lambda x:x.index) # re-sort new genes by index
        print 'Done extracting mutations from VCF files (%0.2f min.).'%((time.time() - time0_mutc)/60)
    
    # log total mutations
    total_nonsilent,total_muts,tot_g_with_nonsilent = 0,0,0
    for g in genes:
        total_muts += g.total_mutations
        total_nonsilent += g.total_nonsilent_mutations
        if g.total_nonsilent_mutations > 0: tot_g_with_nonsilent += 1
    oflog.writelines('%d total non-silent somatic mutations identified in %d genes.\n'%(total_nonsilent,tot_g_with_nonsilent))
    oflog.writelines('%d total somatic mutations identified in gene limits.\n'%(total_muts))

    # close and re-open pool
    pool.close()
    pool = mp.Pool(processes=nprocessors)
     
    # Global background calculations (if necessary)
    print '\n------------------------BACKGROUND CALCULATIONS------------------------'
    time0_gbg = time.time()
    
    # Get background rates according to option selected
    if use_covars:
        print 'Calculating covariate background rates...'
        # Get global bg rates to deal with zero rate cases
        global_bg_rates = get_global_bg_rates(genes,sample_names,bg_vtype) 
        minBG = 1.0
        for s in global_bg_rates:
            bg = global_bg_rates[s]
            if bg > 0 and bg < minBG: minBG = bg
        for s in global_bg_rates:
            if global_bg_rates[s] == 0: global_bg_rates[s] = minBG
        for g in genes: 
            if g.total_mutations == 0: continue
            get_cluster_bg_rates(g,genes,scr,global_bg_rates,bg_vtype) 
    elif use_local:
        print 'Calculating local gene background rates...'
        chunk_size = 100
        genes_for_local = [g for g in genes if g.total_nonsilent_mutations>0]        
        num_chunks = int(ceil(len(genes_for_local)/chunk_size))
        chunks = []
        for i in range(0,num_chunks):
            rr = range(i*chunk_size,min((i+1)*chunk_size,len(genes_for_local)))
            chunks.append(genes_for_local[rr[0]:rr[-1]+1])
        res = [pool.apply_async(get_local_background,args=(chunk,VCFs,sample_names,terms,tType,snps_only,blacklist,False,mapr,bg_vtype)) for chunk in chunks] 
        for r in res:
            rget = r.get()
            for rr in rget:
                assert genes[rr.index].index == rr.index
                genes[rr.index] = rr
        del(res)
        print 'Local backgrounds obtained.'
    else:
        print 'Calculating global per-sample gene background rates...'
        global_bg_rates = get_global_bg_rates(genes,sample_names,bg_vtype)
        minBG = 1.0
        for s in global_bg_rates:
            bg = global_bg_rates[s]
            if bg > 0 and bg < minBG: minBG = bg
        for s in global_bg_rates:
            if global_bg_rates[s] == 0: global_bg_rates[s] = minBG
    print 'Background rates calculated (%0.2f min.).'%((time.time()-time0_gbg)/60)
    
    # find candidate hotspots
    print '\n-----------------------CANDIDATE HOTSPOT FINDING-----------------------'
    print 'Finding candidate somatic hotspots in genes...'
    time0_hs = time.time()
    tot_g_with_hs, tot_hotspots = 0, 0
    for g in genes:
        if g.total_nonsilent_mutations < min_clust_vars: continue
        g.find_hotspots(dist=max_hs_dist,min_clust_vars=min_clust_vars)
        if len(g.clusters) > 0:
            tot_g_with_hs += 1
            tot_hotspots += len(g.clusters)
    oflog.writelines('Identified %d candidate hotspots in %d genes for testing.\n'%(tot_hotspots,tot_g_with_hs))
    print 'Candidate hotspots obtained (%0.2f min.).'%((time.time()-time0_hs)/60)

    # Calculate enrichment p-values
    print '\n--------------------------ENRICHMENT ANALYSES--------------------------'
    print 'Calculating gene and hotspot enrichment p-values with negative binomial tests...'
    print '\nPerforming gene enrichment analysis...'
    time0_enr = time.time()
    if use_covars:
        print '  Using clustered covariate background rates.'
        res = [pool.apply_async(get_gene_enrichments_covar,args=(g,ns)) for g in genes if g.total_nonsilent_mutations>0]
        rget = [r.get() for r in res]
        for r in rget:
            index,pvf = r
            genes[index].gene_pval = pvf
        del(res)
    elif use_local:
        print '  Using local background rates.'
        res = [pool.apply_async(get_gene_enrichments_local_bg,args=(g,ns)) for g in genes if g.total_nonsilent_mutations>0]
        rget = [r.get() for r in res]
        for r in rget:
            index,bgprob,bgtype,pvf = r
            genes[index].background_prob = bgprob
            genes[index].enrichment_bg_type = bgtype
            genes[index].gene_pval = pvf
        del(res)
    else:
        print '  Using global background rates.'
        res = [pool.apply_async(get_gene_enrichments_global_bg,args=(g,global_bg_rates,ns)) for g in genes if g.total_nonsilent_mutations>0]
        rget = [r.get() for r in res]
        for r in rget:
            index,bgprob,bgtype,pvf = r
            genes[index].background_prob = bgprob
            genes[index].enrichment_bg_type = bgtype
            genes[index].gene_pval = pvf
        del(res)
    print '\nGene enrichment p-values obtained (%0.2f min.).'%((time.time() - time0_enr)/60)
     
    # Write enrichment output
    print '\nCorrecting p-values and writing enrichment analyses output...'
    ofr = open(outdir+prefix+'gene_enrichments.txt','w') 
    ofr.writelines('\t'.join(['Gene','coordinates','num_nonsilent','num_bg','full_length','coding_length','bg_type','bg_prob',
                              'gene_pval','FDR_BH','num_samples','nonsilent_position_counts','nonsilent_mutation_counts','samples'])+'\n')
    g_pvals = []
    g_enr = [g for g in genes if g.total_nonsilent_mutations>0]
    g_enr.sort(key=lambda x:x.gene_pval)
    for g in g_enr:
        g_pvals.append(g.gene_pval)
    g_fdrs = fdr_BH(array(g_pvals))
    tot_g_output = 0
    for i,g in enumerate(g_enr):
        g.gene_qval = g_fdrs[i]
        bgp = g.background_prob
        if bgp == None: bgp = -1
       
        coords = '%s:%d-%d'%(g.chrom,g.start,g.stop)
        gname,nnon,nbg,tot_len,code_len = g.name,g.total_nonsilent_mutations,g.total_bg_mutations,g.total_length,g.coding_length
        bg_type = g.enrichment_bg_type
        gpv,fdr = g.gene_pval,g.gene_qval
        samples = []
        for s in g.mutations_by_sample:
            if len(g.mutations_by_sample[s]['nonsilent'])>0: samples.append(s)
        nsamps,sampstr = len(samples),';'.join(sorted(samples))
        pos_counter = Counter(g.positions['nonsilent'])
        pos_counts = []
        for pos in sorted(pos_counter.keys()):
            pos_counts.append('%d_%d'%(pos,pos_counter[pos]))
        pos_counts = ';'.join(pos_counts)
        mut_counts = []
        for mut in sorted(g.samples_by_mutations['nonsilent'].keys()):
            count = len(g.samples_by_mutations['nonsilent'][mut])
            mut_counts.append('%s_%d'%(mut,count))
        mut_counts = ';'.join(mut_counts)

        ol = '%s\t%s\t%d\t%d\t%d\t%d\t%s\t%0.3g\t%0.3g\t%0.3g\t%d\t%s\t%s\t%s\n'%(gname,coords,nnon,nbg,tot_len,code_len,bg_type,bgp,
                                                                                  gpv,fdr,nsamps,pos_counts,mut_counts,sampstr) 
        
        ofr.writelines(ol)
        tot_g_output += 1
    ofr.close()
    oflog.writelines('\n%d gene enrichment results reported.\n'%(tot_g_output))

    # get hotspot enrichments
    print '\nCalculating hotspot enrichment p-values with negative binomial tests...'
    time0_hse = time.time()
    if use_covars:
        print '  Using clustered covariate background rates.'
        clust_output = []
        for g in genes:
            if len(g.clusters) > 0:
                oce = get_hotspot_enrichments_covar(g,genes,ns,global_bg_rates,scr,bg_vtype)
                clust_output.append(oce) 
    elif use_local:
        print '  Using local background rates.'
        clust_results = [pool.apply_async(get_hotspot_enrichments_local,args=(g,ns)) for g in genes if len(g.clusters)>0]
        clust_output = [p.get() for p in clust_results]
    else:
        print '  Using global background rates.'
        clust_results = [pool.apply_async(get_hotspot_enrichments_global,args=(g,global_bg_rates,ns)) for g in genes if len(g.clusters)>0]
        clust_output = [p.get() for p in clust_results]
    clust_out_adjust = []
    for co in clust_output: 
        if co == None: continue
        for c in co: clust_out_adjust.append(c)
    cluster_enrichments = sorted(clust_out_adjust,key=lambda x:x[-5])
    # FDR correct p-values
    pvals = []
    for cr in cluster_enrichments: pvals.append(cr[-5])
    fdrs = fdr_BH(array(pvals))
    for i,cr in enumerate(cluster_enrichments): cr.insert(8,fdrs[i])
    print '\nHotspot enrichment p-values obtained (%0.2f min.).'%((time.time() - time0_hse)/60)
     
    # Assign hotspots to their originating genes
    for i,cr in enumerate(cluster_enrichments):
        gindex = gene2index[cr[0]]
        gene = genes[gindex]
        assert gene.index == gindex
        gene.cluster_enrichments.append(cr)
   
    # Write hotspot enrichment output 
    print 'Writing hotspot enrichment analysis output...'
    of = open(outdir+prefix+'hotspot.txt','w')
    of.writelines('\t'.join(['Gene','hotpsot','num_mutations','hotspot_length','effective_length',
                             'bg_type','bg_prob','pval','BH_qval','num_samples','position_counts','mutation_counts','samples'])+'\n')
    for i,cr in enumerate(cluster_enrichments):
        ostr = '%s\t%s\t%d\t%d\t%d\t%s\t%0.3g\t%0.3g\t%0.3g\t%d\t%s\t%s\t%s\n'%(tuple(cr))
        of.writelines(ostr)
    of.close()
    oflog.writelines('%d gene hotspot enrichment results reported.\n'%(len(cluster_enrichments)))

    # Save genes analysis as python pickle object
    cPickle.dump(genes,open(outdir+prefix+'gene_data.pkl','w'))

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
def load_gtf(gtf,genefield,gene_list,genes_to_use):
    '''
    Load and process genes in GTF file.
    '''
    genes = {}
    bad_genes = set() # list for problematic genes

    for line in open(gtf).readlines():
        if line.startswith('#'): continue # some GTFs have header lines
        l = line.strip().split('\t')
        chrom,typ,start,stop,strand,data = l[0],l[2],int(l[3]),int(l[4]),l[6],l[-1]

        if len(chrom.split('_'))>1: continue # skip complex chromosomes (e.g. randoms)
        
        # parse data string
        gene_id,tx_id = '',''
        for d in data.split(';'):
            if d=='': continue
            elif d[0] == ' ': d = d[1:]
            if d.startswith(genefield):
                d1 = d.split(' ')[1]
                gene_id = d1.replace('"','')
        
        # Check if restriction list supplied
        if gene_list != None:
            if gene_id not in genes_to_use: continue

        # get gene data
        if gene_id not in genes:
            genes[gene_id] = {}
            genes[gene_id]['chrom'] = chrom
            genes[gene_id]['strand'] = strand
            genes[gene_id]['exons'] = []
            genes[gene_id]['CDS'] = []
        if chrom != genes[gene_id]['chrom']: bad_genes.add(gene_id)

        reg = (start,stop)
        if typ == 'exon':
            if reg in genes[gene_id]['exons']: continue
            genes[gene_id]['exons'].append(reg)
        elif typ == 'CDS':
            if reg in genes[gene_id]['CDS']: continue
            genes[gene_id]['CDS'].append(reg)
            
    # delete problematic genes
    print '  Deleting %d genes annotated to multiple chromosomes.'%(len(bad_genes))
    for g in bad_genes:
        del(genes[g])

    # return
    for gene_id in genes:
        genes[gene_id]['exons'].sort()
        genes[gene_id]['CDS'].sort()
    return genes

def merge_list(exlist):
    '''
    Merge overlapping intervals from sorted list.
    '''
    merged = []
    overlapper = []
    for i in range(0,len(exlist)-1):
        ex1,ex2 = exlist[i],exlist[i+1]
        if ex2[0] > ex1[1]:
            if len(overlapper) > 0:
                merged.append(overlapper)
                overlapper = [] 
            else: merged.append((ex1[0],ex1[1]))
            if i == len(exlist)-2:
                merged.append((ex2[0],ex2[1]))

        elif ex2[1] >= ex1[0] and ex2[0] <= ex1[1]:
            if len(overlapper) > 0:
                estart = min(overlapper[0],ex1[0],ex2[0])
                estop = max(overlapper[1],ex1[1],ex2[1])
                overlapper = (estart,estop)
            else:
                estart = min(ex1[0],ex2[0])
                estop = max(ex1[1],ex2[1])
                overlapper = (estart,estop)
            if i == len(exlist)-2:
                merged.append(overlapper)
    return merged

def load_nonsilent_terms(tType):
    '''
    Return dictionary for annotation terms for genes/non-silent mutations.
    '''
    if tType == 'illumina':
        terms = ["missense_variant","splice_region_variant","stop_gained","frameshift_variant","inframe_deletion",
                 "splice_acceptor_variant","splice_donor_variant","inframe_insertion",
                 "stop_lost","start_lost","frameshift_variant","protein_altering_variant"]
        return terms
    elif tType == 'annovar':
        terms = ['frameshift_deletion','frameshift_insertion','frameshift_substitution',
                 'nonframeshift_deletion','nonframeshift_insertion','nonframeshift_substitution',
                 'nonsynonymous_SNV','stopgain','stoploss']
        return terms
    elif tType == 'maf':
        terms =  ["Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                 "Missense_Mutation","Nonsense_Mutation","Splice_Site","Translation_Start_Site","Nonstop_Mutation"]
        return terms
    elif os.path.isfile(tType):
        terms = {}
        terms['terms'] = {}
        for line in open(tType).readlines():
            l = line.strip().split('\t')
            if l[0] == 'Gene':
                terms['gene'] = l[1]
            elif l[0] == 'Effect':
                val,term = l[1],l[2]
                if val not in terms['terms']: terms['terms'][val] = [term]
                else: terms['terms'][val].append(term)

        return terms

def count_mutations_from_vcfs(VCFs,names,genes,terms,tType,snps_only,blacklist=None,mapr=None):
    '''
    Read in somatic VCF files and extract mutations in genes.
    '''
                    
    # Get annotation term types
    anno_val,gname_val = '',''
    if tType == 'illumina':
        anno_val,gname_val = 'CSQT','CSQT'
    elif tType == 'annovar':
        anno_val = 'ExonicFunc.refGene'
        gname_val = 'Gene.refGene'

    # count vars
    for j,vcf_f in enumerate(VCFs):
        vcf = VCF(vcf_f,gts012=True,lazy=False)
        name = names[j]
        tot_added = 0
        for i,g in enumerate(genes):
            chrom,start,stop = g.chrom,g.start,g.stop
            if mapr == None:
                gstr = '%s:%d-%d'%(chrom,start,stop)
                g_strings = [gstr]
            else: g_strings = g.mappable_regions
            
            tot_bg,tot_nonsilent = 0,0 
            for gstr in g_strings:
                prior_var = None
                for v in vcf(gstr): 
                    # variant filtering
                    filt = v.FILTER
                    if filt != None: continue # None is PASS with cyvcf2
                    
                    # Get variant info
                    pos = v.POS
                    ref = v.REF.encode('ascii')
                    for a in v.ALT:
                        alt = a.encode('ascii')
                        var_info = '%d_%s_%s'%(pos,ref,alt)
                        # check if duplicate
                        if prior_var != None:
                            if var_info == prior_var: continue
                        prior_var = var_info
                             
                        # check if black-listed site
                        if blacklist != None:
                            if chrom in blacklist:
                                if var_info in blacklist[chrom]: continue
                    
                        # check if nonsilent
                        nonsilent = False
                        if tType in ['illumina','annovar']:
                            try:
                                found_term = False
                                for term in terms: 
                                    if term in v.INFO[anno_val]: found_term = True
                                if tType == 'annovar': 
                                    if 'splicing' in v.INFO['Func.refGene']: found_term = True
                                if found_term:
                                    for msinfo in v.INFO[gname_val].split(','):
                                        if tType == 'illumina':
                                            ms_split = msinfo.split('|')
                                            ginfo,geffect = ms_split[1],ms_split[3]
                                            if g.name == ginfo and geffect in terms: nonsilent = True
                                        else: 
                                            if g.name in msinfo: nonsilent = True
                            except KeyError: pass
                        else:
                            found_term = False
                            anno_vals = terms['terms'].keys()
                            for av in anno_vals:
                                try: anno = v.INFO[av]
                                except KeyError: continue
                                tlist = terms['terms'][av]
                                for t in tlist:
                                    if t in anno: found_term = True
                            if found_term:
                                try:
                                    gname_val = terms['gene']
                                    gdata = v.INFO[gname_val]
                                    if ',' in gdata or '|' in gdata: # checking for common delimiters
                                        if ',' in gdata: 
                                            if g.name in gdata.split(','): nonsilent = True
                                        elif '|' in gdata: 
                                            if g.name in gdata.split('|'): nonsilent = True
                                    else:
                                        if g.name in gdata: nonsilent = True
                                except KeyError: pass

                        if name not in g.mutations_by_sample: 
                            g.mutations_by_sample[name] = {}
                            g.mutations_by_sample[name]['bg'] = 0
                            g.mutations_by_sample[name]['nonsilent'] = []
                    
                        # Check for duplicate lines
                        # strange examples with Strelka where same exact variant is reported >1 time
                        elif name in g.mutations_by_sample:
                            if nonsilent == False:
                                pass
                                #if var_info in g.mutations_by_sample[name]['bg']: continue
                            elif nonsilent == True:
                                if var_info in g.mutations_by_sample[name]['nonsilent']: continue

                        if nonsilent == True:
                            g.mutations_by_sample[name]['nonsilent'].append(var_info)
                            g.positions['nonsilent'].append(pos)
                            if pos not in g.samples_by_positions['nonsilent']:
                                g.samples_by_positions['nonsilent'][pos] = [name]
                            else: g.samples_by_positions['nonsilent'][pos].append(name)
                            if var_info not in g.samples_by_mutations['nonsilent']:
                                g.samples_by_mutations['nonsilent'][var_info] = [name]
                            else: g.samples_by_mutations['nonsilent'][var_info].append(name)
                            tot_nonsilent += 1
                        elif nonsilent == False:
                            g.mutations_by_sample[name]['bg'] += 1 
                            tot_bg += 1
            
            # update mutation counts for gene
            g.total_mutations += tot_nonsilent + tot_bg
            g.total_bg_mutations += tot_bg
            g.total_nonsilent_mutations += tot_nonsilent
   
    return genes

def count_mutations_from_maf(MAF,genes,gene2index,terms,tType,snps_only,blacklist=None):
    '''
    Count gene mutations from MAF input.
    '''
    # Sample names 
    sample_names = set()

    # new genes list
    new_genes = []

    # Loop over lines in MAF file
    for line in open(MAF).readlines()[1:]: # skip header
        l = line.strip().split('\t')
        try:
            gene,chrom,mstart,mstop,varclass,vartype,ref,alt,name = l[0],l[4],int(l[5]),int(l[6]),l[8],l[9],l[10],l[11],l[15]
        except:
            print 'Could not parse line in MAF: %s'%(line)
            sys.exit()
        if gene not in gene2index: continue
        if not chrom.startswith('chr'): chrom = 'chr'+chrom
        pos = mstart # use mutation start position
        var_info = '%d_%s_%s'%(pos,ref,alt)
            
        # Check mutation type
        nonsilent = False
        if varclass in terms:
            if snps_only:
                if vartype != 'SNP': continue
            else: nonsilent = True

        # get gene of interest
        gind = gene2index[gene]
        g = genes[gind]
        gchrom,gstart,gstop = g.chrom,g.start,g.stop
 
        if nonsilent == True:
            if pos < gstart or pos > gstop:
                print ' %s var %s %s not in gene limits'%(gene,chrom,var_info)
                print gchrom,gstart,gstop
                continue
       
        # Add mutation information
        sample_names.add(name)
        if name not in g.mutations_by_sample: 
            g.mutations_by_sample[name] = {}
            g.mutations_by_sample[name]['bg'] = 0
            g.mutations_by_sample[name]['nonsilent'] = []

        if nonsilent == True:
            g.mutations_by_sample[name]['nonsilent'].append(var_info)
            g.positions['nonsilent'].append(pos)
            if pos not in g.samples_by_positions['nonsilent']:
                g.samples_by_positions['nonsilent'][pos] = [name]
            else: g.samples_by_positions['nonsilent'][pos].append(name)
            if var_info not in g.samples_by_mutations['nonsilent']:
                g.samples_by_mutations['nonsilent'][var_info] = [name]
            else: g.samples_by_mutations['nonsilent'][var_info].append(name)
            g.total_nonsilent_mutations += 1
            g.total_mutations += 1

        elif nonsilent == False:
            g.mutations_by_sample[name]['bg'] += 1 
            g.total_bg_mutations += 1
            g.total_mutations += 1
        
    for g in genes:
        new_genes.append(g)
    return new_genes,list(sample_names)

def get_global_bg_rates(genes,names,bg_vtype):
    '''
    Determine per-sample global background rates by counting background mutations within each gene.
    '''
    # Initialize output dictionary
    bg_rates = {}
    for n in names: bg_rates[n] = [0,0] # Initialize with two-element list for counts and total length

    # Get background rates
    for g in genes:
        tlen = g.total_length
        for s in names:
            try:
                if bg_vtype == 'all':
                    nm = g.mutations_by_sample[s]['bg'] + len(g.mutations_by_sample[s]['nonsilent']) 
                elif bg_vtype == 'silent':
                    nm = g.mutations_by_sample[s]['bg']
            except KeyError: nm = 0
            bg_rates[s][0] += nm
            bg_rates[s][1] += tlen

    # Calculate rates and return
    for s in bg_rates:
        rate = bg_rates[s][0] / bg_rates[s][1]
        bg_rates[s] = rate 

    return bg_rates

def get_local_background(genes,VCFs,names,terms,tType,snps_only,blacklist=None,for_covar=False,mapr=None,bg_vtype='all'): 
    # Windows
    wins = [10e3,50e3,100e3,500e3,1e6] #,2e6]
    
    # load mappable regions file if supplied
    if mapr != None:
        mappable = TabixFile(mapr)
        wins = [100e3,500e3,1e6]

    # Get annotation term types
    anno_val,gname_val = '',''
    if tType == 'illumina':
        anno_val,gname_val = 'CSQT','CSQT'
    elif tType == 'annovar':
        anno_val = 'ExonicFunc.refGene'
        gname_val = 'Gene.refGene'

    # Loop over VCFs
    for j,vcf_f in enumerate(VCFs):
        vcf = VCF(vcf_f,gts012=True,lazy=False)
        name = names[j]
        
        use_local_finder = False
        for g in genes:
            if for_covar: use_local_finder = True
            if use_local_finder == False:
                if name not in g.mutations_by_sample: continue
                if g.total_length >= 10e3 and g.total_nonsilent_mutations>0: 
                    g.local_backgrounds[name] = {}
                    g.local_backgrounds[name]['win_size'] = g.total_length
                    if bg_vtype == 'all':
                        num_bg = g.mutations_by_sample[name]['bg'] + len(g.mutations_by_sample[name]['nonsilent'])
                    elif bg_vtype == 'silent':
                        num_bg = g.mutations_by_sample[name]['bg']
                    g.local_backgrounds[name]['bg_rate'] = num_bg / g.total_length
                else: 
                    if g.total_nonsilent_mutations > 0:
                        use_local_finder = True
                    else: continue
            if use_local_finder == True:
                bgs_per_win = []
                wins_use = wins
                if g.total_length > wins_use[-1]: wins_use.append(g.total_length) # case for very long genes
                for win in wins_use:
                    half_win = win/2
                    midp = int(round((g.start+g.stop)/2))
                    start = max(1,midp-half_win)
                    stop = midp + half_win
                    gstr = '%s:%d-%d'%(g.chrom,start,stop)
                    
                    # Get mappable regions in window (if supplied)
                    g_strings = []
                    win_l = 0
                    if mapr == None:
                        g_strings.append(gstr)
                        win_l = (stop-start)+1
                    else:
                        tbq = [t.split('\t') for t in mappable.fetch(gstr)]
                        for t in tbq:
                            qstart,qstop = int(t[1])+1,int(t[2])
                            if qstart < start: ostart = start
                            else: ostart = qstart
                            if qstop > stop: ostop = stop
                            else: ostop = qstop
                            win_l += ostop - ostart + 1
                            MR = '%s:%d-%d'%(g.chrom,ostart,ostop)
                            g_strings.append(MR)
                
                    if win_l == 0:
                        bgs_per_win.append(0)
                        continue
                    
                    mut_count = 0 
                    for gstr in g_strings:
                        prior_var = None
                        for v in vcf(gstr): 
                            # variant filtering
                            filt = v.FILTER
                            if filt != None: continue # None is PASS with cyvcf2 
                
                            # Get variant info
                            for a in v.ALT:
                                var_info = '%d_%s_%s'%(v.POS,v.REF.encode('ascii'),a.encode('ascii'))
                                # check if black-listed site
                                if blacklist != None:
                                    if g.chrom in blacklist:
                                        if var_info in blacklist[g.chrom]: continue
                            
                                # check if duplicate
                                if prior_var != None:
                                    if var_info == prior_var: continue
                                prior_var = var_info
                            
                                if bg_vtype == 'all': mut_count += 1
                                
                                elif bg_vtype == 'silent':
                                    # check if nonsilent
                                    nonsilent = False
                                    if tType in ['illumina','annovar']:
                                        try:
                                            found_term = False
                                            for term in terms: 
                                                if term in v.INFO[anno_val]: found_term = True
                                            if tType == 'annovar': 
                                                if 'splicing' in v.INFO['Func.refGene']: found_term = True
                                            if found_term:
                                                for msinfo in v.INFO[gname_val].split(','):
                                                    if tType == 'illumina':
                                                        ms_split = msinfo.split('|')
                                                        ginfo,geffect = ms_split[1],ms_split[3]
                                                        if g.name == ginfo and geffect in terms: nonsilent = True
                                                    else: 
                                                        if g.name in msinfo: nonsilent = True
                                        except KeyError: pass
                                    else:
                                        found_term = False
                                        anno_vals = terms['terms'].keys()
                                        for av in anno_vals:
                                            try: anno = v.INFO[av]
                                            except KeyError: continue
                                            tlist = terms['terms'][av]
                                            for t in tlist:
                                                if t in anno: found_term = True
                                        if found_term:
                                            try:
                                                gname_val = terms['gene']
                                                gdata = v.INFO[gname_val]
                                                if ',' in gdata or '|' in gdata: # checking for common delimiters
                                                    if ',' in gdata: 
                                                        if g.name in gdata.split(','): nonsilent = True
                                                    elif '|' in gdata: 
                                                        if g.name in gdata.split('|'): nonsilent = True
                                                else:
                                                    if g.name in gdata: nonsilent = True
                                            except KeyError: pass
                                    
                                    if nonsilent == False:
                                        mut_count += 1
       
                    # update bgs list       
                    bgs_per_win.append(mut_count / win_l)

                # Get max rate
                maxi,maxr = 0,bgs_per_win[0]
                for i,bgr in enumerate(bgs_per_win):
                    if bgr > maxr:
                        maxi,maxr = i,bgr
                
                # Set current gene's local background for current VCF sample
                g.local_backgrounds[name] = {}
                g.local_backgrounds[name]['win_size'] = wins_use[maxi]
                g.local_backgrounds[name]['bg_rate'] = maxr
                
    # impute any remaining 0 backgrounds
    for g in genes:
        szero,minbg = [],1
        for s in g.local_backgrounds:
            bg = g.local_backgrounds[s]['bg_rate']
            if bg == 0: szero.append(s)
            else: 
                if bg < minbg: minbg = bg
        for s in szero:
            if minbg == 1: minbg = 1e-8 # floor value in case zeros still seen
            g.local_backgrounds[s]['bg_rate'] = minbg

    return genes

def get_cluster_bg_rates(g,genes,scr,global_bg,bg_vtype):
    '''
    Determine per-gene background rates from clustered genes. 
    If gene is a member of small cluster, the local background rate is used.
    '''
    # Get appropriate background
    bgpf_l = []
    if g.index in scr:
        csf = [s for s in g.local_backgrounds]
        g.enrichment_bg_type = 'local'
        for s in csf: bgpf_l.append(g.local_backgrounds[s]['bg_rate'])
    else:
        csf = [s for s in g.mutations_by_sample]
        g.enrichment_bg_type = 'clustered_regions'
        members = [g.index] + g.covariate_clusters
        gmems = [genes[i] for i in members]
        bgf_rates = get_global_bg_rates(gmems,csf,bg_vtype)
        for s in csf: 
            if bgf_rates[s] == 0: 
                bg_of_s = global_bg[s] # use global rate for zero cases
            else: bg_of_s = bgf_rates[s]
            bgpf_l.append(bg_of_s)
    bgpf = gmean(bgpf_l) # geometric mean
    g.background_prob = bgpf

def get_gene_enrichments_global_bg(g,bg_rates,ns):
    '''
    Function for determining full gene enrichments using global background rates.
    Return tuple data.
    '''
    # calculate p-value for full region
    xf = g.coding_length * ns
    kf = g.total_nonsilent_mutations #- 1
    bg = []
    samples = g.mutations_by_sample.keys()
    for s in samples: bg.append(bg_rates[s])
    bg = gmean(bg) # take geometric mean
    
    # Calculate full region p-value
    try: 
        pvf = betainc(kf,xf-kf+1,bg)
    except:
        print '  error at gene %s with length: %d, num mutations: %d, and bg: %f'%(g.name,g.coding_length,kf,bg)
        pvf = 1
   
    # Return tuple
    return (g.index,bg,'global',pvf)

def get_gene_enrichments_local_bg(g,ns):
    '''
    Function for determining full gene enrichments using local bacgkround rates.
    Return tuple data.
    '''
    # calculate p-value for full region
    xf = g.coding_length * ns
    kf = g.total_nonsilent_mutations #- 1
    bg = []
    samples = g.local_backgrounds.keys()
    for s in samples: bg.append(g.local_backgrounds[s]['bg_rate'])
    bg = gmean(bg) # take geometric mean
    
    # Calculate full region p-value
    try: 
        pvf = betainc(kf,xf-kf+1,bg)
    except:
        print '  error at gene %s with length: %d, num mutations: %d, and bg: %f'%(g.name,g.coding_length,kf,bg)
        pvf = 1

    # Return tuple
    return (g.index,bg,'local',pvf)

def get_gene_enrichments_covar(g,ns):
    '''
    Function for determining gene enrichments.
    Return tuple data.
    '''
    # calculate p-value for full region
    xf = g.coding_length * ns
    kf = g.total_nonsilent_mutations #- 1
    bg = g.background_prob

    # Calculate full region p-value
    try:
        pvf = betainc(kf,xf-kf+1,bg)
    except:
        print '  error at gene %s with length: %d, num mutations: %d, and bg: %f'%(g.name,g.coding_length,kf,bg)
        pvf = 1
    
    # Return tuple
    return (g.index,pvf)

def get_hotspot_enrichments_covar(g,genes,ns,global_bg,scr,bg_vtype):
    '''
    Compute hotspot enrichment p-values using negative binomial test with covariate cluster background rates.
    '''
    # initialize output
    enrich = []
    # Find clusters and compute p-values
    if len(g.clusters) > 0:
        for clust in g.clusters:
            #if g.clusters[clust]['length'] < 2: continue
                
            # Cluster info
            clust_name = g.clusters[clust]['name']

            # data for NB test 
            k = g.clusters[clust]['count']
            c_len = g.clusters[clust]['length']
            x = c_len * ns # cluster length times number of samples
            cs = g.clusters[clust]['samples']
            
            # Get background            
            bgp_l,bg_type = [],''
            if g.index in scr:
                bg_type = 'local'
                for s in cs: bgp_l.append(g.local_backgrounds[s]['bg_rate']) 
            else:
                bg_type = 'clustered_regions'
                members = [g.index] + g.covariate_clusters
                gmems = [genes[i] for i in members]
                bgf_rates = get_global_bg_rates(gmems,cs,bg_vtype)
                for s in cs: 
                    if bgf_rates[s] == 0: bg_of_s = global_bg[s] # use global rate for zeros
                    else: bg_of_s = bgf_rates[s]
                    bgp_l.append(bg_of_s)
            bgp = gmean(bgp_l) # geometric mean
    
            # Calculate p-values
            pv = betainc(k,x-k+1,bgp) 

            # Get counts for positions and mutations
            counter = Counter(g.clusters[clust]['positions'])
            count_s = []
            for pos in sorted(counter.keys()):
                count_s.append('%d_%d'%(pos,counter[pos]))
            count_s = ';'.join(count_s)
            muts = []
            for mut in sorted(g.samples_by_mutations['nonsilent'].keys()):
                pos = int(mut.split('_')[0])
                if pos in counter:
                    for s in g.samples_by_mutations['nonsilent'][mut]:
                        if s in cs: muts.append(mut)
            mut_counter = Counter(muts)
            count_m = []
            for mut in sorted(mut_counter.keys()):
                count_m.append('%s_%d'%(mut,mut_counter[mut]))
            count_m = ';'.join(count_m)

            # Output info
            ol = [g.name,clust_name,k,c_len,x,bg_type,bgp,pv,len(cs),count_s,count_m,';'.join(sorted(list(cs)))]
            enrich.append(ol)
        
        # return
        if len(enrich) > 0:
            return enrich

def get_hotspot_enrichments_local(g,ns):
    '''
    Compute hotspot enrichment p-values using negative binomial test with local background rates.
    '''
    # initialize output
    enrich = []
    bg_type = 'local'
    # Find clusters and compute p-values
    if len(g.clusters) > 0:
        for clust in g.clusters:
            # Cluster info
            clust_name = g.clusters[clust]['name']
            
            # data for NB test 
            k = g.clusters[clust]['count']
            c_len = g.clusters[clust]['length']
            x = c_len * ns # cluster length times number of samples
            cs = g.clusters[clust]['samples']
                    
            bgp_l = []
            for s in cs: bgp_l.append(g.local_backgrounds[s]['bg_rate'])
            bgp = gmean(bgp_l) # geometric mean
            pv = betainc(k,x-k+1,bgp)

            # Get counts for positions and mutations
            counter = Counter(g.clusters[clust]['positions'])
            count_s = []
            for pos in sorted(counter.keys()):
                count_s.append('%d_%d'%(pos,counter[pos]))
            count_s = ';'.join(count_s)
            muts = []
            for mut in sorted(g.samples_by_mutations['nonsilent'].keys()):
                pos = int(mut.split('_')[0])
                if pos in counter:
                    for s in g.samples_by_mutations['nonsilent'][mut]:
                        if s in cs: muts.append(mut)
            mut_counter = Counter(muts)
            count_m = []
            for mut in sorted(mut_counter.keys()):
                count_m.append('%s_%d'%(mut,mut_counter[mut]))
            count_m = ';'.join(count_m)
  
            # Output info
            ol = [g.name,clust_name,k,c_len,x,bg_type,bgp,pv,len(cs),count_s,count_m,';'.join(sorted(list(cs)))]
            enrich.append(ol)
        
        # return
        if len(enrich) > 0:
            return enrich

def get_hotspot_enrichments_global(g,global_bg_rates,ns):
    '''
    Compute hotspot enrichment p-values using negative binomial test with global background rates.
    '''
    # initialize output
    enrich = []
    bg_type = 'global'
    # Find clusters and compute p-values
    if len(g.clusters) > 0:
        for clust in g.clusters: 
            # Cluster info
            clust_name = g.clusters[clust]['name']

            # data for NB test 
            k = g.clusters[clust]['count']
            c_len = g.clusters[clust]['length']
            x = c_len * ns # cluster length times number of samples
            cs = g.clusters[clust]['samples']
            bgp_l = []
            for s in cs: bgp_l.append(global_bg_rates[s])
            bgp = gmean(bgp_l) # geometric mean
            pv = betainc(k,x-k+1,bgp)

            # Get counts for positions and mutations
            counter = Counter(g.clusters[clust]['positions'])
            count_s = []
            for pos in sorted(counter.keys()):
                count_s.append('%d_%d'%(pos,counter[pos]))
            count_s = ';'.join(count_s)
            muts = []
            for mut in sorted(g.samples_by_mutations['nonsilent'].keys()):
                pos = int(mut.split('_')[0])
                if pos in counter:
                    for s in g.samples_by_mutations['nonsilent'][mut]:
                        if s in cs: muts.append(mut)
            mut_counter = Counter(muts)
            count_m = []
            for mut in sorted(mut_counter.keys()):
                count_m.append('%s_%d'%(mut,mut_counter[mut]))
            count_m = ';'.join(count_m)
  
            # Output info
            ol = [g.name,clust_name,k,c_len,x,bg_type,bgp,pv,len(cs),count_s,count_m,';'.join(sorted(list(cs)))]
            enrich.append(ol)
        
        # return
        if len(enrich) > 0:
            return enrich

#####################
# Class definitions #
#####################
class Gene:
    '''
    Gene class.
    '''
    def __init__(self,gene,name,index,mappable,exome_only):
        '''
        Initialize Gene object.
        '''
        self.name = name
        self.index = index

        # set exons
        if len(gene['exons'])==1: self.exons = gene['exons']
        else: self.exons = merge_list(sorted(gene['exons']))
        
        # set coding regions
        if len(gene['CDS']) == 1: self.coding_regions = gene['CDS']
        else: self.coding_regions = merge_list(sorted(gene['CDS']))
        
        # set introns
        self.introns = []
        for i in range(0,len(self.coding_regions)-1):
            e1,e2 = self.coding_regions[i],self.coding_regions[i+1]
            istart,istop = e1[1]+1,e2[0]-1
            self.introns.append((istart,istop))
        
        # other gene info
        self.chrom = gene['chrom']
        self.strand = gene['strand']
        self.start = self.exons[0][0]
        self.stop = self.exons[-1][1]
        
        # Get lengths and mappable regions (if supplied)
        self.mappable_regions = []
        self.coding_length = 0
        self.exonic_length = 0
        self.total_length = 0
        if mappable != None:
            if not exome_only: # don't need full length if using exome only
                tot_len = 0
                tbq = [t.split('\t') for t in mappable.fetch('%s:%d-%d'%(self.chrom,self.start,self.stop))]
                for t in tbq:
                    qstart,qstop = int(t[1])+1,int(t[2])
                    if qstart < self.start: ostart = self.start
                    else: ostart = qstart
                    if qstop > self.stop: ostop = self.stop
                    else: ostop = qstop
                    tot_len += ostop-ostart+1
                    MR = '%s:%d-%d'%(self.chrom,ostart,ostop)
                    self.mappable_regions.append(MR)
                self.total_length += tot_len

            # now for coding length
            for cr in self.coding_regions:
                crstr = '%s:%d-%d'%(self.chrom,cr[0],cr[1])
                crq = [t.split('\t') for t in mappable.fetch(crstr)]
                for t in crq:
                    qstart,qstop = int(t[1])+1,int(t[2])
                    if qstart < cr[0]: ostart = cr[0]
                    else: ostart = qstart
                    if qstop > cr[1]: ostop = cr[1]
                    else: ostop = qstop
                    self.coding_length += ostop-ostart+1
            
            # compute adjusted exonic length if using exome only
            if exome_only:
                for ex in self.exons:
                    exstr = '%s:%d-%d'%(self.chrom,ex[0],ex[1])
                    exq = [t.split('\t') for t in mappable.fetch(exstr)]
                    for t in exq:
                        qstart,qstop = int(t[1])+1,int(t[2])
                        if qstart < ex[0]: ostart = ex[0]
                        else: ostart = qstart
                        if qstop > ex[1]: ostop = ex[1]
                        else: ostop = qstop
                        self.exonic_length += ostop-ostart+1

                # set total length to exonic length if using exome only option
                self.total_length = self.exonic_length

            else:             
                for ex in self.exons:
                    self.exonic_length += ex[1] - ex[0] + 1

        else:
            for cr in self.coding_regions:
                self.coding_length += cr[1] - cr[0] + 1
            for ex in self.exons:
                self.exonic_length += ex[1] - ex[0] + 1

            # Set total gene length
            #   If only considering exome, use exonic length; otherwise, use total length (exons + introns)
            if exome_only:
                self.total_length = self.exonic_length
            else: self.total_length = self.stop - self.start + 1
        
        # mutation info
        self.total_mutations = 0
        self.total_bg_mutations = 0
        self.total_nonsilent_mutations = 0
        self.mutations_by_sample = {}
        self.samples_by_mutations = {}
        #self.samples_by_mutations['bg'] = {}
        self.samples_by_mutations['nonsilent'] = {}
        self.samples_by_positions = {}
        #self.samples_by_positions['bg'] = {}
        self.samples_by_positions['nonsilent'] = {}
        self.positions = {}
        #self.positions['bg'] = []
        self.positions['nonsilent'] = []

        # hotspots
        self.clusters = {}
        self.cluster_enrichments = []

        # covariate info
        self.covariate_clusters = []
        self.local_backgrounds = {}

        # stats
        self.gene_pval = 1
        self.gene_qval = 1
        self.background_prob = None
        self.enrichment_bg_type = None
 
    def find_hotspots(self,dist=50,min_clust_vars=3):
        '''
        Find candidate non-silent mutation hotspots by merging mutations within 'dist' of each other.
        '''
        self.positions['nonsilent'].sort() # Sort positions
        pos_to_test = self.positions['nonsilent'] # considering all non-silent mutations
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
                for s in self.samples_by_positions['nonsilent'][curr]: c_samples.add(s)

            # Test if neighboring position is within distance
            # If true, add neighbor to current cluster
            if (next_p - curr) <= dist: 
                cluster.append(next_p)
                counts += 1
                for s in self.samples_by_positions['nonsilent'][next_p]: c_samples.add(s)
            
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
                    for s in self.samples_by_positions['nonsilent'][next_p]: c_samples.add(s)
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
                        for s in self.samples_by_positions['nonsilent'][next_p]: c_samples.add(s)
                        self.clusters[clust_num]['samples'] = c_samples
                    
                # Eliminate if not of sufficient size
                else:
                    cluster = [next_p]
                    c_samples = set()
                    for s in self.samples_by_positions['nonsilent'][next_p]: c_samples.add(s)
                    counts = 1
    

##################
# Math functions #
##################
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


