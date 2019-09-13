from __future__ import division, print_function
import sys,os
import math
import numpy as np
from optparse import OptionParser
import subprocess
import multiprocessing as mp
import glob
from math_funcs import run_ap_fast, run_ap_slow

'''
Gene covariate clustering code, part of MutEnricher.
Created by Anthony R. Soltis (anthony.soltis.ctr@usuhs.edu)
'''

########
# MAIN #
########
def covariate_cluster(cfn,weights_fn,sim_prefix,ap_iters,convits,alg,pool,contigs=None):
    '''
    Main driver function for clustering analysis.
    '''
    # Read in covariate information
    print('Loading covariates information...')
    if contigs != None:
        print('  performing by contig analysis.')
        covars,vals_D = read_covariates_file(cfn,contigs)
    else:
        print('  performing analysis using all genes.')
        covars,vals_D = read_covariates_file(cfn)
    print('\nCovariates loaded from input file: ',covars)

    # Load weights
    weights = []
    if weights_fn == None:
        weights = np.array(len(covars)*[1/len(covars)])
    else:
        wd = {}
        for l in open(weights_fn).readlines():
            l = l.strip().split('\t')
            cv,w = l[0],float(l[1])
            wd[cv] = w
        for c in covars:
            weights.append(wd[c])
        weights = np.array(weights)/sum(weights) # assure weights sum to 1
        print('\nUsing user-supplied weights (sum = 1): ',weights)

    # Compute similarities by contig
    print('\nComputing similarities...')
    for contig in vals_D:
        # set up output file
        sfn = sim_prefix + '%s/similarities.txt'%(contig)
        os.system('mkdir -p %s/%s'%(sim_prefix,contig))

        # Get and z-score covariates matrix
        vals = np.array(vals_D[contig]['values']) 
        naninds = np.argwhere(np.isnan(vals))
        nanmeans = np.nanmean(vals,axis=0)
        for ni in naninds:
            i,j = ni
            vals[i,j] = nanmeans[j] # set NaN values to mean of column
        vals = vals - np.mean(vals,axis=0)
        vstd = np.std(vals,axis=0)
        vstd[vstd==0] = 1 # protect against divide by zero
        vals = vals / vstd

        vals_D[contig]['zvals'] = vals
        vals_D[contig]['weights'] = weights
        vals_D[contig]['sim_fn'] = sfn
    
    # call to compute
    res = [pool.apply_async(compute_similarities,args=(vals_D[c]['zvals'],vals_D[c]['sim_fn'],vals_D[c]['weights'])) for c in vals_D]
    resg = [r.get() for r in res]
    print('Similarities computed.')

    # Run AP cluster    
    print('\nRunning affinity propagation on similarities...')
    res = [pool.apply_async(apcluster,args=(vals_D[c]['sim_fn'],sim_prefix+c+'/',None,ap_iters,convits,alg)) for c in vals_D]
    resc = [r.get() for r in res]
    # CODE TO CHECK IF ANY AP RUNS FAILED #
    converged = False
    potentials = np.array([-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1])
    while converged == False:
        unfinished = []
        for c in vals_D:
            if not os.path.isfile(sim_prefix+c+'/idx.txt'):
                unfinished.append(c)
        if len(unfinished) == 0: converged = True
        else: 
            print('  re-running unfinished contigs:',unfinished)
            svals = []
            for c in unfinished:
                s = [float(x.split()[1]) for x in open(sim_prefix+c+'/summary.txt').readlines() if x.startswith('Preferences:')][0]
                valid_potentials = np.array([x for x in potentials if x != s]) # avoid re-using already used value!
                diffs = np.sqrt(pow(valid_potentials - s,2))
                s_new = np.array(valid_potentials[diffs==min(diffs)][0]) # take first case if ties
                svals.append(s_new)    
            resRe = [pool.apply_async(apcluster,args=(vals_D[c]['sim_fn'],sim_prefix+c+'/',svals[i],ap_iters,convits,alg)) for i,c in enumerate(unfinished)]
            resRec = [r.get() for r in resRe] 
    print('Affinity propagation clusters generated.')
    
    # Clean up directories by removing similarities files
    for c in vals_D:
        simfn = vals_D[c]['sim_fn']
        os.system('rm %s'%(simfn))

    print('\nWriting cluster outputs...')
    clusters_D = {}
    for c in vals_D:
        od = sim_prefix+c+'/'
        cl = write_clusters(od,vals_D[c]['ids'])
        clusters_D[c] = cl
    return clusters_D
    print('\nGene clustering finished.')

#####################
# UTILITY FUNCTIONS #
#####################
def read_covariates_file(cfn,contigs=None):
    '''
    Read in covariates file and return matrix of covariates.
    '''
    # Read in lines
    lines = open(cfn).readlines()

    # Read in header line to get covariate names
    covars = [x for x in lines[0].strip().split('\t')[1:]]

    # get covariate values by contig and return as dictionaries 
    VALUES = {}
    for line in lines[1:]:
        l = line.strip().split('\t')
        gene = l[0]
        if contigs != None:
            try: contig = contigs[gene]
            except KeyError: continue
        else: contig = 'ALL'
        if contig not in VALUES:
            VALUES[contig] = {}
            VALUES[contig]['ids'] = [gene]
            VALUES[contig]['values'] = []
        else:
            VALUES[contig]['ids'].append(gene)
        vals = []
        for v in l[1:]:
            try:
                vals.append(float(v))
            except:
                vals.append(np.nan)
        VALUES[contig]['values'].append(vals)

    return covars,VALUES

def compute_similarities(mat,sim_fn,weights):
    '''
    Compute pair-wise similarities using negative Euclidean distance.
    '''
    of = open(sim_fn,'w')
    nr = mat.shape[0]
    for i in range(0,nr-1):
        for j in range(i+1,nr):
            v1 = mat[i,:]
            v2 = mat[j,:]
            sim = dist(v1,v2,weights)
            of.writelines('%d %d %0.3g\n'%(i+1,j+1,sim)) # only need to write s(i->j) in this case
    of.close() 
    print('  %s done.'%(sim_fn.split('/')[-2]))

def dist(v1,v2,weights,metric='negEuclidean'):
    '''
    Function for calculating various distances between two vectors.
    
    Types:
        'negEuclidean' - negative of Euclidean distance
    '''
    if metric == 'negEuclidean':
        diff = (v1-v2)*weights
        d = -(sum(pow(diff,2)))

    return d

def apcluster(sim_fn,outdir,s=None,ap_iters=1000,convits=50,alg='fast'):
    '''
    Call apcluster code
    '''
    if alg == 'fast':
        converged = run_ap_fast(sim_fn,preference=s,maxits=ap_iters,convits=convits,outdir=outdir)
    elif alg == 'slow':
        converged = run_ap_slow(sim_fn,preference=s,maxits=ap_iters,convits=convits,outdir=outdir)
    print('  %s done.'%(sim_fn.split('/')[-2]))
    return converged

def write_clusters(outdir,ID):

    # Get clusters
    cluster_num = 0
    clusters = {}
    for i,idx in enumerate(open(outdir+'idx.txt').readlines()):
        idx = int(idx.strip('\n'))
        exemplar = ID[idx-1]
        member = ID[i]
        if exemplar not in clusters:
            cluster_num += 1
            clusters[exemplar] = {}
            clusters[exemplar]['number'] = cluster_num
            clusters[exemplar]['members'] = [member]
        else:
            clusters[exemplar]['members'].append(member)
    clusters_list = []
    for ex in clusters:
        num = clusters[ex]['number']
        mem = clusters[ex]['members']
        clusters_list.append([num,ex,mem])
    clusters_list.sort(key=lambda x: x[0])
    
    # Write clusters
    ofap = open('%s'%(outdir+'clusters.txt'),'w')
    for c in clusters_list:
        ofap.writelines('%d\t%s\t%s\n'%(c[0],c[1],';'.join(c[2])))
    ofap.close()

    return clusters_list

