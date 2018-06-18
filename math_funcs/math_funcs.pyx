from __future__ import division
from math import exp
from collections import Counter
import numpy as np
import fileinput
import time

##################
# Math functions #
##################

def WAP(x,tau):
    wap = 0
    xC = Counter(x)
    x = list(set(sorted(x)))

    for i in range(0,len(x)-1):
        for j in range(i+1,len(x)):
            wap += xC[x[i]] * xC[x[j]] * exp(-pow(x[i]-x[j],2) / (2*pow(tau,2))) 
    
    return wap

def run_ap_slow(sim_fn,preference=None,convits=200,maxits=1000,dampfact=0.9,outdir='./'):
    '''
    Slow version of affinity propagation code.

    Each line in sim_fn should be space-delimited and of form: i j s, where i and j are indices and s is the similarity value.

    NOTE: this implementation only needs the lower triangle of the similarity matrix. It assumes and sets S(i,j) = S(j,i).

    Affinity propagation citation:
    Frey, B.J. and Dueck, D., "Clustering by Passing Messages Between Data Points", Science Feb. 2007
    '''

    time0 = time.time()

    # load input file as dictionary 
    Sd = {}
    inds = set()
    for line in fileinput.input(sim_fn):
        l = line.strip().split()
        i,j,s = int(l[0])-1,int(l[1])-1,float(l[2])
        if i not in Sd: Sd[i] = {}
        Sd[i][j] = s
        inds.add(i)
        inds.add(j)

    # create matrix for similarities
    inds = sorted(list(inds))
    nsamps = len(inds)
    S = np.zeros((nsamps,nsamps))
    for ii in range(0,len(inds)-1):
        for jj in range(ii+1,len(inds)):
            i = inds[ii]
            j = inds[jj]
            sim = Sd[i][j]
            S[i,j] = sim
            S[j,i] = sim # Set S(i,j) = S(j,i)
    del(Sd)
    del(inds)

    # Set preference
    if preference == None:
        preference = np.median(S[~np.identity(len(S)).astype(bool)])
        S.flat[::(nsamps+1)] = preference
    else:
        S.flat[::(nsamps+1)] = preference
    
    # Degeneracy check
    random_state = np.random.RandomState(0)
    S += ((np.finfo(np.double).eps * S + np.finfo(np.double).tiny * 100) * random_state.randn(nsamps,nsamps))
    realmax = np.finfo(np.double).max
    S[np.isneginf(S)] = -realmax

    # Set up iteration variables
    A = np.zeros((nsamps,nsamps))
    R = np.zeros((nsamps,nsamps))
    e = np.zeros((nsamps,convits))
    ind = np.arange(nsamps)

    # outputs
    converged = False

    for iteration in range(maxits):
       
       # compute responsibilities
        for ii in range(nsamps):
            oldR = R[ii,:].copy()
            AS = A[ii,:] + S[ii,:]
            I = np.argmax(AS)
            Y = AS[I]
            AS[I] = -np.inf
            Y2 = np.max(AS)
            R[ii,:] = S[ii,:]-Y
            R[ii,I] = S[ii,I]-Y2
            R[ii,:] = (1-dampfact)*R[ii,:] + dampfact*oldR
        
        # compute availabilities
        for jj in range(nsamps):
            oldA = A[:,jj].copy()
            Rp = np.zeros(R[:,jj].size)
            np.maximum(R[:,jj],0,Rp)
            Rp[jj] = R[jj,jj]
            A[:,jj] = np.sum(Rp)-Rp
            dA = A[jj,jj].copy()
            A[:,jj] = np.minimum(A[:,jj],0,A[:,jj])
            A[jj,jj] = dA
            A[:,jj] = (1-dampfact)*A[:,jj] + dampfact*oldA
        
        # check convergence
        E = (np.diag(A) + np.diag(R)) > 0
        e[:, iteration % convits] = E
        K = np.sum(E,axis=0)
           
        if iteration >= convits:
            se = np.sum(e,axis=1)
            unconverged = (np.sum((se==convits)+(se==0)) != nsamps)
            if (not unconverged and (K>0)) or (iteration == maxits):
                break

    totIters = iteration + 1
    if not unconverged: converged = True

    # Initialize output summary.txt file
    sout = open(outdir+'summary.txt','w')
    sout.writelines('maxits=%d\n'%(maxits))
    sout.writelines('convits=%d\n'%(convits))
    sout.writelines('dampfact=%s\n'%(str(dampfact)))
    sout.writelines('number of data points: %d\n'%(nsamps))
    sout.writelines('Preferences: %f\n\n'%(preference))

    if converged:
        I = np.where(np.diag(A+R)>0)[0]
        K = I.size
        if K > 0:
            c = np.argmax(S[:,I],axis=1)
            c[I] = np.arange(K)
            for k in range(K):
                ii = np.where(c==k)[0]
                j = np.argmax(np.sum(S[ii[:,np.newaxis],ii],axis=0))
                I[k] = ii[j]
            c = np.argmax(S[:,I],axis=1)
            c[I] = np.arange(K)
            labels = I[c]
            c_inds = np.unique(labels)
            labels = np.searchsorted(c_inds,labels)
            idx = np.array([c_inds[i] for i in labels])+1
         
        else:
            labels = np.empty((nsamps,1))
            c_inds = []
            labels.fill(np.nan)
            idx = np.empty((nsamps,1))
            idx.fill(np.nan)

        # write out idx.txt file
        idxout = open(outdir+'idx.txt','w')
        for i in idx: idxout.writelines('%d\n'%(i))
        idxout.close()

        # add to summary
        sout.writelines('Number of identified clusters: %d\n'%(len(c_inds)))
        sout.writelines('Number of iterations: %d\n'%(totIters))
        sout.writelines('Elapsed time: %0.2f minutes\n'%((time.time() - time0)/60))
        sout.close()

    return converged

def run_ap_fast(sim_fn,preference=None,convits=200,maxits=1000,dampfact=0.9,outdir='./'):
    '''
    Fast version of affinity propagation run code.

    Each line in sim_fn should be space-delimited and of form: i j s, where i and j are indices and s is the similarity value.

    NOTE: this implementation only needs the lower triangle of the similarity matrix. It assumes and sets S(i,j) = S(j,i).

    Affinity propagation citation:
    Frey, B.J. and Dueck, D., "Clustering by Passing Messages Between Data Points", Science Feb. 2007
    '''

    time0 = time.time()

    # load input file as dictionary 
    Sd = {}
    inds = set()
    for line in fileinput.input(sim_fn):
        l = line.strip().split()
        i,j,s = int(l[0])-1,int(l[1])-1,float(l[2])
        if i not in Sd: Sd[i] = {}
        Sd[i][j] = s
        inds.add(i)
        inds.add(j)

    # create matrix for similarities
    inds = sorted(list(inds))
    nsamps = len(inds)
    S = np.zeros((nsamps,nsamps))
    for ii in range(0,len(inds)-1):
        for jj in range(ii+1,len(inds)):
            i = inds[ii]
            j = inds[jj]
            sim = Sd[i][j]
            S[i,j] = sim
            S[j,i] = sim # Set S(i,j) = S(j,i)
    del(Sd)
    del(inds)

    # Set preference
    if preference == None:
        preference = np.median(S[~np.identity(len(S)).astype(bool)])
        S.flat[::(nsamps+1)] = preference
    else:
        S.flat[::(nsamps+1)] = preference
    
    # Degeneracy check
    random_state = np.random.RandomState(0)
    S += ((np.finfo(np.double).eps * S + np.finfo(np.double).tiny * 100) * random_state.randn(nsamps,nsamps))
    realmax = np.finfo(np.double).max
    S[np.isneginf(S)] = -realmax

    # Set up iteration variables
    A = np.zeros((nsamps,nsamps))
    R = np.zeros((nsamps,nsamps))
    e = np.zeros((nsamps,convits))
    old = np.zeros((nsamps,nsamps))
    ind = np.arange(nsamps)

    # outputs
    converged = False

    for iteration in range(maxits):
        
        # compute responsibilities
        np.add(A,S,old)
        I = np.argmax(old,axis=1)
        Y = old[ind,I]
        old[ind,I] = -np.inf
        Y2 = np.max(old,axis=1)
        np.subtract(S,Y[:,None],old)
        old[ind,I] = S[ind, I] - Y2
        R = R*dampfact + old*(1-dampfact)
       
        # compute availabilities
        np.maximum(R,0,old)
        old.flat[::nsamps + 1] = R.flat[::nsamps + 1]
        old -= np.sum(old,axis=0)
        dA = np.diag(old).copy()
        old.clip(0, np.inf, old)
        old.flat[::nsamps+1] = dA
        A = A*dampfact - old*(1-dampfact)
        
        # check convergence
        E = (np.diag(A) + np.diag(R)) > 0
        e[:, iteration % convits] = E
        K = np.sum(E,axis=0)
           
        if iteration >= convits:
            se = np.sum(e,axis=1)
            unconverged = (np.sum((se==convits)+(se==0)) != nsamps)
            if (not unconverged and (K>0)) or (iteration == maxits):
                break

    totIters = iteration + 1
    if not unconverged: converged = True

    # Initialize output summary.txt file
    sout = open(outdir+'summary.txt','w')
    sout.writelines('maxits=%d\n'%(maxits))
    sout.writelines('convits=%d\n'%(convits))
    sout.writelines('dampfact=%s\n'%(str(dampfact)))
    sout.writelines('number of data points: %d\n'%(nsamps))
    sout.writelines('Preferences: %f\n\n'%(preference))

    if converged:
        I = np.where(np.diag(A+R)>0)[0]
        K = I.size
        if K > 0:
            c = np.argmax(S[:,I],axis=1)
            c[I] = np.arange(K)
            for k in range(K):
                ii = np.where(c==k)[0]
                j = np.argmax(np.sum(S[ii[:,np.newaxis],ii],axis=0))
                I[k] = ii[j]
            c = np.argmax(S[:,I],axis=1)
            c[I] = np.arange(K)
            labels = I[c]
            c_inds = np.unique(labels)
            labels = np.searchsorted(c_inds,labels)
            idx = np.array([c_inds[i] for i in labels])+1
         
        else:
            labels = np.empty((nsamps,1))
            c_inds = []
            labels.fill(np.nan)
            idx = np.empty((nsamps,1))
            idx.fill(np.nan)

        # write out idx.txt file
        idxout = open(outdir+'idx.txt','w')
        for i in idx: idxout.writelines('%d\n'%(i))
        idxout.close()

        # add to summary
        sout.writelines('Number of identified clusters: %d\n'%(len(c_inds)))
        sout.writelines('Number of iterations: %d\n'%(totIters))
        sout.writelines('Elapsed time: %0.2f minutes\n'%((time.time() - time0)/60))
        sout.close()

    return converged

