from __future__ import division, print_function
import sys,os
import numpy as np
import argparse
import subprocess
import multiprocessing as mp
from pysam import FastaFile, TabixFile
from math import ceil 

########
# MAIN #
########
def main():
    
    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################
    
    description = 'Create region covariates file for regions in BED file from sequence features and external data'
    usage = 'python %(prog)s <regions.bed> <genome.fa> [options]'
    parser = argparse.ArgumentParser(usage=usage,description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('regions',metavar='regions.bed',type=str,
                        help='Input BED file of regions (Required).')
    parser.add_argument('genome',metavar='genome',type=str,
                        help='Indexed genome fasta file (Required).')

    parser.add_argument('-o','--outname',type=str,default='./region_covariates.txt',dest='outname',help='Provide output filename.')
    parser.add_argument('-t','--table-files',default=[],action='append',dest='tables',
                        help='Provide one or more additional table files with covariate info for regions. Use as -t file1 -t file2 ... -t fileN. \
                        Each table file must be tab-delimited with one column header. First column is reserved for region string, remaining are \
                        for information. Names for each covariate are read from corresponding header line. Provided values are not adjusted and \
                        regions in BED with no information for one or more covariates are set to "NA".')
    parser.add_argument('-i','--interval-files',default=None,dest='interval_fns',
                        help='Provide a text file listing paths to bgzip compressed and tabix-indexed interval files in BedGraph format. \
                        This file should be tab-delimited with no header in the format: file_path interval_dataset_name. This program will \
                        extract intervals from these files for each feature (region) and report the average value (in cases of multiple overlaps)\
                        from these intervals. If no immediate overlapping intervals are found, a wider search window is scanned to find a \
                        proximal interval; if one is not found after this broader search, "None" is reported.\
                        NOTE - numerical score data (e.g. 4th column of BedGraph files) is expected here.')
    parser.add_argument('-p','--processors',type=int,default=1,dest='nprocessors',help='Set number of processors for parallel runs.')
 
    parser.add_argument('--repliseq-fns',type=str,default=None,dest='repliseq_fns',
                        help='DEPRECATED - use -i/--interval-files option instead.\
                        Provide a file with paths to RepliSeq data (for replication timing information). This file should be tab-delimited \
                        with no header in the format: file_path sampleID. The file paths should point to bed/bedgraph files compressed with \
                        bgzip and indexed with tabix. This program will extract a replication timing value from each file for each region by \
                        scanning overlapping intervals in the files. For multiple intersecting intervals, the average value is taken. If no data \
                        exists in the immediate region vicinity, wider windows around each are scanned until a value is determined; otherwise, \
                        "None" is reported.')
    

    # PARSE INPUTS #
    args = parser.parse_args()
    vargs = vars(args)
    regions_fn = vargs['regions']
    if not os.path.isfile(regions_fn): parser.error('BED file does not exist!')
    genome = vargs['genome']
    if not os.path.isfile(genome): parser.error('Genome file not found!')
    
    # Parse options
    ofn = vargs['outname']
    tables = vargs['tables']
    interval_fns = vargs['interval_fns']
    nprocess = vargs['nprocessors']
    repliseq_fns = vargs['repliseq_fns']
    if repliseq_fns != None:
        if interval_fns == None:
            print('--repliseq-fns option is deprecated, setting -i option to this input.')
            interval_fns = repliseq_fns
        else:
            parser.error('Using deprecated --repliseq-fns option with -i/--interval-files option is not permitted. Exiting.') 

    # check tables
    if len(tables) > 0:
        for t in tables:
            if not os.path.isfile(t): parser.error('Table file not found: %s'%(t))

    # get interval files
    int_fns = []
    if interval_fns != None:
        for x in open(interval_fns).readlines():
            x = x.strip('\n').split('\t')
            rsfn,rsn = x[0],x[1]
            if not os.path.isfile(rsfn): parser.error('File for %s not found!'%(rsn))
            int_fns.append([x[0],x[1]])
       
    # set up multiprocessing pool
    pool = mp.Pool(processes=nprocess) 

    # Load in regions
    print('Loading regions...')
    regions = load_regions(regions_fn)
    print('Loaded %d regions from input BED file.'%(len(regions)))

    # load table data
    if len(tables) > 0:
        EXP = load_tables(tables)

    # Get sequence data
    GENOME = FastaFile(genome)
    for r in regions:
        r.get_sequence(GENOME)
        r.get_seq_gc_cont()
        r.get_seq_CpG()
        if len(tables) > 0:
            samples = EXP['samples']
            if r.region_string in EXP['values']:
                r.get_exp(samples,EXP['values'][r.region_string])
            else: 
                for s in samples: r.extra[s] = 'NA'
   
    # If replication timing data provided
    if len(int_fns) > 0:
        dones,was_done = [],0
        chunk_size = 1000
        num_chunks = int(ceil(len(regions)/chunk_size))
        chunks = []
        for i in range(0,num_chunks):
            rr = range(i*chunk_size,min((i+1)*chunk_size,len(regions)))
            chunks.append(regions[rr[0]:rr[-1]+1])
        print('  Divided %d regions into %d region chunks.'%(len(regions),num_chunks))
        res = [pool.apply_async(get_interval_data,args=(c,int_fns),callback=dones.append) for c in chunks]
        while len(dones) != num_chunks:
            if len(dones)%10==0 and was_done != len(dones):
                was_done = len(dones)
                print('    %d of %d region chunks complete.'%(len(dones),num_chunks))
        regions = []
        for r in res:
            rget = r.get()
            for rr in rget: regions.append(rr)
    
    # write output
    of = open(ofn,'w')
    header = ['Region','length','GC','CpG']
    if len(tables) > 0: 
        for s in EXP['samples']: header.append('%s'%(s)) 
    for ints in int_fns: header.append(ints[1])
    of.writelines('\t'.join(header)+'\n')
    for r in regions:
        ostr = [r.region_string,str(r.length),'%0.3f'%(r.GC),'%0.3f'%(r.CpG)]
        if len(tables)>0:
            for s in EXP['samples']: ostr.append('%s'%(r.extra[s]))
        for ints in int_fns:
            samp = ints[1]
            if r.intervalData[samp] != None:
                ostr.append('%0.3f'%(r.intervalData[samp]))
            else: ostr.append(str(r.intervalData[samp]))
        ostr = '\t'.join(ostr)+'\n'
        of.writelines(ostr)
    of.close()

#####################
# UTILITY FUNCTIONS #
#####################
def load_regions(regions_fn):
    '''
    Load in regions BED file. 
    Returns list of Region class variables.
    '''
    regions = []
    for i,r in enumerate(open(regions_fn).readlines()): # testing
        region = Region(r,i)
        regions.append(region)
    return regions

def load_tables(tables):
    ''' 
    Load covariate data from tables.
    
    Expect one header line, first column with gene name and remaining columns with sample IDs.
    '''
 
    EXP = {}
    for table in tables:
        lines = open(table).readlines()
        samples = lines[0].strip().split('\t')[1:]
        if 'samples' not in EXP: EXP['samples'] = samples
        else: EXP['samples'] += samples
        
        if 'values' not in EXP: EXP['values'] = {}
        for line in lines[1:]:
            l = line.strip().split('\t')
            g,vals = l[0],l[1:]
            if g not in EXP['values']: EXP['values'][g] = []
            for v in vals:
                EXP['values'][g].append(v)
   
    return EXP

def get_interval_data(regions,INT):
    '''
    Get interval data for region for input files in int_fns. 
    Computes mean data value in 100 kb window around region midpoint.
    '''
    for fn,name in INT:
        tb = TabixFile(fn)

        for r in regions:
        
            # Get region for searching replication timing data
            r_len = r.length
            midp = round((r.start+r.stop)/2)
            min_width = 10e3 # search region at least 1 kb
            if r_len < min_width:
                start = midp - round(min_width/2)
                stop = midp + round(min_width/2)
                rstr = '%s:%d-%d'%(r.chrom,start,stop)
            else: rstr = r.region_string

            it_regions = tb.fetch(rstr)
            intData = []
            for rtr in it_regions:
                if rtr == '': continue
                intData.append(float(rtr.split('\t')[-1]))

            if len(intData)>0:
                r.intervalData[name] = np.mean(intData)
                continue
            else: 
                # Extend search if value not found
                extends = [50e3,100e3,500e3,1e6]
                found = False
                for e in extends:
                    start = max(1,midp - round(e/2))
                    stop = midp + round(e/2)
                    rstr = '%s:%d-%d'%(r.chrom,start,stop)
                                       
                    it_regions = tb.fetch(rstr)
                    for rtr in it_regions:
                        if rtr == '': continue
                        intData.append(float(rtr.split('\t')[-1]))
                        found = True
                    if found == True:
                        r.intervalData[name] = np.mean(intData)
                        break
                
                if found == False:
                    r.intervalData[name] = None

    return regions

#####################
# Class definitions #
#####################
class Region:
    '''
    Region class.
    '''
    def __init__(self,region,index):
        '''
        Initialize region. Parse BED file line for region.
        '''
        # Data from BED file
        self.index = index
        r = region.strip('\n').split('\t')
        self.chrom = r[0]
        self.start = int(r[1])+1 # Add one to change 0-indexing in BED
        self.stop = int(r[2]) 
        self.length = self.stop-self.start+1
        self.region_string = '%s:%d-%d'%(self.chrom,self.start,self.stop)
        self.name = r[3]
    
        # set up other values
        self.seq = None
        self.GC = None
        self.CpG = None
        self.extra = {}
        self.intervalData = {}

    def get_sequence(self,genome):
        ''' '''
        rstr = self.region_string 
        seq = genome.fetch(region=rstr)
        self.seq = seq

    def get_seq_gc_cont(self):
        ''' Determine GC content of sequence. ''' 
    
        seq = self.seq
        
        numC = len([x for x in seq if x == 'C'])
        numG = len([x for x in seq if x == 'G'])
        tot = len(seq)

        GCcont = (numC+numG) / tot
        self.GC = GCcont

    def get_seq_CpG(self):
        ''' Determine normalized CpG content of sequence. '''
        seq = self.seq
        tot = len(seq)
    
        # Get C and G counts
        numC = len([x for x in seq if x == 'C'])
        numG = len([x for x in seq if x == 'G'])    

        # Count CpGs
        CpG = 0
        for (a,b) in zip(seq[0:-1],seq[1:]):
            if a == 'C' and b == 'G': CpG += 1
        
        # Compute normalized CpG content
        # norm CpG = (# CpGs) / [(numC * numG) / (length of seq)]
        Exp_CpG = numC * numG
        CpG = CpG*tot
        try: CpG_cont = CpG/Exp_CpG
        except ZeroDivisionError: CpG_cont = 0.0
        self.CpG = CpG_cont

    def get_exp(self,samples,vals):
        ''' Assign expression.'''
        for i,s in enumerate(samples): self.extra[s] = vals[i]


############
# Main call #
#############
if __name__ == '__main__': main()
