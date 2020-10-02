from __future__ import division, print_function
import sys,os
import numpy as np
import argparse
import multiprocessing as mp
from pysam import FastaFile, TabixFile
from math import ceil
import gzip

########
# MAIN #
########
def main():
    
    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################
    
    description = 'Create gene covariates file for genes in GTF from sequence features and external data'
    usage = 'python %(prog)s <GTF> <genome.fa> [options]'
    parser = argparse.ArgumentParser(usage=usage,description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('GTF',metavar='genes.gtf',type=str,
                        help='Input GTF file (Required).')
    parser.add_argument('genome',metavar='genome',type=str,
                        help='Indexed genome fasta file (Required).')

    parser.add_argument('--gene-field',type=str,default='gene_id',dest='genefield',
                        help='Provide field name from input GTF containing gene name/id information.')
    parser.add_argument('-o','--outname',type=str,default='./gene_covariates.txt',dest='outname',help='Provide output filename.')
    parser.add_argument('-t','--table-files',default=[],action='append',dest='tables',
                        help='Provide one or more additional table files with covariate info for genes. Use as -t file1 -t file2 ... -t fileN. \
                        Each table file must be tab-delimited with one column header. First column is reserved for gene name, remaining are \
                        for information. Names for each covariate are read from corresponding header line. Provided values are not adjusted and \
                        genes in GTF with no information for one or more covariates are set to "NA".')
    parser.add_argument('-i','--interval-files',default=None,dest='interval_fns',
                        help='Provide a text file listing paths to bgzip compressed and tabix-indexed interval files in BedGraph format. \
                        This file should be tab-delimited with no header in the format: file_path interval_dataset_name. This program will \
                        extract intervals from these files for each feature (gene) and report the average value (in cases of multiple overlaps)\
                        from these intervals. If no immediate overlapping intervals are found, a wider search window is scanned to find a \
                        proximal interval; if one is not found after this broader search, "None" is reported.\
                        NOTE - numerical score data (e.g. 4th column of BedGraph files) is expected here.')
    parser.add_argument('-p','--processors',type=int,default=1,dest='nprocessors',help='Set number of processors for parallel runs.')
    parser.add_argument('-g','--gene-list',type=str,default=None,dest='gene_list',
                        help='Provide list of genes to which analysis should be restricted (one gene per-line in text file). \
                        Analysis will only considers genes from GTF file that are present in this list. \
                        Default behavior is to query all coding genes present in input GTF.')

    parser.add_argument('--repliseq-fns',type=str,default=None,dest='repliseq_fns',
                        help='DEPRECATED - use -i/--interval-files option instead.\
                        Provide a file with paths to RepliSeq data (for replication timing information). This file should be tab-delimited \
                        with no header in the format: file_path sampleID. The file paths should point to bed/bedgraph files compressed with \
                        bgzip and indexed with tabix. This program will extract a replication timing value from each file for each gene by \
                        scanning overlapping intervals in the files. For multiple intersecting intervals, the average value is taken. If no data \
                        exists in the immediate gene vicinity, wider windows around the gene are scanned until a value is determined; otherwise, \
                        "None" is reported.')
    

    # PARSE INPUTS #
    args = parser.parse_args()
    vargs = vars(args)
    gtf_fn = vargs['GTF']
    if not os.path.isfile(gtf_fn): parser.error('GTF file does not exist!')
    genome = vargs['genome']
    if not os.path.isfile(genome): parser.error('Genome file not found!')
    
    # Parse options
    genefield = vargs['genefield']
    ofn = vargs['outname']
    tables = vargs['tables']
    interval_fns = vargs['interval_fns']
    nprocess = vargs['nprocessors']
    gene_list = vargs['gene_list']
    genes_to_use = []
    if gene_list != None:
        if not os.path.isfile(gene_list): parser.error('Supplied gene list file does not exist!')
        else:
            genes_to_use = [x.strip() for x in open(gene_list).readlines()]
            genes_to_use = set(genes_to_use)
            print('Only considering %d genes from input gene list in analysis.'%(len(genes_to_use)))
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

    # Load genes 
    print('Loading GTF...')
    GTF = load_gtf(gtf_fn,genefield,gene_list,genes_to_use)
    print('GTF loaded.')
    print('Loading genes...')
    genes = []
    g_index_counter = 0
    for g in sorted(GTF.keys()):
        gene = Gene(GTF[g],g,g_index_counter)
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
    print('Loaded %d genes from input GTF file.'%(len(genes)))
    
    # load table data
    if len(tables) > 0:
        EXP = load_tables(tables)

    # Get sequence and expression data
    GENOME = FastaFile(genome)
    for g in genes:
        g.get_sequence(GENOME)
        g.get_seq_gc_cont()
        g.get_seq_CpG()
        
        if len(tables) > 0:
            samples = EXP['samples']
            if g.name in EXP['values']:
                g.get_exp(samples,EXP['values'][g.name])
            else: 
                for s in samples: g.expression[s] = 'NA'
    
    # If interval data is provided
    if len(int_fns) > 0:
        dones,was_done = [],0
        chunk_size = 1000
        num_chunks = int(ceil(len(genes)/chunk_size))
        chunks = []
        for i in range(0,num_chunks):
            rr = range(i*chunk_size,min((i+1)*chunk_size,len(genes)))
            chunks.append(genes[rr[0]:rr[-1]+1])
        print('  Divided %d genes into %d gene chunks.'%(len(genes),num_chunks))
        res = [pool.apply_async(get_interval_data,args=(c,int_fns),callback=dones.append) for c in chunks]
        while len(dones) != num_chunks:
            if len(dones)%10==0 and was_done != len(dones):
                was_done = len(dones)
                print('    %d of %d gene chunks complete.'%(len(dones),num_chunks))
        genes = []
        for gene in res:
            gget = gene.get()
            for g in gget: genes.append(g)
   
    # write output
    of = open(ofn,'w')
    header = ['Gene','full_length','coding_length','GC','CpG']
    if len(tables)>0:
        for s in EXP['samples']: header.append('%s'%(s))
    for ints in int_fns: header.append(ints[1])
    of.writelines('\t'.join(header)+'\n')
    for g in genes:
        tot_len = np.log2(g.total_length * 1e-3) # log2 of length in kb
        cod_len = np.log2(g.coding_length * 1e-3) # log2 of length in kb
        ostr = [g.name,str(tot_len),str(cod_len),'%0.3f'%(g.GC),'%0.3f'%(g.CpG)]
        if len(tables)>0:
            for s in EXP['samples']: ostr.append('%s'%(g.expression[s]))
        for ints in int_fns:
            samp = ints[1]
            if g.intervalData[samp] != None:
                ostr.append('%0.3f'%(g.intervalData[samp]))
            else: ostr.append(str(g.intervalData[samp]))
        ostr = '\t'.join(ostr)+'\n'
        of.writelines(ostr)
    of.close()

#####################
# UTILITY FUNCTIONS #
#####################
def load_gtf(gtf,genefield,gene_list,genes_to_use):
    '''
    Load and process genes in GTF file.
    '''
    genes = {}
    bad_genes = set() # list for problematic genes

    if gtf.endswith('.gz'): FH = gzip.open(gtf, 'rt')
    else: FH = open(gtf)

    for line in FH:
        if line.startswith('#'): continue # some GTFs have header lines
        l = line.strip().split('\t')
        chrom,typ,start,stop,strand,data = l[0],l[2],int(l[3]),int(l[4]),l[6],l[-1]

        if len(chrom.split('_'))>1: continue
        #if chrom == 'chrY' or chrom == 'Y' or chrom == 'chrM' or chrom == 'M': continue 

        # parse data string
        gene_id,tx_id = '',''
        for d in data.split(';'):
            if d=='': continue
            elif d[0] == ' ': d = d[1:]
            if d.startswith(genefield):
                d1 = d.split(' ')[1]
                gene_id = d1.replace('"','')
        
        # check if gene list provided
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
    print('  Deleting %d genes annotated to multiple chromosomes.'%(len(bad_genes)))
    for g in bad_genes:
        del(genes[g])

    # return
    return genes

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

def run_merge(exlist):
    '''
    Function to recursively run merge_list script until all overlapping intervals accounted for.
    '''

    mergelen = len(exlist)
    mergelist = merge_list(sorted(exlist))
    while len(mergelist) < mergelen and len(mergelist) > 1:
        mergelen = len(mergelist)
        mergelist = merge_list(mergelist) 
    return mergelist
     
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
                #continue # BUG
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

def get_interval_data(genes,INT):
    '''
    Get interval data for gene from input files in int_fns. 
    '''
    for fn,name in INT:
        tb = TabixFile(fn)

        for g in genes:
        
            # Get region for searching replication timing data
            g_len = g.total_length
            midp = round((g.start+g.stop)/2)
            min_width = 10e3 # search region at least 1 kb
            if g_len < min_width:
                start = midp - round(min_width/2)
                stop = midp + round(min_width/2)
                gstr = '%s:%d-%d'%(g.chrom,start,stop)
            else: gstr = '%s:%d-%d'%(g.chrom,g.start,g.stop)

            # Call to tabix to get dat from bedGraph 
            it_genes = tb.fetch(gstr)
            intData = []
            for itr in it_genes:
                if itr == '': continue
                itr = itr.split('\t')
                intData.append(float(itr[-1]))
            if len(intData)>0:
                g.intervalData[name] = np.mean(intData)
                continue
            else: 
                # Extend search if value not found
                extends0 = [50e3,100e3,500e3,1e6]
                extends = []
                for e in extends0: 
                    if e > g_len: extends.append(e)
                found = False
                for e in extends:
                    start = max(1,midp - round(e/2))
                    stop = midp + round(e/2)
                    gstr = '%s:%d-%d'%(g.chrom,start,stop)
                    
                    it_genes = tb.fetch(gstr)
                    for itr in it_genes:
                        if itr == '': continue
                        itr = itr.split('\t')
                        intData.append(float(itr[-1]))
                        found = True
                    if found == True:
                        g.intervalData[name] = np.mean(intData)
                        break
                
                if found == False:
                    g.intervalData[name] = None

    return genes

#####################
# Class definitions #
#####################
class Gene:
    '''
    Gene class.
    '''
    def __init__(self,gene,name,index):
        '''
        Initialize Gene object.
        '''
        self.name = name
        self.index = index

        # set exons
        if len(gene['exons'])==1: self.exons = gene['exons']
        else: self.exons = run_merge(sorted(gene['exons']))
        
        # set coding regions
        if len(gene['CDS']) == 1: self.coding_regions = gene['CDS']
        else: self.coding_regions = run_merge(sorted(gene['CDS']))
        
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
        self.total_length = self.stop - self.start + 1
        self.coding_length = 0
        for cr in self.coding_regions:
            self.coding_length += cr[1] - cr[0] + 1
        self.exonic_length = 0
        for ex in self.exons:
            self.exonic_length += ex[1] - ex[0] + 1
   
        # set up covariate values
        self.seq = None
        self.GC = None
        self.CpG = None
        self.expression = {}
        self.intervalData = {}

    def get_sequence(self,genome):
        ''' Get coding sequence. '''
        seq = ''
        for cr in self.coding_regions:
            gstr = '%s:%d-%d'%(self.chrom,cr[0],cr[1])
            seq += genome.fetch(region=gstr)
        self.seq = seq.upper()
    
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
        for i,s in enumerate(samples): self.expression[s] = vals[i]

#############
# Main call #
#############
if __name__ == '__main__': main()
