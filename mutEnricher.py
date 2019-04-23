'''
Main module for running somatic mutation enrichment analyses. 
'''

import sys, os
import argparse

def main():

    #########
    # PATHS #
    #########
    codeDir = os.path.abspath(os.path.dirname(__file__))
    mfDir = os.path.abspath(os.path.join(codeDir,'math_funcs'))
    sys.path.insert(0,mfDir)
    
    ###############################
    # ARGUMENT AND OPTION PARSING #
    ###############################

    usage = 'python %(prog)s'
    description = 'Perform somatic coding or non-coding analysis on sets of somatic mutation calls.'
    epilog = 'For command line options of sub-commands, type: %(prog)s COMMAND -h'
    version = "1.1.3"

    # set up parser and sub-parsers
    parser = argparse.ArgumentParser(usage=usage,description=description,epilog=epilog,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version',action='version',version='%(prog)s '+version)
    subparsers = parser.add_subparsers(dest='subcommand_name')

    # Add sub-parsers
    add_coding_parser(subparsers)
    add_noncoding_parser(subparsers)

    # Get args
    if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    
    # Run appropriate command
    subcommand = args.subcommand_name
    
    if subcommand == 'coding':
        from coding_enrichment import run
        run(parser, args, version)
    elif subcommand == 'noncoding':
        from noncoding_enrichment import run
        run(parser, args, version)

def add_coding_parser(subparsers):
    coding_parser = subparsers.add_parser('coding',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                    help="Perform somatic coding-gene analyses on non-silent mutations in genes.")

    coding_parser.add_argument('GTF',metavar='genes.gtf',type=str,
                        help='Input GTF file (Required). Can be provided as plain text or gzip-compressed file.')
    coding_parser.add_argument('vcfs',metavar='vcfs_list.txt',type=str,
                        help='Input VCFs list file (Required). Required columns: file path, sample name.\
                        NOTE: sample names must be unique for each sample!')

    coding_parser.add_argument('-o','--outdir',type=str,default='./',dest='outdir',
                        help='Provide output directory for analysis.')
    coding_parser.add_argument('--prefix',type=str,default='mutation_enrichment',dest='prefix',
                        help='Provide prefix for analysis.')
    
    coding_parser.add_argument('--gene-field',type=str,default='gene_id',dest='genefield',
                        help='Provide field name from input GTF containing gene name/id information.')
    coding_parser.add_argument('-g','--gene-list',type=str,default=None,dest='gene_list',
                               help='Provide list of genes to which analysis should be restricted (one gene per-line in text file). \
                               Analysis will only considers genes from GTF file that are present in this list. \
                               Default behavior is to query all coding genes present in input GTF.')
    
    coding_parser.add_argument('--bg-vars-type',dest='bg_vars_type',type=str,default='all',
                        help="Select which variants should be counted in background rate calculations. Choices are: 'all' and 'silent'. \
                        If 'all' is selected, all variants (silent + non-silent) are counted in background calculations. If 'silent' is \
                        selected, only silent mutations count towards background.")
 
    coding_parser.add_argument('--maf',type=str,dest='maf',default=None,
                        help='Instead of VCF list file, provide MAF (mutation annotation format) file with mutation information. \
                        To use, provide a dummy character (e.g. "-") for the VCFs argument and provide a MAF file with this option. Gene \
                        information (e.g. lengths) are computed from input GTF. Genes not present by genefield in GTF (read from first column \
                        of MAF) are skipped. Input MAF can be provided as plain text of gzip-compressed file.')
    coding_parser.add_argument('--exome-only',action='store_true',dest='exome_only',
                        help='If using exome-based data, choose this flag to only consider exonic coordinates of genes for background estimates. \
                        Default behavior is to consider full gene length (exons + introns) in calculations.')
    
    coding_parser.add_argument('--anno-type',type=str,default='illumina',dest='tType',
                        help='Select annotation type for determining non-silent somatic variants. Alternatively, provide tab-delimited input text \
                        file describing terms for use. Valid default options are: illumina, annovar. If providing text file, must include one \
                        term per row with 3 columns: 1) String that is either "Gene" or "Effect" to denote field with gene name or gene effect,\
                        respectively; 2) value from VCF INFO field for code to search for matching gene name or non-silent effect;\
                        3) valid terms (can be left blank for "Gene" row). \
                        If MAF input is used, this option is ignored and default MAF terms are used.')
    
    coding_parser.add_argument('-m','--mappable-regions',dest='map_regions',default=None,type=str,
                        help='Provide BED file of mappable genomic regions (sorted and tabix-indexed). If provided, only portions of regions \
                        from input file overlapping these mappable regions will be used in analsyis. Region lengths are also adjusted for \
                        enrichment calculations.')
    coding_parser.add_argument('-p','--processors',type=int,default=1,dest='nprocessors',
                        help='Set number of processors for parallel runs.')
    coding_parser.add_argument('--snps-only',dest='snps_only',action='store_true',
                        help='Set this flag to tell program to only consider SNPs in analysis. Default is to consider all variant types.')
    coding_parser.add_argument('-c','--covariates-file',type=str,default=None,dest='cov_fn',
                        help='Provide covariates file. Format is tab-delimited text file, with first column listing gene name according to \
                        gene_id field in input GTF. Header should contain covariate names in columns 2 to end.')
    coding_parser.add_argument('-w','--covariate-weights',type=str,default=None,dest='weights_fn',
                        help='Provide covariates weight file. Format is tab-delimited file (no header) with: covariate name, weight. \
                        Weights are normalized to sum=1. If not provided, uniform weighting of covariates is assumed.') 
    coding_parser.add_argument('--by-contig',dest='by_contig',action='store_true',
                        help='Use this flag to perform clustering on genes by contig (i.e. by chromosome). This speeds computation of gene \
                        clusters. If not set, clusters are computed using all genes in same run.')
    coding_parser.add_argument('--use-local',dest='use_local',action='store_true',
                        help='Use this flag to tell program to use local gene background rate instead of global background rate. \
                        If covariate files or pre-computed covariates supplied, this option is ignored.')
    coding_parser.add_argument('--min-clust-size',type=int,default=3,dest='min_clust_size',
                        help='Set minimum number of covariate cluster members. Regions belonging to a cluster with only itself or less than\
                        this value are flagged and a local background around the region is calculated and used instead.')
    coding_parser.add_argument('--precomputed-covars',type=str,default=None,dest='cov_precomp_dir',
                        help='Provide path to pre-computed covariate clusters for regions in input BED file.')
    coding_parser.add_argument('-d','--hotspot-distance',type=int,default=50,dest='max_hs_dist',
                        help='Set maximum distance between mutations for candidate hotspot discovery.')
    coding_parser.add_argument('--min-hs-vars',type=int,default=3,dest='min_hs_vars',
                        help='Set minimum number of mutations that must be present for a valid candidate hotspot.')
    coding_parser.add_argument('--min-hs-samps',type=int,default=2,dest='min_hs_samps',
                        help='Set minimum number of samples that must contain mutations to inform a valid candidate hotspot.')
    coding_parser.add_argument('--blacklist',type=str,default=None,dest='blacklist_fn',
                        help='Provide a blacklist of specific variants to exclude from analysis. Blacklist file format is tab-delimited\
                        text file with four required columns: contig (chromosome), position (1-indexed), reference base, alternate base.')
    
    coding_parser.add_argument('--ap-iters',type=int,default=1000,dest='ap_iters',
                        help='Set maximum number of AP iterations before re-computing with alternate self-similarity.')
    coding_parser.add_argument('--ap-convits',type=int,default=50,dest='ap_convits',
                        help='Set number of convergence iterations for AP runs (i.e. if exemplars remain constant for this many iterations,\
                              terminate early). This value MUST be smaller than the total number of iterations.')
    coding_parser.add_argument('--ap-algorithm',type=str,default='fast',dest='ap_alg',
                        help="Select between one of two versions of AP clustering algorithm: 'slow' or 'fast'. The 'fast' version is faster\
                              in terms of runtime but consumes more memory than 'slow'.")

def add_noncoding_parser(subparsers):
    noncoding_parser = subparsers.add_parser('noncoding',
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                       help="Perform noncoding regional and hotspot somatic enrichment analysis.")  
    
    noncoding_parser.add_argument('regions',metavar='regions.bed',type=str,
                        help='Input regions BED file (Required). Can be provided as plain text or gzip-compressed file. \
                        Required columns: contig, 0-based start, 1-based end. A name for the region can be supplied in the 4th column.')
    noncoding_parser.add_argument('vcfs',metavar='vcfs_list.txt',type=str,
                        help='Input VCFs list file (Required). Required columns: file path, sample name.\
                        NOTE: sample names must be unique for each sample!')

    noncoding_parser.add_argument('-o','--outdir',type=str,default='./',dest='outdir',
                        help='Provide output directory for analysis.')
    noncoding_parser.add_argument('--prefix',type=str,default='mutation_enrichment',dest='prefix',
                        help='Provide prefix for analysis.')
    noncoding_parser.add_argument('-m','--mappable-regions',dest='map_regions',default=None,type=str,
                        help='Provide BED file of mappable genomic regions (sorted and tabix-indexed). If provided, only portions of regions \
                        from input file overlapping these mappable regions will be used in analsyis. Region lengths are also adjusted for \
                        enrichment calculations.') 
    noncoding_parser.add_argument('-p','--processors',type=int,default=1,dest='nprocessors',
                        help='Set number of processors for parallel runs.')
    noncoding_parser.add_argument('--snps-only',dest='snps_only',action='store_true',
                            help='Set this flag to tell program to only consider SNPs in analysis. Default is to consider all variant types.')
    noncoding_parser.add_argument('-c','--covariates-file',type=str,default=None,dest='cov_fn',
                        help='Provide covariates file. Format is tab-delimited text file, with first column listing regions in format \
                        <contig>:<1-based start>-<1-based end> and remaining columns with values. Header should contain covariate names \
                        in columns 2 to end.')
    noncoding_parser.add_argument('-w','--covariate-weights',type=str,default=None,dest='weights_fn',
                        help='Provide covariates weight file. Format is tab-delimited file (no header) with: covariate name, weight. \
                        Weights are normalized to sum=1. If not provided, uniform weighting of covariates is assumed.')
    noncoding_parser.add_argument('--use-local',dest='use_local',action='store_true',
                        help='Use this flag to tell program to use local region background rate instead of global background rate. \
                        If covariate files or pre-computed covariates supplied, this option is ignored.')
    noncoding_parser.add_argument('--min-rclust-size',type=int,default=3,dest='min_rclust_size',
                        help='Set minimum number of covariate cluster members. Regions belonging to a cluster with only itself or less than\
                        this value are flagged and a local background around the region is calculated and used instead.')
    noncoding_parser.add_argument('--precomputed-covars',type=str,default=None,dest='cov_precomp_dir',
                        help='Provide path to pre-computed covariate clusters for regions in input BED file.')
    noncoding_parser.add_argument('-d','--hotspot-distance',type=int,default=50,dest='max_hs_dist',
                        help='Set maximum distance between mutations for candidate hotspot discovery.')
    noncoding_parser.add_argument('--min-hs-vars',type=int,default=3,dest='min_hs_vars',
                        help='Set minimum number of mutations that must be present for a valid candidate hotspot.')
    noncoding_parser.add_argument('--min-hs-samps',type=int,default=2,dest='min_hs_samps',
                        help='Set minimum number of samples that must contain mutations to inform a valid candidate hotspot.')
    noncoding_parser.add_argument('--blacklist',type=str,default=None,dest='blacklist_fn',
                        help='Provide a blacklist of specific variants to exclude from analysis. Blacklist file format is tab-delimited\
                        text file with four required columns: contig (chromosome), position (1-indexed), reference base, alternate base.')
    
    noncoding_parser.add_argument('--no-wap',action='store_true',dest='no_wap',
                        help='Select flag to skip weighted average proximity (WAP) procedure.')
    
    noncoding_parser.add_argument('--ap-iters',type=int,default=1000,dest='ap_iters',
                        help='Set maximum number of AP iterations before re-computing with alternate self-similarity.')
    noncoding_parser.add_argument('--ap-convits',type=int,default=50,dest='ap_convits',
                        help='Set number of convergence iterations for AP runs (i.e. if exemplars remain constant for this many iterations,\
                              terminate early). This value MUST be smaller than the total number of iterations.')
    noncoding_parser.add_argument('--ap-algorithm',type=str,default='fast',dest='ap_alg',
                        help="Select between one of two versions of AP clustering algorithm: 'slow' or 'fast'. The 'fast' version is faster\
                              in terms of runtime but consumes more memory than 'slow'.")

if __name__ == '__main__': main()
