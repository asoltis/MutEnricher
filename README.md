###############
# MutEnricher #
###############

Author: Anthony R. Soltis (anthony.soltis.ctr@usuhs.edu, anthonyrsoltis@gmail.com)
Institution: Uniformed Services University of the Health Sciences, Bethesda, MD
License: MIT License, see LICENSE.txt
Version: 1.0.0

################
# Introduction #
################

MutEnricher is a flexible toolset that performs somatic mutation enrichment analysis of both 
protein-coding and non-coding genomic loci from whole genome sequencing (WGS) data. 

MutEnricher contains two distinct modules:
    1. coding - for performing somatic enrichment analysis of non-silent variation in protein-coding genes
    2. noncoding - for performing enrichment analysis of non-coding regions

The main driver script is mutEnricher.py and each tool can be evoked from here, i.e.:
    1. python mutEnricher coding ...
    2. python mutEnricher noncoding ...

See help pages and associated documentation for methodological and run details. 

################
# Installation #
################

1. Python packages

    MutEnricher has been tested with various Python 2.7.x versions. 
    Beyond the standard Python libraries, MutEnricher requires:

        1. NumPy
        2. SciPy
        3. Cython
        4. cyvcf2
        5. pysam

    These modules are readily available through the Python package manager pip. If pip is not already installed on your system, 
    see https://pip.pypa.io/en/stable/installing/ for instructions on how to obtain it. Alternatively, each module can be manually 
    downloaded and installed on your system.

2. Cythonize helper math functions

    MutEnricher calls several custom mathematical functions that need to be "cythonized" prior to use. To do this:

        1. From the main install directory (i.e. this directory), cd into the math_funcs sub-directory
        2. Run the following command: 
            python setup.py build_ext --inplace

    Successful completion of the above steps will generate math_funcs.so, which is used by MutEnricher. 

---
You can test the install by executing the help commands for the various tools, e.g.:
    python mutEnricher.py -h
    python mutEnricher.py coding -h
    python mutEnricher.py noncoding -h

########################
# Additional utilities #
########################

In the "utilities" sub-directory, we include two helper functions for generating covariate files for use with MutEnricher's 
covariate clustering functions:

    1. get_gene_covariates.py  
    2. get_region_covariates.py

See the help pages for example usage. (1) above requires GTF input (as for the coding module) and (2) requires and input BED (as for 
the noncoding module). Both also require a copy of an indexed genome FASTA file (e.g. for hg19/hg38 human genomes) as input.

################
# Example data #
################

We include various example files for testing MutEnricher on synthetic somatic data. See the "example_data" sub-folder.

Files/folders:
    
    1. example_data/annotation_files

        Contains example GTF and BED files for running MutEnricher's coding and noncoding modules. 

           (NOTE: uncompress with gunzip prior to use)
        2. ucsc.refFlat.20170829.promoters_up2kb_downUTR.no_chrMY.bed

    2. example_data/covariates

        Contains example covariate and covariate weights files for running the covariate clustering background method:

        For coding:
            1. ucsc.refFlat.20170829.no_chrMY.covariates.txt.gz
               (NOTE: uncompress with gunzip prior to use)
            2. ucsc.refFlat.20170829.no_chrMY.covariate_weights.txt
        For noncoding:
            1. ucsc.refFlat.20170829.promoters_up1kb_down200.no_chrMY.covariates.txt.gz
               (NOTE: uncompress with gunzip prior to use)
            2. ucsc.refFlat.20170829.promoters_up1kb_down200.no_chrMY.covariate_weights.txt

    3. nonsilent_terms.txt

        Example non-silent terms file for use with coding module. This example is applicable to VCFs annotated with ANNOVAR refGene models
        (the sample VCFs are annotated in this way). Use with the --anno-type option in the coding module.

        NOTE: These same terms will be used if "annovar" is passed to the --anno-type option. 

    4. precomputed_apcluster

        This folder provides pre-computed affinity propagation results for the datasets in (1) and (2) above. These directories can be
        supplied to MutEnricher via the --precomputed-covars option. 

        For coding:
            - coding.ucsc.refFlat.20170829.no_chrMY
                - for all gene clustering results: all_genes
        For noncoding:
            - noncoding.ucsc.refFlat.20170829.promoters_up2kb_downUTR.no_chrMY/apcluster_regions

    5. quickstart_commands.txt
        
        Sample execution commands (associated with quickstart guide).

    6. vcf_files.txt

        Sample VCF input files list file. This file contains local paths and assumes working directory is "example_data" sub-directory.

    7. vcfs

        Sub-directory containing 100 synthetic somatic VCF files (compressed with index .tbi files). These files were generated by randomly
        inserting "somatic mutations" at positions in the hg19 genome at a target rate of ~2 mutations/Mb. Three true positive cases 
        are included, two coding and one non-coding, whereby non-silent mutations were inserted into the TP53 and KRAS genes and somatic
        mutations were inserted into the TERT gene promoter region. 

    8. doc

        Tutorial and quick start documentation are included is this sub-directory. 

##############
# Change log #
##############
# 06-15-2018 
- Initial release; The development of this Software was sponsored by the Uniformed Services University of the Health Sciences (USU); however, the information or content and conclusions do not necessarily represent the official position or policy of, nor should any official endorsement be inferred on the part of, USU, the Department of Defense, or the U.S. Government. 

