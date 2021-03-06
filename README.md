# MutEnricher #
----------------

<img src='images/MutEnricher_Fig1_20190422_for_GitHub.png'>

Author: Anthony R. Soltis (anthony.soltis.ctr@usuhs.edu, anthonyrsoltis@gmail.com)

Institution: Uniformed Services University of the Health Sciences, Bethesda, MD

License: MIT License, see [License](https://github.com/asoltis/MutEnricher/blob/master/LICENSE.txt)

Version: 1.3.3

Introduction:
---------------

MutEnricher is a flexible toolset that performs somatic mutation enrichment analysis of both 
protein-coding and non-coding genomic loci from whole genome sequencing (WGS) data, implemented
in Python and **usable with Python 2 and 3.** 

**MutEnricher is now also available as a** [Docker image](https://hub.docker.com/r/asoltis/mutenricher).

MutEnricher contains two distinct modules:
1. coding - for performing somatic enrichment analysis of non-silent variation in protein-coding genes
2. noncoding - for performing enrichment analysis of non-coding regions

The main driver script is mutEnricher.py and each tool can be evoked from here, i.e.:
1. python mutEnricher coding ...
2. python mutEnricher noncoding ...

See help pages and associated documentation for methodological and run details. 

Citation:
---------
A [MutEnricher manuscript](https://rdcu.be/b51ka) is now published in BMC Bioinformatics. Please cite if using this software:

Soltis, A.R., Dalgard, C.L., Pollard, H.B., & Wilkerson, M.D. MutEnricher: a flexible toolset for somatic mutation enrichment analysis of tumor whole genomes. BMC Bioinformatics (2020). 20(1).

Info and User Guides:
---------------------

[Wiki](https://github.com/asoltis/MutEnricher/wiki)

[Quickstart guide](https://github.com/asoltis/MutEnricher/wiki/Quickstart-guide)

[Tutorial](https://github.com/asoltis/MutEnricher/wiki/Tutorial)

[Output file descriptions](https://github.com/asoltis/MutEnricher/wiki/Output-file-descriptions)

Installation:
---------------

See [Installation Guide](https://github.com/asoltis/MutEnricher/wiki/Installation-Guide) section on Wiki.

Additional utilities
----------------------

In the "utilities" sub-directory, we include two helper functions for generating covariate files for use with MutEnricher's 
covariate clustering functions:

    1. get_gene_covariates.py  
    2. get_region_covariates.py

See the help pages for example usage. (1) above requires GTF input (as for the coding module) and (2) requires and input BED (as for 
the noncoding module). Both also require a copy of an indexed genome FASTA file (e.g. for hg19/hg38 human genomes) as input.

Example data
--------------

We include various example files for testing MutEnricher on synthetic somatic data. See the "example_data" sub-folder. 

Several quickstart commands are provided in example_data/quickstart_commands.txt file. A sample quickstart command for coding analysis:

```
cd example_data
python ../mutEnricher.py coding annotation_files/ucsc.refFlat.20170829.no_chrMY.gtf.gz vcf_files.txt --anno-type nonsilent_terms.txt -o test_out_coding --prefix test_global
```

Files/folders contained in example_data:
    
1. example_data/annotation_files

    Contains example GTF and BED files for running MutEnricher's coding and noncoding modules. 
    - ucsc.refFlat.20170829.no_chrMY.gtf.gz
    - ucsc.refFlat.20170829.promoters_up2kb_downUTR.no_chrMY.bed

    NOTE: Input GTF (coding analysis) and BED files (noncoding analysis) can be gzip compressed or not. 

2. example_data/covariates

    Contains example covariate and covariate weights files for running the covariate clustering background method:

    For coding:
    - ucsc.refFlat.20170829.no_chrMY.covariates.txt
    - ucsc.refFlat.20170829.no_chrMY.covariate_weights.txt
    
    For noncoding:
    - ucsc.refFlat.20170829.promoters_up1kb_down200.no_chrMY.covariates.txt
    - ucsc.refFlat.20170829.promoters_up1kb_down200.no_chrMY.covariate_weights.txt

3. nonsilent_terms.txt

    Example non-silent terms file for use with coding module. This example is applicable to VCFs annotated with ANNOVAR refGene models
    (the sample VCFs are annotated in this way). Use with the --anno-type option in the coding module.

    NOTE: These same terms will be used if "annovar" is passed to the --anno-type option. 

4. precomputed_apcluster

    This folder provides pre-computed affinity propagation results for the datasets in (1) and (2) above. These directories can be
    supplied to MutEnricher via the --precomputed-covars option. 

    For coding (all genes):
    - coding.ucsc.refFlat.20170829.no_chrMY/all_genes
    
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

# Change log #
---------------
06-15-2021
----------
- Version 1.3.3
- Updates:
    - Include VEP annotation parsing capabilities (via "CSQ" field) in coding module. 
    - Included missing function in coding analysis code to parse blacklist variant input file.

05-11-2021
----------
- Version 1.3.2
- Updates:
    - Included SnpEff annotation parsing capabilities (via "ANN" INFO field) in coding module. Set --anno-type options to 'SnpEff' to 
      use pre-set annotations compatible with this tool.
    - Improved error handling for interval files and regions in covariate utility scripts.

10-01-2020
----------
- Version 1.3.1
- Bug fix:
    - Update to coding module and gene covariate code to address incomplete merging of overlapping gene feature intervals (exons, CDS).

06-10-2020
----------
- Dockerfile added for creation of Docker image.
- No code updates.

10-23-2019
----------
- Version 1.3.0
- Major updates:
    - 'nsamples' (binomial testing method) is now default statistical testing (--stat-type) option.
    - Combined covariate clustering plus local background rate method implemented. When covariates are supplied and --use-local is also set,
      programs compute local backgrounds around features part of clusters during background calculations. 

10-10-2019
----------
- Version 1.2.1
- Minor update to local background method, whereby minimum search window is increased to 1 Mb.

09-13-2019
----------
- Version 1.2.0
- Major updates:
    - Code updated for compatibility with Python 3.
    - Included --stat-type option to select between original negative binomial test based on mutation counts (nmutations, default) or
      binomial test on number of mutated samples (nsamples).
- Minor updates:
    - Updated --anno-type preset options to better reflect various ANNOVAR gene annotations. 
    - Deprecated --repliseq-fns option in utilities code and updated to -i/--interval-files option

03-25-2019
----------
- Version 1.1.3
- Updates:
    - Noncoding code now produces <prefix>_region_WAP_hotspot_Fisher_enrichments.txt output file, which includes an overall combined
      Fisher's combined p-value for the overall region, WAP, and hotspot (if present) p-values.

02-12-2019
----------
- Version 1.1.2
- Updates:
    - In both coding and noncoding modules, new option --min-hs-samps included for setting minimum number of samples that must contain 
      mutations in a candidate hotspot region for subsequent testing. Default is set to 2; setting to 1 is equivalent to prior default 
      behavior. 

01-15-2019
-----------
- Version 1.1.1
- Updates/bug fixes:
    - Coding analysis code now produces output file with combined Fisher p-value for overall gene and hotspot(s) enrichments.
    - Updated method used to compute Fisher p-values for better numerical accuracy.
    - utilities/get_gene_covariates.py updated to read gzipped GTF files.
    - Fixed minor bug in coding analysis code associated with local background rate calculation method.
    - Updated coding analysis code to calculate gene background mutation rate from samples possessing at least one non-silent mutation.

06-15-2018 
-----------
- Initial release; The development of this Software was sponsored by the Uniformed Services University of the Health Sciences (USU); however, the information or content and conclusions do not necessarily represent the official position or policy of, nor should any official endorsement be inferred on the part of, USU, the Department of Defense, or the U.S. Government. 

