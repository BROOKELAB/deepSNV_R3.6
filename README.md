# deepSNV_R3.6
viral variant calling with deepSNV for R 3.6.

This is a nextflow pipeline that implements the steps of the lauringlab variant calling workflow.
See this site: https://github.com/lauringlab/variant_pipeline/tree/master/tutorial

Our nextflow pipeline uses deepSNV version 1.30 for R/3.6 and Bioconductor/3.9.
It has been tested on IGV's biocluster with this software stack:

nextflow/19.07.0-Java-1.8.0_152
MultiQC/1.7-IGB-gcc-4.9.4-Python-3.6.1
fastp/0.19.5-IGB-gcc-4.9.4
SAMtools/1.9-IGB-gcc-4.9.4
picard/2.10.1-Java-1.8.0_152
Bowtie2/2.3.2-IGB-gcc-4.9.4
R-lib/3.6-deepSNV
MUSCLE/3.8.31-IGB-gcc-4.9.4
Python/2.7.13-IGB-gcc-4.9.4
Biopython/1.68-IGB-gcc-4.9.4-Python-2.7.13
VarScan/2.3.9-Java-1.8.0_152
snpEff/4.3t-Java-1.8.0_152
snpEffDB/modified_PR8

