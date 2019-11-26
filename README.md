# viral variant calling with deepSNV for R 3.6.

This is a nextflow pipeline that implements the steps of the lauringlab variant calling workflow.
See this site: https://github.com/lauringlab/variant_pipeline/tree/master/tutorial

# DESCRIPTION

insert text here

# WORKFLOW OF THE PIPELINE

insert text here

# DEPENDENCIES

Our nextflow pipeline uses deepSNV version 1.30 for R/3.6 and Bioconductor/3.9.

It has been tested on IGV's biocluster with this software stack:

- <b>Nextflow</b>    tested with version 19.07.0 ( download page https://github.com/nextflow-io/nextflow/releases )
- <b>MultiQC</b>     tested with version 1.7 ( download page https://github.com/ewels/MultiQC )
- <b>fastp</b>      tested with version 0.19.5 ( download page https://github.com/OpenGene/fastp )
- <b>samtools</b>      tested with version 1.9 ( download page https://sourceforge.net/projects/samtools/files/samtools/ ) 
- <b>PicardTools</b>   tested with version 2.10.1 ( download page https://broadinstitute.github.io/picard/ )
- <b>Bowtie2</b>     tested with version 2.3.2 ( download page  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml )
- <b>deepSNV</b>     tested with R-lib 3.6 ( download page https://bioconductor.org/packages/release/bioc/html/deepSNV.html )
- <b>MUSCLE</b>      tested with version 3.8.31 ( download page https://www.drive5.com/muscle/downloads.htm )
- <b>Python</b>      tested with  version  2.7.13 ( download page https://www.python.org/downloads/ )
- <b>BioPython</b>      tested with  version 1.68 ( download page https://biopython.org/wiki/Download )
- <b>VarScan</b>      tested with version 2.3.9 ( download page http://varscan.sourceforge.net/ )
- <b>SnpEff</b>      tested with version 4.3t ( download page http://snpeff.sourceforge.net/download.html )
- <b>SnpEff database</b>      custome-built database for avial flu virus sp. PR8


# INSTALLATION INSTRUCTIONS

- Install all dependencies first. You may need to have root access to install some of these tools/languages on a cluster.
- Install NextFlow by running the following commands:

<pre>

# Make sure that Java v1.7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin

# OR make a system-wide installation (root user)
sudo mv nextflow /usr/local/bin

</pre>

- Do not forget to launch the 'hello world' nextflow pipeline (as per https://www.nextflow.io/) to make sure it works fine.
- Install this pipeline: The pipeline itself does not need installation. Simply copy this repo to a local folder and nextflow should be able to run it


# RUNNING THE PIPELINE

- This pipeline expects TWO samples as input: one is the control sample and the other one is treatment sample.
- Each sample of viral RNA is made up of short Illumina reads (single-end or paired-end). We tested the pipeline with avian flu strains.
- The sample(s) to be analyzed by this pipeline must be placed together in the same folder and must have a similar file naming patter; for example, all files should end in fastq | fq | fastq.gz | fq.gz
- Prepare a configuration file.  Some examples of configuration files are provided in the folder <b>config-files/</b>
- To run the pipeline type this command at the prompt: 

<pre>
nextflow -c config.file deepSNV_R3.6_pipeline_v2.nf
</pre>

# OUTPUTS

Nextflow generates two folders to keep track of execution progress. You can delete them once the execution ends successfully. They are called <i>.nextflow/ </i> and <i>work/ </i>

The actual results of the pipeline are placed in these folders:

- <b> Prepared-Reads /</b>  contains the results of QC, filter and trim of the raw reads with fastp
- <b> Bowtie2-align-n-sort /</b> contains the results of the alignment step wwith Bowtie2
- <b> MarkDup /</b> contains the results of the deduplication step with Picard's MarkDuplicates
- <b> ReadyBAMs /</b> contains the results of the filtering step with Samtools
- <b> deepSNV_variants /</b> contains the results of variant discovery with deepSNV
- <b> VarScan_SnpEff_variants /</b> contains the results of variant discovery with VarScan and annotation with SnpEff
- <b> MultiQC /</b> contains the results of read tracking with MultiQC

The final results are files with extension *.final.annotated.vcf* inside the <b>VarScan_SnpEff_variants</b> folder.
Each treatment sample should have one such file in this folder, even if it is empty of variants.


# LICENSE


This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>

# CITATIONS






