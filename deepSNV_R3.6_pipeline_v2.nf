#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
               deepSNV variant pipeline part1
========================================================================================
*/

// version
version                     = 0.12

// Credit to Phil Ewels for this segment to generate the help message at the command-line

def helpMessage() {
    log.info"""
    =========================================
     deepSNV variant pipeline  part1 v${version} 
    =========================================

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run -c <conf>  <deepSNV_pipeline_part1.nf>
    
    where conf is the configuration file for a particular experiment
    
    To override existing values from the command line, please use these parameters:

    Required options
    --reads                 Input data( this string must be surrounded with quotes and should contain a full path and regular expression). PE reads expected
    --genome                Genome in fasta format. (must be surrounded by quotes)
    --outdir                Full path of output folder. (must be surrounded by quotes)
    
    Read preparation options
    --singleEnd             options: true|false. true for single aka SE reads; false for paired aka PE reads. Default: false.
    --guess_adapter         options: true|false. auto-detect adapter from input file. Default: false.
    --trim_lowC             options: true|false. Filter low complexity reads. Default: false.
    --trim_polyG            options: true|false. Trim polyGs from the end of the read. Default: false.
    --min_read_length       positive integer. minimum read length. (must be surrounded by quotes)
    --min_BQ                positive integer [0-40] minimum base quality. Default: 20. (must be surrounded by quotes)
    --forward_adaptor       string with sequence of forward adaptor for R1. (must be surrounded by quotes)
    --reverse_adaptor       string with sequence of reverse adapter for R2. (must be surrounded by quotes)
    
    Alignment and BAM filtering options
    --bowtie2params         Bowtie2 alignment parameters. (must be surrounded by double quotes)
    --min_MQ                minimum mapping quality. Default: 0 to disable this filter.
    --RGLB                  LB field of RG line. Default: TruSeq
    --RGPL                  PL field of RG line. Default: illumina
    --RGPU                  PU field of RG line. Default: Novaseq6000
    --RGCN                  CN field of RG line. Default: CBC

    Variant discovery and annotation options
    --skipAnnotation        options: true|false. Skip this part of the workflow. Default: false.
    --method                adjustment method for multiple testing corrections. deepSNV option. Options: BH | bonferroni. Default: bonferroni
    --P_CUT                 p value cutoff. deepSNV option. Default: 0.01
    --combinePvalMethod     method for combining two p values. deepSNV option. Options: fisher | average | max. Default: fisher 
    --dispersion            alternative dispersion. deepSNV option. Options: two.sided | one.sided | bin. Default: two.sided 
    --stringentFreq         stringency frequency. deepSNV option. Default: 0.01
    --minCoverage           minimum coverage. VarScan option. Default: 3
    --maxDepth              maximum read depth. Default: 0 to disable this filter 
    --snpEffDB              name of snpEff database to use for annotation. Default: modified_PR8

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Override these values from the command line

// Required parameters
params.reads                = false
params.genome               = false
params.outdir               = false
params.config               = false

// read prep options
params.singleEnd            = false              /* true for SE read. false for PE reads */
params.guess_adapter        = false              /* if true, fatsp will guess adaptors from the first 1M reads */
params.min_read_length      = false              /* minimum read length after trimming */
params.min_BQ               = '20'               /* minimum base quality */
params.forward_adaptor      = false              /* sequence of forward adaptor surrounded by quotes */
params.reverse_adaptor      = false              /* sequence of reverse adaptor surrounded by quotes */
params.trim_lowC            = true               /* Filter low complexity reads */
params.trim_polyG           = false              /* Trim polyGs from  read */
params.RGLB                 = "TruSeq"           /* LB field of RG line */
params.RGPL                 = "illumina"         /* PL field of RG line */
params.RGPU                 = "Novaseq6000"      /* PU field of RG line */
params.RGCN                 = "CBC"              /* CN field of RG line */
params.RGSN                 = false              /* SN field of RG line */
params.RGID                 = false              /* ID field of RG line */

// alignment options
params.bowtie2params        = false
params.min_MQ               = '0'                /* minimum mapping quality */

// deepSNV options
params.skipAnnotation       = false
params.method               = "bonferroni"
params.P_CUT                = '0.01'
params.combinePvalMethod    = "fisher"
params.dispersion           = "two.sided"
params.stringentFreq        = '0.01'

// pileup, VarScan options
params.minCoverage          = '3'
params.maxDepth             = '0'               /* maximum read depth. Set it to 0 to disable the filter */

// other options
params.email                = false
params.pythonScriptDir      = false

// variables for indices
params.genomesDir           = "/igbgroup/groups/hpcbio/projects/cbrooke/2018-Nov-VariantPipeline/data/reference/nf-produced"
genomeFile                  = file(params.genome)
genomePrefix                = genomeFile.getBaseName()
genomeStore                 = "${params.genomesDir}/${genomePrefix}"

// cluster-specific variables
params.executor             = "slurm"
params.Queue                = 'normal'
params.trimTime             = "06:00:00"
params.trimThreads          = '4'
params.trimMem              = '20'
params.alignTime            = '06:00:00'
params.alignThreads         = '10'
params.alignMem             = '250'
params.dedupTime            = '06:00:00'
params.dedupThreads         = '4'
params.dedupMem             = '50'
params.vcallTime            = '06:00:00'
params.vcallThreads         = '4'
params.vcallMem             = '100'
params.RMem                 = '100'
params.RTime                = '06:00:00'
params.RThreads             = '4'
params.mem                  = '10'
params.cpus                 = '1'

// Software stack
params.multiqcMod           = 'MultiQC/1.7-IGB-gcc-4.9.4-Python-3.6.1'
params.fastpMod             = 'fastp/0.19.5-IGB-gcc-4.9.4'
params.samtoolsMod          = 'SAMtools/1.9-IGB-gcc-4.9.4'
params.picardMod            = 'picard/2.10.1-Java-1.8.0_152'
params.bowtie2Mod           = 'Bowtie2/2.3.2-IGB-gcc-4.9.4'
params.deepSNV_RMod         = 'R-lib/3.6-deepSNV'
params.R_lib                = "/igbgroup/groups/hpcbio/projects/cbrooke/2018-Nov-VariantPipeline/src/R3.6_lib"
params.scriptdir            = "/home/groups/hpcbio/projects/cbrooke/2018-Nov-VariantPipeline/src"
params.MuscleMod            = "MUSCLE/3.8.31-IGB-gcc-4.9.4"
params.PythonMod            = "Python/2.7.13-IGB-gcc-4.9.4"
params.BiopythonMod         = "Biopython/1.68-IGB-gcc-4.9.4-Python-2.7.13"
params.VarScanMod           = "VarScan/2.3.9-Java-1.8.0_152"
params.snpEffPath           = "/home/groups/hpcbio/projects/cbrooke/2018-Nov-VariantPipeline/data/snpEff/"
EBROOTSNPEFF                = "/home/groups/hpcbio/projects/cbrooke/2018-Nov-VariantPipeline/data/snpEff/"
params.snpEffDB             = "modified_PR8"
params.JavaOptions          = " -Xmx4g -XX:ParallelGCThreads=2 "

// Sanity checks
if (!params.outdir)                         exit 1, "Must set --outdir with path to results"
if (!params.reads)                          exit 1, "PE reads expected. Must set --reads with a valid regexp"
if (!params.genome)                         exit 1, "Must set --genome file for indexing and alignment"
if (!params.config)                         exit 1, "Must set --config  tsv file with three columns control, input, yaml-options"

if (!params.guess_adapter && !params.singleEnd && !params.forward_adaptor &&  !params.reverse_adaptor)  exit 1, "Must set --forward_adaptor and --reverse_adaptor "
if (!params.guess_adapter &&  params.singleEnd && !params.forward_adaptor) exit 1, "Must set --forward_adaptor"
if (!params.guess_adapter && !params.singleEnd && !params.forward_adaptor &&  params.reverse_adaptor)  exit 1, "Must set --forward_adaptor"
if (!params.guess_adapter && !params.singleEnd &&  params.forward_adaptor && !params.reverse_adaptor)  exit 1, "Must set --reverse_adaptor"

// variables for fastp
adapterOptionsPE          = params.guess_adapter ? " --detect_adapter_for_pe " : " --detect_adapter_for_pe --adapter_sequence=${params.forward_adaptor}  --adapter_sequence_r2=${params.reverse_adaptor}  "
adapterOptionsSE          = params.guess_adapter ? " " : " --adapter_sequence=${params.forward_adaptor} "
trimOptions2              = params.trim_lowC     ? ' --low_complexity_filter ' : ' ' 
trimOptions3              = params.trim_polyG    ? ' --trim_poly_g ' : ' --disable_trim_poly_g ' 
params.trimOptions        = " -5 -3 --cut_mean_quality ${params.min_BQ} --length_required ${params.min_read_length} "

// Validate inputs
config = file(params.config)
if( !config.exists() ) exit 1, "Missing config file: '$params.config'. Specify file with --config"

/*
 * Create channels for config file
 */
Channel
    .from(config.readLines())
    .map { line ->
        list         = line.split(',')
        test_id      = list[0]
        control_id   = list[1]
        options_file = list[2]
        [ test_id, control_id, options_file ]
    }
    .into{ deepSNV_param; recipVar_param; posStat_param  }

    
// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Header log info

runInfo = """

------------------------------------------------------------------------- 
deepSNV pipeline  part1 v${version}
------------------------------------------------------------------------- 
Pipeline ver               : ${version}

REQUIRED ARGUMENTS
Reads                      : ${params.reads}
Genome                     : ${params.genome}
Genome Store               : ${genomeStore}

TOOLS
Bowtie2                    : ${params.bowtie2Mod}
fastp                      : ${params.fastpMod}
samtools                   : ${params.samtoolsMod}
Picard                     : ${params.picardMod}
deepSNV                    : ${params.deepSNV_RMod}
R_lib                      : ${params.R_lib}
Python                     : ${params.PythonMod}
Biopython                  : ${params.BiopythonMod}
Muscle                     : ${params.MuscleMod}

READ PREPARATION OPTIONS
read length                : ${params.read_length}
min read len after trim    : ${params.min_read_length}
if(!params.guess_adapter) {
guess.adaptor              : false
forward adaptor            : ${params.forward_adaptor}
reverse adaptor            : ${params.reverse_adaptor}
}
if(params.guess_adapter) {
guess.adaptor              : true
}
if(params.singleEnd) {
Read type                  : SE
adaptor trimming options   : ${adapterOptionsSE}
}
if(!params.singleEnd) {
Read type                  : PE
adaptor trimming options   : ${adapterOptionsPE}
}
other filter/trimming options  : ${params.trimOptions} ${trimOptions2} ${trimOptions2}


ALIGNMENT OPTIONS
Bowtie2 options            : ${params.bowtie2params}
samtools filter options    : -q${params.min_MQ} -F 1024

DEEPSNV OPTIONS
method                     : ${params.method}
p valuue cutoff            : ${params.P_CUT}
combine Pval method        : ${params.combinePvalMethod}
dispersion                 : ${params.dispersion}
stringent freq cutoff      : ${params.stringentFreq}

VARIANT FILTERING OPTIONS
snpEff database            : ${params.snpEffDB}
min coverage               : ${params.minCoverage}
p value cutoff             : ${params.P_CUT}
min_var_freq               : ${params.stringentFreq}
minimum base quality       : ${params.min_BQ}
mimumum mapping quality    : ${params.min_MQ}

OTHER OPTIONS
Current home               : $HOME
Current user               : $USER
Current path               : $PWD
Script dir                 : ${params.scriptdir}
Working dir                : $workDir
Output dir                 : ${params.outdir}
------------------------------------------------------------------------- 
"""

println(runInfo)

/*
 *
 * Phase 1: Genome preparation
 *
 * Index genome if needed, make bed file
 *
 */

 
process Prepare_Genome {
    tag                    { gf }
    executor               params.executor
    cpus                   params.cpus
    queue                  params.Queue
    memory                 "${params.mem} GB"
    module                 params.samtoolsMod,params.picardMod 
    storeDir               genomeStore
    
    input:
    file gf from genomeFile

    output:
    file "*.bed" into refBED
    file "*" into ref4VarScan
    

    script:
    """
    samtools faidx $gf
    cat ${gf}.fai | awk -v OFS='\t' '{chr = \$1; len = \$2; print chr, 1, len }' > ${gf}.bed
    java ${params.JavaOptions} -jar \${EBROOTPICARD}/picard.jar CreateSequenceDictionary R=${gf} O=${gf}.dict
    """
}


process Bowtie2_Index_Genome {
    tag                    { gf }
    executor               params.executor
    cpus                   params.alignThreads
    queue                  params.Queue
    memory                 "${params.alignMem} GB"
    module                 params.bowtie2Mod
    storeDir               genomeStore
    validExitStatus        0

    input:
    file gf from genomeFile

    output:
    file "*.bt2" into bowtie2Index

    script:
    """
    bowtie2-build --threads ${task.cpus} ${gf} ${gf.getBaseName()}
    """
}



/*
 *
 * Phase 2: Read preprocessing/QC
 *
 * fastp  performs QC_PreTrim, Trim+Filter, QC-PostTrim
 *
 */

Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_trimming; reads2Trim2 }


process Read_preparation {
    tag                    { name }
    executor               params.executor
    cpus                   params.trimThreads
    queue                  params.Queue
    memory                 "${params.trimMem} GB"
    module                 params.fastpMod
    publishDir             "${params.outdir}/Prepared-Reads", mode: 'copy'
    validExitStatus        0,1

    input:
    set val(name), file(reads) from read_files_trimming

    output:
    set val(name), file('*.trimmed.fq') optional true into trimmed2Align
    set val(name), file('*.json') optional true into  fastpstat
    file '*'
    
    script: 	    
    if(params.singleEnd){
	"""
	fastp --in1 ${reads[0]} --out1 "${name}.R1.trimmed.fq" \
	${adapterOptionsSE} ${params.trimOptions} ${trimOptions2} ${trimOptions3} \
	--thread ${task.cpus} -w ${task.cpus} --html "${name}"_fastp.html --json "${name}"_fastp.json

	"""
    } else {
	"""
	fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${name}.R1.trimmed.fq" --out2 "${name}.R2.trimmed.fq"  \
	${adapterOptionsPE}  ${params.trimOptions} ${trimOptions2} ${trimOptions3} \
	 --thread ${task.cpus} --html "${name}"_fastp.html --json "${name}"_fastp.json
	"""
    }  /* end if  singleEnd*/
    
}


/*
 *
 * Phase 3: BAM preparation
 *
 * Bowtie2_alignment, sort, dedup, filter
 *
 */


process Bowtie2_aln_n_sort {
    tag                    { id }
    executor               params.executor
    cpus                   params.alignThreads
    queue                  params.Queue
    memory                 "${params.alignMem} GB"
    errorStrategy          'finish'
    module                 params.bowtie2Mod,params.samtoolsMod,params.picardMod
    publishDir             "${params.outdir}/Bowtie2-align-n-sort", mode: 'copy'

    input:
    set val(id), file(reads) from trimmed2Align
    file idx from bowtie2Index

    output:
    set val(id), file('*.sorted.bam') optional true into bowtie2SortedBam
    file "*.stat"  optional true into alnstat      
    file '*.bai'   optional true into bowtie2SortedBai
    file '*'

    script:
    if(params.singleEnd){    
    """
    bowtie2 -p ${task.cpus} \
        ${params.bowtie2params} \
        -x ${genomeFile.getParent()}/${genomeFile.getBaseName()} \
        -1 ${reads[0]}  2> ${id}.aln.runlog \
        | samtools view -bhS - > ${id}.aln.bam       

    java ${params.JavaOptions} -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I=${id}.aln.bam \
      O=${id}.aln.wrg.bam \
      RGID=${id} \
      RGSM=${id} \
      RGLB=${params.RGLB} \
      RGPL=${params.RGPL} \
      RGPU=${params.RGPU} \
      RGCN=${params.RGCN}
     
    java ${params.JavaOptions} -jar \$EBROOTPICARD/picard.jar SortSam \
        SO=coordinate \
        I=${id}.aln.wrg.bam \
        O=${id}.sorted.bam \
        VALIDATION_STRINGENCY=LENIENT \
        CREATE_INDEX=true

    samtools flagstat ${id}.sorted.bam > ${id}.sorted.bam.stat         
    """
    } else {
    """
    bowtie2 -p ${task.cpus} \
        ${params.bowtie2params} \
        -x ${genomeFile.getParent()}/${genomeFile.getBaseName()} \
        -1 ${reads[0]} -2 ${reads[1]}  2> ${id}.aln.runlog \
        | samtools view -bhS - > ${id}.aln.bam       

    java ${params.JavaOptions} -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I=${id}.aln.bam \
      O=${id}.aln.wrg.bam \
      RGID=${id} \
      RGSM=${id} \
      RGLB=${params.RGLB} \
      RGPL=${params.RGPL} \
      RGPU=${params.RGPU} \
      RGCN=${params.RGCN}
     
    java ${params.JavaOptions} -jar \$EBROOTPICARD/picard.jar SortSam \
        SO=coordinate \
        I=${id}.aln.wrg.bam \
        O=${id}.sorted.bam \
        VALIDATION_STRINGENCY=LENIENT \
        CREATE_INDEX=true

    samtools flagstat ${id}.sorted.bam > ${id}.sorted.bam.stat         
    """
    } /* end if */
}

process MarkDupOnly {
    tag                    { id }
    executor               params.executor
    cpus                   params.dedupThreads
    queue                  params.Queue
    memory                 "${params.dedupMem} GB"
    module                 params.picardMod,params.samtoolsMod
    publishDir             "${params.outdir}/MarkDup", mode: 'copy'

    input:
    set val(id), file(sortedBam) from bowtie2SortedBam
    file sortedBai               from bowtie2SortedBai    

    output:
    set val(id), file('*.dedup.bam') optional true into bowtie2dupBam
    file "*.stat"  optional true into dedupstat     
    file '*.dedup.bai'  optional true into bowtie2dupBai
    file '*'


    script:
    """
    java  ${params.JavaOptions} -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        I=$sortedBam \
        O=${id}.dedup.bam \
        M=${id}.metrics.txt \
        CREATE_INDEX=true
        
    samtools flagstat ${id}.dedup.bam > ${id}.dedup.bam.stat     
    """
}


process RemoveDup_and_Filter {
    tag                    { id }
    executor               params.executor
    cpus                   params.dedupThreads
    queue                  params.Queue
    memory                 "${params.dedupMem} GB"
    module                 params.samtoolsMod
    publishDir             "${params.outdir}/ReadyBAMs", mode: 'copy'

    input:
    set val(id), file(dedupBam) from bowtie2dupBam
    file dedupBai from bowtie2dupBai
    
    output:
    set val(id), file('*.ready.bam') optional true into  readyBam1, readyBam2, readyBam3, readyBam4
    file '*.ready*bai'  optional true into  readyBai1, readyBai2, readyBai3, readyBai4  
    file "*.stat"  optional true into filterstat 
   
    file '*'

    script:
    """
    samtools view -bh -q${params.min_MQ} -F 1024 ${dedupBam} > ${id}.ready.bam
    samtools index ${id}.ready.bam   
    samtools flagstat ${id}.ready.bam > ${id}.ready.bam.stat
    """
}

process MultiQC {
    executor         params.executor
    cpus             2
    queue            params.Queue
    memory           "${params.trimMem} GB"
    module           params.multiqcMod
    publishDir       "${params.outdir}/MultiQC", mode: 'copy', overwrite: true

    input:
    file('/readprep/*') from fastpstat.collect()
    file('/aln/*')      from alnstat.flatten().toList()
    file('/dedup/*')    from dedupstat.flatten().toList()
    file('/filter/*')   from filterstat.flatten().toList()

    output:
    file "*_report.html" into multiqc_report_post
    file "*_data"

    script:
    """
    multiqc -d -f  .
    """
}

/*
 *
 * Phase 4: Variant discovery and annotation
 *
 * deepSNV variants, VarScan variants, SnpEff annotation
 *
 */


process deepSNV_RawVariants {
	tag { test_id }
	executor         params.executor
	cpus             params.RThreads
	queue            params.Queue
	memory           "${params.RMem} GB"
	module           params.deepSNV_RMod,params.PythonMod,params.BiopythonMod,params.MuscleMod 
	publishDir       "${params.outdir}/deepSNV_variants/${test_id}", mode: 'copy', overwrite: true
	validExitStatus  0,1
	errorStrategy   'finish'

	input:
	file refgenome from genomeFile
	file readybam from readyBam1.collect()
	file readybai from readyBai1.collect()
	set  test_id, control_id, optionsfile from deepSNV_param

	output:
	set val(test_id), file("${test_id}.rawVariantsOnly.csv")  optional true into VariantsCSVChannel
	file '*' 

	script:
	def input = "${test_id}.ready.bam"
	def control = "${control_id}.ready.bam"
	def optionf = "${optionsfile}"

	"""

	echo part1 run deepSNV to generate raw variants

	Rscript --vanilla --slave ${params.scriptdir}/deepSNV_R $refgenome \
	$input \
	$control \
	${params.method} \
	${params.P_CUT} \
	${params.combinePvalMethod} \
	${params.dispersion} \
	${params.stringentFreq} \
	${test_id}.rawVariantsOnly.csv \
	${test_id}.output.fa \
	${params.R_lib}

	"""

}  /* end process */

if (!params.skipAnnotation) {

	process VarScan_SnpEff_variants {
		tag { id }
		executor         params.executor
		cpus             params.cpus
		queue            params.Queue
		memory           "${params.mem} GB"
		module           params.VarScanMod,params.samtoolsMod
		publishDir       "${params.outdir}/VarScan_SnpEff_variants", mode: 'copy', overwrite: true
		validExitStatus  0,1
		errorStrategy   'finish'

		input:
		file refgenome from genomeFile
		file refidx from ref4VarScan
		set  val(id), file(readybam) from readyBam2
		file readybai from readyBai2
		file csvFile from VariantsCSVChannel.collect()

		output:
		set val(id), file("*.final.annotated.vcf")  optional true into CombineVariantsVCF
		file '*' 

		script:
		"""
		samtools mpileup -d ${params.maxDepth} -C ${params.min_MQ} -B -f ${refgenome} ${readybam} | java ${params.JavaOptions} -jar \${EBROOTVARSCAN}/VarScan.v2.3.9.jar mpileup2cns  \
		--p-value ${params.P_CUT} --min-var-freq ${params.stringentFreq} --min-coverage ${params.minCoverage} \
		--variants --output-vcf > ${id}.rawVariantsOnly.vcf

		java ${params.JavaOptions} -jar ${params.snpEffPath}/snpEff.jar \
		-config ${params.snpEffPath}/snpEff.config \
		-s ${id}-snpEff-summary-report.html \
		${params.snpEffDB} ${id}.rawVariantsOnly.vcf > ${id}.raw.snpEff.annotated.vcf


		${params.scriptdir}/common_variants.pl -invcf ${id}.raw.snpEff.annotated.vcf -incsv ${id}.rawVariantsOnly.csv -out ${id}.final.annotated.vcf
		"""	
	}  /* end process */	


} /* end if */

workflow.onComplete {

      def subject = "[deepSNV pipeline] Successful: $workflow.runName"
      if(!workflow.success){
          subject = "[deepSNV pipeline] FAILED: $workflow.runName"
      }
      
    finalLog = """
-------------------------------------------------------------------------    
[deepSNV Pipeline part 1] execution summary
-------------------------------------------------------------------------
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
workDir      : ${workflow.workDir}
exit status  : ${workflow.exitStatus}
Error report : ${workflow.errorReport ?: '-'}
-------------------------------------------------------------------------

"""

    ['mail', '-s', subject, params.email ].execute() << runInfo
    
    log.info "[deepSNV pipeline part1] COMPLETED. Sent summary e-mail to $params.email"

}

