#! /usr/bin/env nextflow

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "    gatk4-DataPreProcessing-nf v1: From BAM to hg38 Callable BAM       "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/gatk4-DataPreProcessing-nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input                   BAM FILE              BAM file (between quotes for BAMs)"
    log.info "--output_dir             OUTPUT FOLDER         Output for gVCF file"
    log.info "--ref_fasta              FASTA FILE            Reference FASTA file"
    log.info "--dbsnp                  VCF FILE              dbSNP VCF file"
    log.info "--onekg                  VCF FILE              1000 Genomes High Confidence SNV VCF file"
    log.info "--mills                  VCF FILE              Mills and 1000 Genomes Gold Standard SID VCF file"
    log.info "--gatk_jar               JAR PATH              Full path to GATK4 binary"
    exit 1
}


//
// Parameters Init
//
params.input      = null
params.output_dir = "."
params.gatk_jar   = null
params.ref_fasta  = null
params.dbsnp      = null
params.onekg      = null
params.mills      = null


//
// Parse Input Parameters
//
ubam_ch   = Channel
			.fromPath(params.input)
			.map { input -> tuple(file.baseName, input) }
output    = file(params.output_dir)
GATK      = params.gatk_jar
ref       = file(params.ref_fasta)
ref_dbsnp = file(params.dbsnp)
ref_1kg   = file(params.onekg)
ref_mills = file(params.mills)
ref_dict  = ref.parent / ref.baseName + ".dict"
ref_in    = ref.parent / ref.baseName + ".fasta.fai"
ref_amb   = ref.parent / ref.baseName + ".fasta.amb"
ref_ann   = ref.parent / ref.baseName + ".fasta.ann"
ref_bwt   = ref.parent / ref.baseName + ".fasta.bwt"
ref_pac   = ref.parent / ref.baseName + ".fasta.pac"
ref_sa    = ref.parent / ref.baseName + ".fasta.sa"


//
// Process uBAM to FASTQ
//
process PICARD_SamToFastq{
  tag "$sampleID"

  cpus 8
  memory '16 GB'
  time '2h'
    
  input: 
      set sampleID, file(ubam) from ubam_ch

  output: 
      set sampleID, file("${sampleID}_1.fastq"), file("${sampleID}_2.fastq") into fastq_ch

	script:
	"""
	java -Xmx4g -XX:ParallelGCThreads=8 -jar \$PICARD_TOOLS_LIBDIR/picard.jar \
		SamToFastq \
		VALIDATION_STRINGENCY=LENIENT \
		I=${ubam} \
		F=${sampleID}_1.fastq \
		F2=${sampleID}_2.fastq
	
	"""
}


//
// Process Mapping With BWA MEM (AND MORE)
//
process BWA_mapping{
  tag "$sampleID"

  cpus 8
  memory '64 GB'
  time '12h'
  
  input: 
      file genome from ref
      file "${genome.baseName}.dict" from ref_dict
      file "${genome.baseName}.fasta.fai" from ref_in
      file "${genome.baseName}.fasta.amb" from ref_amb
      file "${genome.baseName}.fasta.ann" from ref_ann
      file "${genome.baseName}.fasta.bwt" from ref_bwt
      file "${genome.baseName}.fasta.pac" from ref_pac
      file "${genome.baseName}.fasta.sa" from ref_sa
	  set sampleID, file(fastq1), file(fastq2) from fastq_ch
  
  output: 
      set sampleID, file("${sampleID}.aln.bam"), file("${sampleID}.aln.bam.bai") into aln_bam_ch

  script:    
  """
	bwa mem -t 8 -R '@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:Illumina' $genome $fastq1 $fastq2 | \
		k8 \$BWAKIT_EXEDIR/bwa-postalt.js \$POSTALT_REF | \
		sambamba view -t 8 -S -f bam -l 0 /dev/stdin | \
		sambamba sort -t 8 -m 6G --tmpdir=/tmp /dev/stdin -o ${sampleID}.aln.bam		
	
  """

}


//
// Process Mark Duplicates
//
process SAMBAMBA_markdup{
  tag "$sampleID"

  cpus 8
  memory '64 GB'
  time '12h'
  
  module 'sambamba/0.6.5'
  
  input: 
	  set sampleID, file(bam), file(bai) from aln_bam_ch
  
  output: 
      set sampleID, file("${sampleID}.aln.dup.bam"), file("${sampleID}.aln.dup.bam.bai") into dup_bam_ch

  script:    
  """
	sambamba markdup -t 8 --tmpdir=/tmp $bam ${sampleID}.aln.dup.bam
	
  """

}


//
// Process BaseRecalibrator + ApplyBQSR
//
process GATK_BaseRecalibrator_ApplyBQSR{
  tag "$sampleID"

  cpus 1
  memory '32 GB'
  time '4h'
    
  publishDir "$output/$sampleID", mode: 'copy'
  
  input: 
      file genome from ref 
      file "${genome.baseName}.dict" from ref_dict
      file "${genome.baseName}.fasta.fai" from ref_in
      set sampleID, file(bam), file(index) from dup_bam_ch

  output: 
      set sampleID, file("${sampleID}.recal.bam"), file("${sampleID}.recal.bai") into recal_bam_ch

  script:    
  """
	java -XX:ParallelGCThreads=1 -Xmx4G -Djava.io.tmpdir=/tmp -jar ${GATK} \
    		BaseRecalibrator \
    		-R ${genome} \
    		-I ${bam} \
            --use-original-qualities \
            -O ${sampleID}.recal_data.table \
            -known-sites ${ref_dbsnp} \
            -known-sites ${ref_1kg} \
            -known-sites ${ref_mills}
			
	java -XX:ParallelGCThreads=1 -Xmx4G -Djava.io.tmpdir=/tmp -jar ${GATK} \
    		ApplyBQSR \
    		-R ${genome} \
    		-I ${bam} \
            --use-original-qualities \
            -O ${sampleID}.recal.bam \
            --bqsr-recal-file ${sampleID}.recal_data.table
			
  """

}


//
// Process Post-Alignment QC
//
process QUALIMAP_BamQC{
  tag "$sampleID"

  cpus 8
  memory '32 GB'
  time '4h'
      
  publishDir "$output/$sampleID", mode: 'copy'
  
  input: 
	  set sampleID, file(bam), file(bai) from recal_bam_ch

  output: 
      set sampleID, file("${sampleID}.recal_stats") into recal_bam_stats_ch

  script:    
  """
	qualimap bamqc \
		--java-mem-size=24G \
		-c \
		-nt 8 \
		-bam $bam

  """

}





