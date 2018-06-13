# gatk4-DataPreProcessing-nf
Nextflow pipeline for pre-process BAM(s) with hg38 and GATK4, following GATK [Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145).

<div style="text-align:center"><img src="https://us.v-cdn.net/5019796/uploads/editor/3o/dznasg7toiq1.png" width="200" /></div>


## Description

Tailored to fit the need of re-analyzing BAM files under new GATK4 Best Practices, and with all hg38 databases.

## Dependencies 

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. [GATK4 executables](https://software.broadinstitute.org/gatk/download/)
3. [Picard tools](https://broadinstitute.github.io/picard/)
4. [BWA](https://github.com/lh3/bwa/tree/master/bwakit), especially BWAKIT, because a post-alignment treatment is required ([more info](https://github.com/lh3/bwa/blob/master/bwakit/bwa-postalt.js)).
5. [Sambamba](http://lomereiter.github.io/sambamba/).
6. [Qualimap](http://qualimap.bioinfo.cipf.es/) binary in your PATH (for a nice QC per BAM).
7. References (genome in fasta, dbSNP vcf, 1000 Genomes vcf, Mills and 1000 Genomes Gold Standard vcf), available in [GATK Bundle](https://software.broadinstitute.org/gatk/download/bundle).

## Input

- `--input` : your intput BAM file(s) (do not forget the quotes for multiple BAM files e.g. `--input "test_*.bam"`)
- `--output_dir` : the folder that will contain your aligned, recalibrated, analysis-ready BAM file(s).
- `--ref_fasta` : your reference in FASTA. 
- `--dbsnp` : dbSNP VCF file. 
- `--onekg` : 1000 Genomes High Confidence SNV VCF file. 
- `--mills` : Mills and 1000 Genomes Gold Standard SID VCF file. 
- `--gatk_jar` : the full path to your GATK4 jar file.

A nextflow.config is also included, please modify it for suitability outside our pre-configured clusters ([see Nexflow configuration](https://www.nextflow.io/docs/latest/config.html#configuration-file)).

## Usage for Cobalt cluster
```
nextflow run iarcbioinfo/gatk4-DataPreProcessing.nf -profile cobalt --input "/data/test_*.bam" --output_dir /data/myRecalBAMs --ref_fasta /ref/Homo_sapiens_assembly38.fasta --gatk_jar /bin/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar --dbsnp /ref/dbsnp_146.hg38.vcf.gz --onekg /ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz --mills Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

