#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=merge/call
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=3
#SBATCH --mem=385092M
#SBATCH --time=15-00:00:00
#SBATCH --mail-user=wbranne1@charlotte.edu
#SBATCH --mail-type=END,FAIL

module load samtools

### establish file path locations

ref_genome="/projects/cooper_research2/Whitney/RefGenome/Culex-tarsalis_knwr_CONTIGS_CtarK1.fa"
bam_output="/projects/cooper_research2/Whitney/WGS_bam"

### merge bams 
samtools merge --threads 16 $bam_output/merged_WGS.bam $bam_output/*.rg.bam 

echo Progress: Merge done.

### sort merged file by query names for fixmate
samtools sort -o $bam_output/name_sorted.bam -n --threads 16 $bam_output/merged_WGS.bam 
### fixmate adds flags for markdup
samtools fixmate --threads 16 -m $bam_output/name_sorted.bam $bam_output/fixmate.bam 
### resort by coordinates for markdup
samtools sort -o $bam_output/positionsort.bam $bam_output/fixmate.bam 
### remove duplicates
samtools markdup --threads 16 -r -f $bam_output/markdups_info.txt $bam_output/positionsort.bam $bam_output/culex_WGS_rm_dups.bam 

echo Progress: Duplicates removed.

## create index file 
samtools index -b $bam_output/merged_WGS.bam
echo Progress: File indexed.

### bam -> bcf | call snps and indels (-v) multiallelic (-m) to vcf.gz file (-Ov)
bcftools mpileup -Ou -f $ref_genome $bam_output/merged_WGS.bam \
	| bcftools call -Oz -vm -o $bam_output/culex_WGS_variants.vcf.gz -

echo Progress: VCF file completed.