#!/bin/bash
#SBATCH --job-name="bwa/samtools"
#SBATCH --time=15-00:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem=385092M
#SBATCH --cpus-per-task=12
#SBATCH --partition=Orion

###-------------------------------------------------
# script for bwa alignment of paired end fastq files
###-------------------------------------------------

module load bwa
module load samtools

### Set up path shortcuts
fastq_input=/nobackup/cooper_research/Whitney/Culex_WGS_TrimmedFiltered_Renamed

output_dir=/scratch/wbranne1/Outputs/BWA_culex

ref_genome=/projects/cooper_research/Ref_Genomes/Mosquito/culex_assembly.fasta

rg_sorted=/scratch/wbranne1/Outputs/sorted_with_readgroups

for file in $fastq_input/*_1.fq.gz; do
    ### establish name variables for each file set
    full_name=$(basename "$file")
    read_pair=$(echo $full_name | cut -d _ -f 1-5)
    sample_id=$(echo $read_pair | cut -d _ -f 2)
    pop_id=${sample_id::-1}
    
    ### run BWA and pipe to samtools sort bam file
    bwa mem -t 16 $ref_genome $fastq_input/${read_pair}_1.fq.gz $fastq_input/${read_pair}_2.fq.gz \
    | samtools sort -@16 -o $output_dir/${read_pair}.sorted.bam -
 
    ### add read groups
    samtools addreplacerg -r "@RG\tID:${read_pair}\tSM:${sample_id}" \
    -o $rg_sorted/${read_pair}.rg.bam $output_dir/${read_pair}.sorted.bam

    ### keep track of progress in slurm out file
    echo $read_pair completed.

done
    
