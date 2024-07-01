#/usr/bin/bash

### to activate "htseq" conda env to work in the project environment where the hisat2 tool is installed. Also for reporducibility of my work
### to check if the environment is first present
env_name="htseq"
if  ! conda info --env | grep "htseq"
    then
        echo Error: the conda environment $env_name does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
fi

## dry-run the 'source' command, then run the code
# echo "source '$(conda info --base)/etc/profile.d/conda.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$env_name"

SECONDS=0
echo  "Performing read quantification using HTSeq - read_quant_htseq_samout_union.sh"

### PATH TO THE... ANNOTATION FILE <-- INPUT FILE
# external_drive_path_gff="/media/user/HildaWN/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27"
external_drive_path_gtf="/mnt/d/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27"
external_drive_path_gtf="/mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27"

### ... SORTED .BAM FILES <-- INPUT FILES
# external_drive_path_sorted_bam="/media/user/HildaWN/FASTQ_files/HISAT_alignment/results/aligned_reads/sorted_bam"
# external_drive_path_sorted_bam="/mnt/d/htseq_trials"
external_drive_path_sorted_bam="/mnt/d/FASTQ_files/HISAT_alignment/results/aligned_reads/sorted_bam/"
external_drive_path_sorted_bam="/mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/hisat_sorted_indexed_bam"
cd $external_drive_path_sorted_bam

### ... HTSEQ OUTPUT FILES
# external_drive_path_htseq="/media/user/HildaWN/FASTQ_files/HTSeq"
external_drive_path_htseq="/mnt/d/FASTQ_files/HTSeq/counts_samout_union"
external_drive_path_htseq="/mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/HTSeq_rerun"


### ... HTSEQ SAMOUT FILES
external_drive_path_samout="/mnt/d/FASTQ_files/HTSeq/counts_samout_union"
external_drive_path_samout="/mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/samout_output"

## (E) quantifying reads while writing out all SAM alignments again, this time with an addition field with the tag "XF"
# that identifies each alignment with its assigned feature 
## NB (B) is what was originally ran on the data.


for bam_file in *.bam
do 
    prefix=`echo $bam_file | grep -o "sample_[0-9]*"`
    echo $bam_file
    echo $prefix
    htseq-count -f bam \
    -s yes \
    -t exon \
    -i gene_id \
    --additional-attr=gene_name \
    -m union \
    --nonunique none \
    --version \ 
    --samout=$external_drive_path_samout/"$prefix"_samout.sam \
    $bam_file $external_drive_path_gtf/genomic.gtf > $external_drive_path_htseq/"$prefix"_counts.txt

    samtools 
done

echo " "
echo Finished quantifying reads
duration=$SECONDS
echo "$duration" seconds
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
echo " "
echo "##############################################################################################################"

# htseq-count -f bam -s yes -t exon -i gene_id --additional-attr=gene_name -m union --nonunique none --samout=sample_1_samout.sam ../hisat_sorted_indexed_bam/sample_1_sorted.bam /mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27/genomic.gtf > sample_1_counts.txt
