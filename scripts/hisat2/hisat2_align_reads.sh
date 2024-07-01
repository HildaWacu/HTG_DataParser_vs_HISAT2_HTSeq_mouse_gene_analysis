#!/bin/bash

### to activate "mscproject" conda env to work in the project environment where the hisat2 tool is installed. Also for reporducibility of my work
### to check if the environment is present first
env_name="mscproject"
if  ! conda info --env | grep "mscproject"
    then
        echo Error: the conda environment $env_name does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
fi

## dry-run the 'source' command, then run the code
# echo "source '$(conda info --base)/etc/profile.d/conda.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$env_name"

echo  "Running the alignment with indexed ref gen GRCm39 indexed genome"
# path_to_dir2="/Volumes/HildaWN/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences"
external_drive_path_hisat2="/media/user/HildaWN/FASTQ_files/HISAT_alignment"
external_drive_path_fastq="/media/user/HildaWN/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences"

## create a dir for the output
mkdir -p $external_drive_path_hisat2/results/aligned_reads $external_drive_path_hisat2/results/summary_files

## run the alignment

cd $external_drive_path_fastq

SECONDS=0

for i in *.fastq
    do 
        #echo $i
        prefix=$(echo $i | grep -o "^[0-9]*")
        # echo $i
        ## dry run the conmfirm the paths are correct
        # echo $external_drive_path_hisat2/results/sample_$prefix.sam --summary-file $external_drive_path_hisat2/results/summary_sample$prefix.txt

        hisat2 -q -x $external_drive_path_hisat2/ref_genome_grcm39/indexed_ref/genome \
        -U $i \
        -S $external_drive_path_hisat2/results/aligned_reads/sample_$prefix.sam \
        --summary-file $external_drive_path_hisat2/results/summary_files/summary_sample$prefix.txt

done

echo " "
echo " Finished aligning the reads"
duration=$SECONDS
echo $duration
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"

echo " "
echo "###############################################"



