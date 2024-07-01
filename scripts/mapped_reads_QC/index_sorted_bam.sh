#!/bin/bash

## create 'mscproject' conda environment
# conda create --name mscproject -y

### to activate "mscproject" conda env to work in the project environment for reporducibility of my work
### to check if the environment is present first
env_name="mscproject"
if  ! conda info --env | grep "mscproject"
    then
        echo Error: the conda environment $env_name does not exits.
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
fi

# ## dry-run the 'source' command, then run the code
# # echo "source '$(conda info --base)/etc/profile.d/conda.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$env_name"

## to indexed the sorted bam files
SECONDS=0
## Run samtools to index .bam files 
external_drive_path_hisat2="/media/user/HildaWN/FASTQ_files/HISAT_alignment/results/aligned_reads/sorted_bam"

cd $external_drive_path_hisat2
# echo `pwd`

# for i in *.bam
# do
#     # echo $i
#     # echo "samtools index $i"
#     samtools index $i
# done


echo Finished sorting the sam files
duration=$SECONDS
echo "$duration" seconds
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
echo " "
echo "##############################################################################################################"

## compress the indexed .bam files 
# zip -r hisat_sorted_indexed_bam_files.zip *.bam*
if ! ls | grep -o "hisat_sorted_indexed_bam_files.zip"
    then
        echo "Indexed bam files not compressed"
        exit 1
    else
        echo "Indexed bam files successfully compressed!"
fi

echo " "
echo "##############################################################################################################"
