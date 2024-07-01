#!/bin/bash

## create 'mscproject' conda environment
# conda create --name mscproject -y

### to activate "mscproject" conda env to work in the project environment for reporducibility of my work
### to check if the environment is present first
env_name="mscproject"
if  ! conda info --env | grep "mscproject"
    then
        echo Error: the conda environment $env_name does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
fi

# ## dry-run the 'source' command, then run the code
# # echo "source '$(conda info --base)/etc/profile.d/conda.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$env_name"

SECONDS=0
## Run samtools to convert the .sam files to .bam
external_drive_path_samout="/mnt/d/FASTQ_files/HTSeq/counts_samout_intersection-empty"
mkdir $external_drive_path_samout/bamout

cd $external_drive_path_samout/bamout
# echo `pwd`

for i in `ls ../*.sam`
do
    # echo $i
    prefix=`echo $i | grep -o "sample_[0-9]*"`
    # echo $prefix
    # echo "samtools sort -O BAM -o "$prefix"_sorted.bam $i"
    samtools view -O BAM -o "$prefix"_bamout.bam $i
done

echo Finished sorting the sam files
duration=$SECONDS
echo "$duration" seconds
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
echo " "
echo "##############################################################################################################"