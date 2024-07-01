#!/bin/bash

# # create 'mscproject' conda environment
# conda create --name mscproject -y

## to activate "mscproject" conda env to work in the project environment for reporducibility of my work
## to check if the environment is present first
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

## install HISAT2 from Bioconda in the mscproject
conda install -c bioconda hisat2 -y

## check if HISAT2 has been successfully installed
if  ! conda list | grep "hisat2"
    then
        echo Error: hisat2 package does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
    else 
        echo HISAT2 successfully installed
fi

echo " "
echo "################################################################################################################################################################"