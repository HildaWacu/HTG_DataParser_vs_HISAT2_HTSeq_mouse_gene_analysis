#!/bin/bash

## create 'htseq' conda environment
# conda create --name htseq -y

### to activate "htseq" conda env to work in the project environment for reporducibility of my work
### to check if the environment is present first
env_name="htseq"
if  ! conda info --env | grep "htseq"
    then
        echo Error: the conda environment $env_name does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
fi

# ## dry-run the 'source' command, then run the code
# # echo "source '$(conda info --base)/etc/profile.d/conda.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$env_name"

## install HTSeq from bioconda
echo "Installing HTSeq"
conda install -c bioconda htseq -y

## check if HTSeq has been succesfully installed
if ! conda list | grep "htseq"
then 
    echo Error: HTSeq package does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
    else 
        echo " "
        echo "HTSeq successfully installed"
fi
        

echo " "
echo "################################################################################################################################################################"