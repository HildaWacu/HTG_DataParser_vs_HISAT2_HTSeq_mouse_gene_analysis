#!/usr/bin/bash

### to activate conda env to work in the project environment where the hisat2 tool is installed. Also for reporducibility of my work

# env_name="htseq"
env_name=$1
echo creating conda environment $env_name

if  ! conda info --env | grep "$env_name"
    then
        #source "$(conda info --base)/etc/profile.d/conda.sh"
        conda create -n "$env_name" -y
        # echo Conda environment "$env_name" created and activated successfully!
    else
        echo Error: the conda environment "$env_name" already exits.
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
        # echo $?        
fi

if conda info --env | grep "$env_name"
    then 
        echo Conda environment "$env_name" created and activated successfully!
    else
        echo Error: the conda environment "$env_name" not successfuly created.
        exit 1
fi


## https://www.baeldung.com/linux/use-command-line-arguments-in-bash-script