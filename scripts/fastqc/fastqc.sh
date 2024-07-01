#!/bin/bash

# to activate "mscproject" conda env to work in the project environment for reporducibility of my work
### to check if the environment is present first
env_name="mscproject"
if ! conda info --env | grep "mscproject"
    then
        echo Error: the conda environment $env_name does not exits.
        # echo $?
        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
fi
## dry-run the 'source' command, then run the code
# echo "source '$(conda info --base)/etc/profile.d/conda.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$env_name"

# FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences ## this reads into the test folders's directories (output is the .gz files in each)
# rm -r /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC
# mkdir /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC
# cd $FOLDERS
# for file in *
# do 
#     # echo $file
#     # echo "fastq $file -o /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC/"
#     fastqc $file -o /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC/
#     echo FASTQC report for $file is ready for review

# done
# echo done running the quality check

# to run multiqc on the fastqc results
# rm -r /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/multiQC
# mkdir /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/multiQC
# multiqc /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC -o /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/multiQC

echo done aggregating fastqc reports to a single MultiQC report
## confirms that indeed the script is run in the mscproject conda env
conda info --env