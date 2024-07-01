#!/bin/bash

# to work in the project environment fo ensure reporducibility of my work
# conda activate mscproject ### find out how to add a script to activate conda activate
#to unzip the fastqc reports inthe .zip format
FOLDERS=/mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC
cd $FOLDERS

for file in *.zip
do 
    #to dry-run the code
    echo "unzip $file"
    # unzip $file
    echo $file unzipped
done
echo done unipping the fastqc.zip files

