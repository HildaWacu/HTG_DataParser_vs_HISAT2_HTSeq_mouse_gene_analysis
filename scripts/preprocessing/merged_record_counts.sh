#!/bin/bash
##### count the records(reads) for each of the merged fastq.gz files 
FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences ## this reads into the test folders's directories (output is the .gz files in each)
##### remove the record_counts directory if present to append data afresh each time a run the code
rm -r /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_record_counts
mkdir /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_record_counts
echo created a merged_record_counts directory
cd $FOLDERS
for folder in [0-9]*;
do 
    echo $folder
    ##### count the number of records(reads) in the merged_fastq.gz file
    echo $folder >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_record_counts/total_sample_reads.txt
    # echo "less $folder| grep '^+$' | wc -l  >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_record_counts/total_sample_reads.txt"
    less $folder| grep '^+$' | wc -l  >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_record_counts/total_sample_reads.txt
    echo appended reads in ${folder} to total_sample_reads.txt
done