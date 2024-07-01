#!/bin/bash
##### To count the number of reads in each lane file and compare the sum of these reads with that of the merged file

### sum the records in each .fastq.gz file from the lanes and compare the output with that from the merged files with the `if` command?

#### on the test data
# FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/ ## this rreads into the test folders's directories (output is the .gz files in each)
# ##### remove the record_counts directory if present to append data afresh each time a run the code
# rm -r /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/record_counts
# mkdir /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/record_counts
# echo created a record_counts directory
# cd $FOLDERS
# for folder in [0-9]*/;
# do 
#     echo $folder
#     prefix=$(echo $folder | grep -o "^[0-9]*")
#     ##### count the number of records in the fastq.gz file
#     echo $folder >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/record_counts/$prefix.record_count.txt
#     ### running the echo "grep..." to dry-run the code to test it
#     # echo "less $folder*| grep '^+$' | wc -l  >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/record_counts/$prefix.record_count.txt"
#     less $folder*| grep '^+$' | wc -l >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/record_counts/$prefix.record_count.txt
#     echo appended records in ${folder} into $prefix.record_count.txt
# done

#### on the actual data
FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/ ## this rreads into the test folders's directories (output is the .gz files in each)
##### remove the record_counts directory if present to append data afresh each time a run the code
rm -r /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/record_counts
mkdir /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/record_counts
echo created a record_counts directory
cd $FOLDERS
for folder in [0-9]*/;
do 
    echo $folder
    prefix=$(echo $folder | grep -o "^[0-9]*")
    ##### count the number of records in the fastq.gz file
    echo $folder >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/record_counts/$prefix.record_count.txt
    ### running the echo "grep..." to dry-run the code to test it
    # echo "less $folder*| grep '^+$' | wc -l  >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/record_counts/$prefix.record_count.txt"
    less $folder*| grep '^+$' | wc -l >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/record_counts/$prefix.record_count.txt
    echo appended records in ${folder} into $prefix.record_count.txt
done