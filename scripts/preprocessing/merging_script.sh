#!/bin/bash
#FILES=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/*
# cd /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/
# for file in *
# do
# echo $file merged_$file
# done


#### Started with an "if" statement . 
# for file in *
# do
#     sample=$(echo $file | cut -d "_" -f1)
#     # echo $sample $file
#     x=10
#     # echo $x
#     if [ "$sample" = "$x" ]; then
#         echo "yes"

#      else
#         echo "need to add 1"
#     fi
# x=$(expr $sample + 1)
# # echo $x
# done



# VAR1="Linuxize"
# VAR2="Linuxize"

# if [ "$VAR1" = "$VAR2" ]; then
#     echo "Strings are equal."
# else
#     echo "Strings are not equal."
# fi


### will try a "while" statement
# for file in *
# do
#     sample=$(echo $file | cut -d "_" -f1)
#     # echo $sample $file
#     x=10
#     # echo $x
#     while [ "$sample" = "$x" ]; 
#     do
#         echo "$file"
#         echo $x
    
# x=$(expr $sample + 1)
# echo $x
# done
# done

#cd /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/
#cd /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/data/raw_data/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078
# for file in *
# do
#     sample=$(echo $file | cut -d "_" -f1)  ##this was so as to get the number before the underscore, hence echo the file, then cut the output
# #     echo $sample $file
# #     x=1
# #     echo $x
# #     if [ "$sample" = "$x" ]; 
# #     then
# #         echo $sample
# #         echo "$file"
# #         # echo $x
# #     fi
# #     x=$(expr $sample + 1)
# #     x=$(expr $x + 1)
# # x=$((x+1))
# echo $sample|sort -n | uniq
# # echo "$file"
# done
# bash merging_script.sh | sort -n | uniq
#######################################################################################################################
### Let's give this code another try!
### This code now works

# FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/

# # ls $FOLDERS
# rm /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list.txt    #### need to learn how to catch exceptions on bash!!! 
# echo removed folder sample_list.txt if present
# touch /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list.txt
# for folder in `ls $FOLDERS`
# do 
#     # i=$(echo $folder | cut -d "_" -f1)

#     echo $folder | cut -d "_" -f1 >> /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list.txt
    
    
# done
# echo done creating a sample_list.txt
# # less sample_list.txt| sort -n| uniq > /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/trial_merged/sample_list.txt ### not sure why this did not work
# # less sample_list.txt|sort -n|uniq > /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/trial_merged/sample_list.txt  ## or this
# # sample_list_var=$(less sample_list.txt|sort -n|uniq)
# less /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list.txt| sort -n| uniq > /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list_sorted.txt ### this creates the list of sample no. 1 to 64
# echo done creating a sample_list_sorted.txt

# sample_list_sorted=$(less /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list_sorted.txt)
# for i in `less /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list_sorted.txt`  ### this and the code below do the same thing
# for i in $(less /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/sample_list_sorted.txt)
# for i in $sample_list_sorted
# do  
# #   echo $i

#     for folder in $FOLDERS
#     do
#         sample=$(ls $folder | cut -d "_" -f1)
      
#         echo $sample
#         # if [ "$sample" = "$i" ]
#         # then
#         #     echo$1
#         #     # echo $sample
#         # fi
#         # then
#     # echo zcat 
# done
# done
# echo done 

#### Trying Manaseh's approach
# FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/
# mkdir /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences
# echo created a merged_sequences directory
# for folder in $FOLDERS
# do 
#     echo $(ls $folder | grep -o "^[0-9]*")
#     prefix=$(ls $folder | grep -o "^[0-9]*")
#     # echo $prefix
#     ### running the echo "zcat..." to dry-run the code to test it
#     # echo "zcat $folder/* >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences/$prefix.merged_fastq.gz"
#     echo zcat "${folder}"/* >> "/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences/${prefix}.merged_fastq.gz"
#     zcat $folder/* >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences/$prefix.merged_fastq.gz

#     # echo merged fastq.gz files in ${folder}
# done

##### Trial II
#### Trying Manaseh's approach
# FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/ ## this rreads into the test folders's directories (output is the .gz files in each)
# # FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/ ## finds no folders matching the pattern [0-9]*

# mkdir /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences
# echo created a merged_sequences directory
# cd $FOLDERS
# for folder in [0-9]*/;
# do 
#     # echo $folder
#     prefix=$(echo $folder | grep -o "^[0-9]*")
#     # echo $prefix
#     ### running the echo "zcat..." to dry-run the code to test it
#     # echo "zcat $folder/* >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences/$prefix.merged_fastq.gz"
#     zcat $folder/* >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences/$prefix.merged_fastq.gz

#     echo merged fastq.gz files in ${folder}
# done
# ####+++ the above code works perfectly thus far!+++

####### Now modifying the script to run on the actual data
FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/ ## this rreads into the test folders's directories (output is the .gz files in each)
rm -r /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences
mkdir /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences
echo created a merged_sequences directory
cd $FOLDERS
for folder in [0-9]*/;
do 
    # echo $folder
    prefix=$(echo $folder | grep -o "^[0-9]*")
    # echo $prefix
    ### running the echo "zcat..." to dry-run the code to test it
    # echo "zcat $folder/* >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/test_data/merged_sequences/$prefix.merged.fastq"
    zcat $folder/* >> /mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences/$prefix.merged.fastq

    echo merged fastq files in ${folder}
done
echo done merging

FOLDERS=/mnt/f/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences ## this reads into the test folders's directories (output is the .gz files in each)
rm -r /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC
mkdir /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC
cd $FOLDERS
for file in *
do 
    # echo $file
    # echo "fastq $file -o /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC/"
    fastqc $file -o /mnt/d/MSc_Bioif/EANBIT_cohort4/MSc_project/MScProject/results/raw_reads_QC/
    echo FASTQC report for $file is ready for review

done
echo done running the quality check
