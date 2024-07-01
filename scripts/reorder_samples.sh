#!/usr/bin/bash

sample_list="/mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/sample_list.csv"
count_dir="/mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/HTSeq/"
# rm -r /mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/HTSeq_reordered_files/
mkdir /mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/reorderedHTSeq_files/
reordered_sampleFiles="/mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/reorderedHTSeq_files/"

# for line in `less $sample_list`
for line in `cut -d "," -f1 $sample_list`
do
    for file in `ls $count_dir`
do
    # echo $line
    # echo "running please wait"
    # echo $file
    prefix=`echo $file | grep -o "sample_[0-9]*"`
    # echo $prefix 
    echo $line : $prefix
    # echo `cut -d "," -f2 $file`
    # cut_f1=`cut -d "," -f1 $line `
    # echo $cut_f1
    # if [ $line != $prefix ] ;
    # then
    #     echo "not matching"
    #     echo $line : $prefix
    #     # `cp $file $reordered_sampleFiles`
    # continue

    # else 

    # echo $line : $prefix
    # echo "found a match"
    # cp "$count_dir$file" "$reordered_sampleFiles"
    
    # fi
    if [ $line == $prefix ] ;
    then
        echo "matching"
        echo $line : $prefix
        cp "$count_dir$file" "$reordered_sampleFiles"
    continue
    
    fi
done
done

# cp /mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/HTSeq/sample_5_counts.txt /mnt/c/hne_files/MSc_Bioinf/hisat2htseq_dataset/HTSeq_reordered_files/