#!/usr/bin/bash


samout_files=/mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/samout_output_sample
cd $samout_files
# echo $samout_files
# less -S sample_1_samout.sam |rev | cut -f1 |rev| uniq |head
# `less -S $samout_files |rev | cut -f1 |rev| uniq |head`
# for samout_file in *.sam
# do
#     genename=`cat "$samout_file" | rev | cut -f1 | rev | uniq | cut -d ":" -f3`
#     echo "$genename"
#     sequence=`cut $samout_file -f10| uniq`
#     # if 

# done

htseq_genelist_output=/mnt/c/hne_files/MSc_Bioinf/rerunning_HTSeq/htseq_genelist_output

for samout_file in *.sam
do
    genename=`cat "$samout_file" | rev | cut -f1 | rev | cut -d ":" -f3`
    # echo "$genename"
    sequence=`cut $samout_file -f10`
    # echo "$genename,$sequence" > $htseq_genelist_output/sample_1_htseq_genelist.csv
    for i in $genename
    do
    for j in $sequence
    do 
        echo $i,$j >> $htseq_genelist_output/sample_1_htseq_genelist.txt
done 
done
done

less -S $htseq_genelist_output/htseq_genelist.txt |sort |uniq > htseq_final_genelist.csv
