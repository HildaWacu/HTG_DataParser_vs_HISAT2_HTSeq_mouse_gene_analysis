#!/bin/bash

path=/mnt/c/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts
list1=`cut -d "," ${path}/exclusive_HTG_HTSeq_filtered.csv -f2`
list2=`cut -d "," ${path}/exclusive_HTG_HTSeq_filtered.csv -f3`

for string1 in ${list1};
do 
for string2 in ${list2};
    do
    # echo "${string1:0:4},${string2}"

    if [ "${string1:0:4}" = "${string2:0:4}" ]; then
    echo "${string1},${string2}" 
fi
    done
done > ${path}/present_in_both_P_3char.csv
# done > ${path}/trial_file.txt

for string1 in ${list1};
do 
for string2 in ${list2};
    do
    # echo "${string1:0:4},${string2}"

    if [ "${string1:0:5}" = "${string2:0:5}" ]; then
    echo "${string1},${string2}" 
fi
    done
done > ${path}/present_in_both_P_4char.csv
# done >> ${path}/trial_file.txt


# string1="H2-L_H2-D1"
# string2="H2-T-ps"
# if [ "${string1:0:4}" = "${string2:0:4}" ]; then
#     echo "yes"
# else
#     echo "no"
# fi



