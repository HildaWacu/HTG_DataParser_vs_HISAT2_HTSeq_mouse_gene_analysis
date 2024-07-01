#!/bin/bash

# ## create 'mscproject' conda environment
# conda create --name mscproject -y

# # to activate "mscproject" conda env to work in the project environment for reporducibility of my work
# ### to check if the environment is present first
# env_name="mscproject"
# if  ! conda info --env | grep "mscproject"
#     then
#         echo Error: the conda environment $env_name does not exits.
#         # echo $?
#         exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
# fi

# ## dry-run the 'source' command, then run the code
# # echo "source '$(conda info --base)/etc/profile.d/conda.sh"
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate "$env_name"

# ## install HISAT2 from Bioconda in the mscproject
# conda install -c bioconda hisat2 -y

# ## check if HISAT2 has been successfully installed
# if  ! conda list | grep "hisat2"
#     then
#         echo Error: hisat2 package does not exits.
#         # echo $?
#         exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
#     else 
#         echo HISAT2 successfully installed
# fi

#echo " "
#echo "################################################################################################################################################################"


## dowmload the reference genome
### first create the ref  dir
#mkdir -p /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome
#cd /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome

## ref genome
# curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001635.27.zip" -H "Accept: application/zip"
#curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001635.27.zip" -H "Accept: application/zip"
#unzip GCF_000001635.27.zip

## to install datasets package
#conda install -c conda-forge datasets -y

## check if datasets has been successfully installed
#if  ! conda list | grep "datasets"
#    then
#        echo Error: datasets package does not exits.
#        # echo $?
#        exit 1 ## once it encounters a situation where this condition is met, the Bash script exits and commands cease to be executed
#    else 
#        echo 'datasets' successfully installed
#fi

## Download the ref genome dataset
# datasets download genome accession GCF_000001635.27 --include gff3,rna,cds,protein,genome,seq-report --filename GCF_000001635.27.zip

#datasets download genome accession GCF_000001635.27 --include genome,seq-report --filename GCF_000001635.27.zip
#echo " "
#echo "########################################################################################################################################################################################3"

# ## confirm is ref file has been successfully downloaded
# if [ -z "$(ls -A . )" ]; then
#    echo "Error: ref genome not downloaded"
#    exit 1
# else
#    echo "Ref genome downloaded successfully!"
# fi

# ## build an index of the reference
# mkdir /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome/indexed_ref

# hisat2-build /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome/indexed_ref/genome  # where genome is the basename
#mkdir indexed_ref
### rename the .fna files, update the direcotries to the data .fna files
path_to_dir="/Users/eht/Downloads/hilda_ref/HISAT_alignment"
#echo $path_to_dir/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna 
#cp $path_to_dir/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna  $path_to_dir/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fa
#hisat2-build $path_to_dir/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fa /Users/eht/Downloads/hilda_ref/HISAT_alignment/ref_genome_grcm39/indexed_ref/genome 

# if [ -z "$(ls -A /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome/indexed_ref)" ]; then
#    echo "Error: ref indexes not built"
#    exit 1
# else
#    echo "Buiding of ref indexes complete"
# fi


# ## perform alignment of the reads and output to a SAM file
# mkdir /media/user/HildaWN/FASTQ_files/HISAT_alignment/results

# hisat2 -q -x /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome/indexed_ref/genome -U /media/user/HildaWN/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences/*.merged.fastq -S /media/user/HildaWN/FASTQ_file/HISAT_alignment/results/hisat2_align.sam 
# if [ -z "$(ls -A /media/user/HildaWN/FASTQ_files/HISAT_alignment/results)" ]; then
#    echo "Error: alignment incomplete"
#    exit 1
# else
#    echo "HISAT2 alignment complete!"
# fi  

echo " "
echo "########################################################################################################################################################################################3"
#echo  "Running the alignment with indexed ref gen GRCm38 indexed genome"
#zSECONDS=0
## STEP 1 - download the indexed ref from the HISAT2 website
## create the working directory and sub-directories
#mkdir -p /Users/eht/Downloads/hilda_ref/HISAT_alignment/results
#cd /Users/eht/Downloads/hilda_ref/HISAT_alignment

#wget "https://cloud.biohpc.swmed.edu/index.php/s/grcm38/download"


echo  "Running the alignment with indexed ref gen GRCm39 indexed genome"
SECONDS=0
path_to_dir2="/Volumes/HildaWN/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences"
cd $path_to_dir2
#astq_files=$(ls $path_to_dir2)
#echo $fastq_files

## creating an array/list of my fastq files
#fastq_files_list= { $(ls $path_to_dir2) }

#echo $fastq_files_list
#echo $path_to_dir2
#hisat2 -q -x /Users/eht/Downloads/hilda_ref/indexed_ref/grcm38/genome -U /Volumes/HildaWN/FASTQ_files/210616_201118_POk_Hultgren_EdgeSeq-13320318/FASTQ_Generation_2021-06-16_19_46_15Z-19096078/merged_sequences/*.merged.fastq -S /Users/eht/Downloads/hilda_ref/HISAT_alignment/results/hisat2_align.sam --summary-file /Users/eht/Downloads/hilda_ref/HISAT_alignment/results/
#hisat2 -q -x $path_to_dir/ref_genome_grcm39/indexed_ref/genome -U $fastq_files -S $path_to_dir/results/hisat2_align.sam --summary-file $path_to_dir/results/

## running with the input as list
#hisat2 -q -x $path_to_dir/ref_genome_grcm39/indexed_ref/genome -U $fastq_files_list -S $path_to_dir/results/hisat2_align.sam --summary-file $path_to_dir/results/

## running sample one only
#hisat2 -q -x $path_to_dir/ref_genome_grcm39/indexed_ref/genome -U 1.merged.fastq -S $path_to_dir/results/hisat2_align.sam --summary-file $path_to_dir/results/
## running with two saamples only
#hisat2 -q -x $path_to_dir/ref_genome_grcm39/indexed_ref/genome -U 1.merged.fastq 2.merged.fastq -S $path_to_dir/results/hisat2_align.sam --summary-file $path_to_dir/results/

for i in *.fastq
    do 
        #echo $i
        hisat2 -q -x $path_to_dir/ref_genome_grcm39/indexed_ref/genome -U $i -S $path_to_dir/results/hisat2_align.sam--summary-file $path_to_dir/results/

done

echo " "
echo " Finished aligning the reads"
duration=$SECONDS
echo $duration
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"

echo " "
echo "###############################################"



