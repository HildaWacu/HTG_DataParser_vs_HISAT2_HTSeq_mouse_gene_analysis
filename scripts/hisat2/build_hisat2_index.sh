#!/bin/bash

# ## build an index of the reference
external_drive_path="/media/user/HildaWN/FASTQ_files/HISAT_alignment"
# external_drive_path="/Users/eht/Downloads/hilda_ref/HISAT_alignment"
# mkdir /media/user/HildaWN/FASTQ_files/HISAT_alignment/ndexed_ref

## extract exons and splice-sites

# hisat2-build $external_drive_path/ref_genome $external_drive_path/indexed_ref/genome  # where genome is the basename
#mkdir indexed_ref
### rename the .fna files, update the direcotries to the data .fna files

#echo $external_drive_path/ref_genome/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna 
#cp $external_drive_path/ref_genome/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna  $external_drive_path/ref_genome/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fa
#hisat2-build $path_to_dir/ref_genome_grcm39/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fa /Users/eht/Downloads/hilda_ref/HISAT_alignment/ref_genome_grcm39/indexed_ref/genome 

# if [ -z "$(ls -A /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome/indexed_ref)" ]; then
#    echo "Error: ref indexes not built"
#    exit 1
# else
#    echo "Buiding of ref indexes complete"
# fi
