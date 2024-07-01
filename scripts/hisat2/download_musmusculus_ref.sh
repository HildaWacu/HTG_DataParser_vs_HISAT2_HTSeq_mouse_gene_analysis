#!/bin/bash

## To download the reference genome
### first create the ref  dir
#mkdir -p /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome

### move to the directory
#cd /media/user/HildaWN/FASTQ_files/HISAT_alignment/ref_genome

## ref genome
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000001635.27.zip" -H "Accept: application/zip"
unzip GCF_000001635.27.zip


# ## confirm is ref file has been successfully downloaded
# if [ -z "$(ls -A . )" ]; then
#    echo "Error: ref genome not downloaded"
#    exit 1
# else
#    echo "Ref genome downloaded successfully!"
# fi