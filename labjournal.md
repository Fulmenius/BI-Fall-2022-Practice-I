# BI-Fall-2022-Practice-I: determining antibiotic resistance mechanism in E.coli. In collaboration with @pavlovanadia #

## Step 1. Downloading the data ##

The reference sequence of *E. coli* was downloaded from NCBI FTP with command

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz`

The annotation of the reference sequence was downloaded from NCBI FTP with command

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz`


Raw Illumina sequencing reads (paired end run) of *E. coli* were downloaded with commands:

`wget https://figshare.com/ndownloader/files/23769689` 
*forward reads*

and

`wget https://figshare.com/ndownloader/files/23769692`
*reverse reads*

Then both files with reads were renamed so that the name would contain the .fastq.gz format with commands:

`mv 23769689 amp_res_1.fastq.gz` *# forward reads*

`mv 23769692 amp_res_2.fastq.gz` *# reverse reads*

All the downloaded files were copied from data directory to working directory:

`cp /home/nadia/BI/projects/project_1/BI-Fall-2022-Practice-I/rawdata/* /home/nadia/BI/projects/project_1/BI-Fall-2022-Practice-I/working_directory`