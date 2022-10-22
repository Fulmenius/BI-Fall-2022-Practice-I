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

## Step 2. Inspecting raw sequences manually ##

Raw sequence files, reference genome file and reference genome annotation file were opened without unpacking them and first 20 lines of each file were inspected for verification the accuracy of the file format with following commands:

`zcat GCF_000005845.2_ASM584v2_genomic.fna.gz | head -n 20` *# reference genome if in fasta format#

`zcat GCF_000005845.2_ASM584v2_genomic.gff.gz | head -n 20` *# command to open the reference genome annotation file*

`zcat amp_res_1.fastq.gz | head -n 20` *# command to open the forward raw Illumina reads which is in fastq format*

`zcat amp_res_2.fastq.gz | head -n 20` *# command to open the reverse raw Illumina reads which is in fastq format*


**Task from the materials**
*Use ‘cat’ to open the entire fasta reference file. Do you notice anything different about it?*

We used `zcat` so that we would not need to gunzip the reference genome. We used the command:

`zcat GCF_000005845.2_ASM584v2_genomic.fna.gz`

The contig is too big to fit in the terminal window, but it seems like there is only one contig in the file. It seems quite reasonable as long as bacterial genome is represented by a single circular chromosome. To make sure that in .fasta file of *E. coli* reference genome there is a single contig, we decided to count the number of symbols ">" in the file. Every single contig in .fasta files starts with the symbol ">" so the number of ">" equals to the number of contigs. We used the command:

`zcat GCF_000005845.2_ASM584v2_genomic.fna.gz`

The result is 1, so reference genome is a single contig.

The naxt task is to count the number of reads in raw Illumina data. Now we can not avoid gunzipping the raw data files so we execute the commands:

`gunzip amp_res_1.fastq.gz` *# to unzip forward reads*

`gunzip amp_res_2.fastq.gz` *# to unzip reverse reads*

Now we can count the number of lines in each file by executing commands:

`wc -l amp_res_1.fastq`

`wc -l amp_res_2.fastq`

The number of lines in each file is 1823504. Each fastq read consists of four lines, so, the number of reads in each file is 1823504/4 = 455876. 

## Step 3. FASTQC: reads inspecting. ##

The fastqc utilite was used to evaluate the quality of the readings.
We executed the command:

`fastqc -o . amp_res_1.fastq amp_res_2.fastq` 

and it created two .html files. The files *amp_res_1_fastqc.html* and *amp_res_2_fastqc.html* could be found in the working directory if you would like to see them.

# WE NEED TO DESCRIBE THE FASTQC RESULTS HERE THOROUGHLY!!! #

## Step 4. Filtering the reads ##

I have trimmomatic installed by executing the following command:

`sudo apt install trimmomatic`


