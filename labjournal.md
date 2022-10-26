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

The quality of raw forward readings (file *amp_res_1_fastqc.html*) is following:
**Normal**: Basic Statistics, Per sequence quality scores, Per base N content, Sequence Length Distribution, Sequence Duplication Levels, Overrepresented Sequences, Adapter Content.
**Slightly abnormal**: Per base sequence content, per sequence GC content
**Very unusual**: per base sequence quality, per tile sequence quality

The quality of raw reverse readings (file *amp_res_2_fastqc.html*) is following:
**Normal**: Basic Statistics, Per sequence quality scores, Per base N content, Sequence Length Distribution, Sequence Duplication Levels, Overrepresented Sequences, Adapter Content.
**Slightly abnormal**: per tile sequence quality, Per base sequence content, per sequence GC content
**Very unusual**: per base sequence quality

## Step 4. Filtering the reads ##

I have trimmomatic installed by executing the following command:

`sudo apt install trimmomatic`

To run trimmomatic we need to know the path to trimmomatic.jar file. To find it we used command

`dpkg -L trimmomatic`

and in its output found the path needed. In our case it is is */usr/share/java/trimmomatic.jar*.

Then we filtered the reads with trimmomatic.

The parameters of the trimmomatic command were following:
- input forward readings file *amp_res_1.fastq*
- input reverse readings file *amp_res_2.fastq*
- output forward paired reads file name *forward_paired.fq*
- output forward unpaired reads file name *forward_unpaired.fq*
- output reverse paired reads file name *reverse_paired.fq*
- output reverse unpaired reads file name *reverse_unpaired.fq*
- from the start of each reading bases are cut if quality below 20
- from the end of each read bases are cut if quality below 20
- trim reads with sliding window of size 10 and average quality within the window 20
- minimun read length required is 20


`java -jar /usr/share/java/trimmomatic.jar PE amp_res_1.fastq amp_res_2.fastq forward_20_paired.fq forward_20_unpaired.fq reverse_20_paired.fq reverse_20_unpaired.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20`

We filteres the raw data second time with quality threshold in all filters 30 using the command:

`java -jar /usr/share/java/trimmomatic.jar PE amp_res_1.fastq amp_res_2.fastq forward_30_paired.fq forward_30_unpaired.fq reverse_30_paired.fq reverse_30_unpaired.fq LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:20`

After filtering to compare the readings quality before and after filtering the data we visualized forward and reverse reads with fastqc with the command:

`fastqc -o . forward_20_paired.fq reverse_20_paired.fq forward_30_paired.fq reverse_30_paired.fq`

To compare the readings quality before and after filtering we suggest to compare three .html files for forward and reverse readings:

**Forward readings**:
*Without filtering*  Per base sequence quality and Per tile sequence quality were very unusual.
*After filtering with quality threshold 20* Per tile sequence quality was very unusual.
*After filtering with quality threshold 30* Per tile sequence quality was very unusual.

Anyway, the main aim of filtering which was to filter the reads with low quality, is achieved, because the plot Per base sequence quality became normal.

**Reverse readings**:
*Without filtering*  Per base sequence quality was very unusual.
*After filtering with quality threshold 20* none of the parameters were very unusual.
*After filtering with quality threshold 30* none of the parameters were very unusual.

**We decided to choose for ruther analysis the result of the filtering with lower quality threshold which were the files forward_20_paired.fq and reverse_20_paired.fq**



## Step 5. Aligning sequences to reference ##

Firslty, we neede to indexate the reference genome.
For indexation we used `bwa index`. The command was

` bwa index GCF_000005845.2_ASM584v2_genomic.fna.gz`

After indexation the mentioned below files were created:

*GCF_000005845.2_ASM584v2_genomic.fna.gz.amb*
*GCF_000005845.2_ASM584v2_genomic.fna.gz.ann*
*GCF_000005845.2_ASM584v2_genomic.fna.gz.bwt*
*GCF_000005845.2_ASM584v2_genomic.fna.gz.pac*
*GCF_000005845.2_ASM584v2_genomic.fna.gz.sa*


After reference genome indezation we performed the alignment of out reads to the reference genome using `bwa mem`. We used the command:

`bwa mem GCF_000005845.2_ASM584v2_genomic.fna.gz forward_20_paired.fq reverse_20_paired.fq > alignment.sam`

Now the result of alignment is in the file *alignment.sam*.

Then we compressed the .sam file to .bam file with the command:

`samtools view -S -b alignment.sam > alignment.bam`

and ran the command:

`samtools flagstat alignment.bam` to see the basic alignment statistics.

99.87% of reads were mapped.

**Then we nedded to sort and index alignment.bam file**

To sort the .bam file we performed the command:

`samtools sort alignment.bam -o alignment_sorted.bam`

To index the sorted .bam file we performed the command:

`samtools index alignment_sorted.bam`

Then we visualized the mapping in IGV. For reference genome, we had to gunzip it before uploading to IGV:

`gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz`

## Step 6. Variant calling ##

For variant calling we performed the command:

`samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna alignment_sorted.bam > my.mpileup`

VarScan installation:

`sudo apt install varscan`

Then we did variant calling using the following command with allene frequence 50%:

`java -jar VarScan.v2.3.4.jar mpileup2snp my.mpileup --min-vae-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf`

## Step 7. Variant annotation ##

Database creation
First, we created empty text file snpEff.config vith one string:

`k12.genome : ecoli_K12`

Then we created folder for the database:

`mkdir -p data/k12`

Then we put there the .gbk file and renamed it (file was taken from [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz)): 

`gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz`
`cp GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk`

Then we created database

`java -jar snpEff.jar build -genbank -v k12`

After this we performed the annotation and redirected the StdOut to another .vcf file:

`java -jar snpEff.jar ann k12 VarScan_results.vcf > VarScan_results_annotated.vcf`

After that in our directory a new file *VarScan_results_annotated.vcf* appeared. We can upload it to IGV for obtaining more variant information.

## Step 8. List of variants ##

We found 6 variants:

1. position 93043 C -> G (forward strand)

GCC -> GGC

Ala -> Gly

gene name ftsI

missense variant

2. position 482698 A -> T (reverse strand)

CAG -> CTG

Gln -> Leu

gene name acrB

missense variant

3. position 852762 A -> G 

this variant is in non-coding region

4. position 1905761 G -> A (forward strand)

GGT -> GAT

Gly -> Asp

gene name mntP

missense variant

5. position 3535147 T -> G (reverse strand)

GTA -> GGA

Val -> Gly

gene name envZ

missense variant

6. position 4390754 C -> A (reverse strand)

GCC -> GCA

Ala -> Ala 

synonymous variant

gene name rsgA