# Bracken Lab RNA-Seq tutorial example

This notebook can be followed to get to grips with the basics of running a simple RNA-Seq analysis and to create bigwig tracks for use in UCSC or IGV.
The data used in this example comes from Evans 2019 PRC2 paper and uses the mutant day 2 and day 8 samples as an example. The end result of this code should result in similar to findings to those shown in figure S4 - with increased expression of differentiation genes in the D8 samples such as Gata4, Fgf5 and T.

## Step 1 - Set up
``` bash
#The example data and genome indices are already stored on the server in /data/evanRnaExample/ and /index/mm10/ respectively.
#We therefore only need to first create a directory in which we will work to carry out this analysis

#The following line of code will create a directory at ~/exampleRnaSeq/analysis (in your home directory) and includes the -p flag to ensure parent directories (exampleRnaSeq) are created also
mkdir -p ~/exampleRnaSeq/analysis
```


## Step 2 - Quality control checks
``` bash
#First we create a directory in which we will store the fastqc results
mkdir ~/exampleRnaSeq/analysis/fastqc

#Then we use fastQC to generate html quality summary reports for each file
fastqc -t 4 -o ~/exampleRnaSeq/analysis/fastqc /data/evanRnaExample/*

#The -t 4 flag here specifies to run on 4 computer threads to run on all 4 files at once.
#The -o flag tells the program to output results to ~/exampleRnaSeq/analysis/fastqc
#The final part of the command '/data/evanRnaExample/*' specifies the input files where a * denotes all files in that directory
```

## Step 3 - Inspecting FastQC reports
``` bash
#Open the html files in firefox. The server should forward the graphics of an internet browser window to your local computer
firefox ~/exampleRnaSeq/analysis/fastqc/KO_EB_D*.html

#These reports summarise the basic quality stats of a fastq file. The main plots to look at here are the first onewhich reports the per nucleotide sequence quality (in which you will see all nucleotide positions are of high quality in the green region of the plot) and the adapter sequence plot (in which you will see a tiny amount of adapter in the 3' end)

#Overall the sequence quality of this data is excellent and trimming the data isnt necessarily required however is included as the next step for demonstration purposes.
```

## Step 4 - Quality and adapter trimming
``` bash
#The program cutadapt can be used to quickly remove adapters and quality trim data in a single command. In this example we will run the program on each of the 4 files individually however I'd typically use a loop to do this on all the files in one line of code. 

mkdir ~/exampleRnaSeq/analysis/cutadapt

cutadapt -a AGATCGGAAGAG -q 10 -j 40 -m 20 -o ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D2_rep1.trim.fq.gz /data/evanRnaExample/KO_EB_D2_rep1.fq.gz 
cutadapt -a AGATCGGAAGAG -q 10 -j 40 -m 20 -o ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D2_rep2.trim.fq.gz /data/evanRnaExample/KO_EB_D2_rep2.fq.gz 
cutadapt -a AGATCGGAAGAG -q 10 -j 40 -m 20 -o ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D8_rep1.trim.fq.gz /data/evanRnaExample/KO_EB_D8_rep2.fq.gz 
cutadapt -a AGATCGGAAGAG -q 10 -j 40 -m 20 -o ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D8_rep2.trim.fq.gz /data/evanRnaExample/KO_EB_D8_rep2.fq.gz 

#The -a flag specifies which adapter to use. The example adapter here is the start of the universal Illumina adapter. If working with paired-end data you can specify the 5' adapter with -A.
#The -q flag says trim bases from the ends of reads with a quality below 10
#-j 40 tells the program to multi-thread and run on 40 CPU cores
#-m 20 specifys a minimum length post-trimming of 21 nucleotides. Reads below this are discarded as cannot be accurately aligned
#-o specifys the output file and the last parameter is the input.

#If working with paired data the read-pairs should be trimmed together by specifying a second output with -p and including both files as the last input. 

#After this stage you can rerun FastQC to see an improvement to the data quality using the trimmed files in ~/exampleRnaSeq/analysis/cutadapt/ as input.
```

## Step 5 - Alignment to the genome
``` bash
#A genome index has already been created for the mm10 genome so the alignment can be run without doing so. Equivalent indices for the human genome are also generated in the /index directory on the server
mkdir ~/exampleRnaSeq/analysis/alignments

#First unzip the compressed fastq files
unpigz -p 40 ~/exampleRnaSeq/analysis/cutadapt/*

#Align each fastq file to the mm10 genome
STAR --runMode alignReads --runThreadN 40 --readFilesIn ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D2_rep1.trim.fq --genomeDir /index/mm10 --outFileNamePrefix ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep1_ --outSAMtype BAM Unsorted 
STAR --runMode alignReads --runThreadN 40 --readFilesIn ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D2_rep2.trim.fq --genomeDir /index/mm10 --outFileNamePrefix ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep2_ --outSAMtype BAM Unsorted 
STAR --runMode alignReads --runThreadN 40 --readFilesIn ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D8_rep1.trim.fq --genomeDir /index/mm10 --outFileNamePrefix ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep1_ --outSAMtype BAM Unsorted 
STAR --runMode alignReads --runThreadN 40 --readFilesIn ~/exampleRnaSeq/analysis/cutadapt/KO_EB_D8_rep2.trim.fq --genomeDir /index/mm10 --outFileNamePrefix ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep2_ --outSAMtype BAM Unsorted

#Here we use the STAR aligner to align each sample to the genome. Most of the paramenters should be self explanatory based on the flag name. --outSAMtype BAM Unsorted skips a stage of analysis and saves disk space by allowing us to bypass saving alignments in the uncompressed SAM format. 

#We now no longer need the trimmed fastq files so delete them
rm ~/exampleRnaSeq/analysis/cutadapt/*
```

## Step 6 - Sorting BAM files
``` bash
samtools sort -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep1_Aligned.out.bam > ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep1_Aligned.out.sort.bam
samtools sort -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep2_Aligned.out.bam > ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep2_Aligned.out.sort.bam
samtools sort -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep1_Aligned.out.bam > ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep1_Aligned.out.sort.bam
samtools sort -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep2_Aligned.out.bam > ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep2_Aligned.out.sort.bam

#Here we use samtools to sort the bam files. STAR can output sorted BAM files using the command --outSAMtype BAM Sorted however it crashes when using lots of CPU threads to align the data and it is therefore quicker to do the analysis in this way. Here we use -@ 40 to use 40 cpu threads.

#index the BAM files for bigwig generation
samtools index -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep1_Aligned.out.sort.bam
samtools index -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep2_Aligned.out.sort.bam
samtools index -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep1_Aligned.out.sort.bam
samtools index -@ 40 ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep2_Aligned.out.sort.bam

#Remove the non sorted files
rm ~/exampleRnaSeq/analysis/alignments/*Aligned.out.bam§
```

## Step 7 - Quantifying Gene expression
``` bash
mkdir ~/exampleRnaSeq/analysis/featureCounts
featureCounts -T 40  -a /index/mm10/mm10.gtf -g gene_name -o ~/exampleRnaSeq/analysis/featureCounts/counts.txt ~/exampleRnaSeq/analysis/alignments/KO_*_Aligned.out.sort.bam 

#This single command  will quantify the overlaps of alignments in the BAM file with the gene features annotated in the GTF genome format. -T specifies computer processors, -a for the annotation in gtf format, -o the output file, and the last parameter is the input sorted bam files for each sample, using * to specify all files which match the pattern
* In this example we use -g gene_name to summarise counts to the gene names annotated in the GTF file rather than ensembl gene_ids. I would typically use gene ids and convert them to gene names in downstream analyses but for ease of interpretation we will use the gene_names in this example.

#In order to get this into the right format for differential expression analysis in R, we can use the following one-liner which will remove genome annotation columns, header lines and replace the name of the first coloum with nothing (required for matrix input in R)
cut -f 1,7,8,9,10 ~/exampleRnaSeq/analysis/featureCounts/counts.txt |sed 1d|sed 's/Geneid//g' > ~/exampleRnaSeq/analysis/featureCounts/counts.tsv
```

## Step 8 - Differential Expression analysis in R
``` bash
#Change our working directory to where the output files are saved
cd ~/exampleRnaSeq/analysis/featureCounts
#Start the R program
R
```
``` R
#load the DESeq2 library
library(DESeq2)

#Read the expression count data into a table called counts
counts<-read.table("counts.tsv")

#Rename the col headers for ease of interpretability (you could do this outside R also)
colnames(counts)<-c("KO_EB_D2_rep1","KO_EB_D2_rep2","KO_EB_D8_rep1","KO_EB_D8_rep2")
#In R we can specify a list of character strings using c("","") and assign objects such as lists into variables using <- or = 

#You can check the count table using
head(counts)
#Your output should look something like
#              KO_EB_D2_rep1 KO_EB_D2_rep2 KO_EB_D8_rep1 KO_EB_D8_rep2
#4933401J01Rik             0             0             0             0
#Gm26206                   0             0             0             0
#Xkr4                     94            44            39            39
#Gm18956                   0             0             0             0
#Gm37180                   0             0             0             0
#Gm37363                   0             2             0             0

#Create an experimental design matrix for DESeq
design<- cbind(condition=c('D2','D2','D8','D8'),samples=c("KO_EB_D2_rep1","KO_EB_D2_rep2","KO_EB_D8_rep1","KO_EB_D8_rep2"))
#Take a look at the matrix by typing the name of the object
design
#It will look like this
#     condition samples        
#[1,] "D2"      "KO_EB_D2_rep1"
#[2,] "D2"      "KO_EB_D2_rep2"
#[3,] "D8"      "KO_EB_D8_rep1"
#[4,] "D8"      "KO_EB_D8_rep2"
#This experimental design matrix is used by DESeq2 to know which samples are in which replicate condition. More complex experimental designs can be created at this stage such as multi-group analyses, removal of confounding batch effect factors and time series experiments but this is the most simple and most popular of designs. It is important the names of the samples are exactly the same as the col headers in the count matrix.
#In this case we can leave the design as it is as our control conditions name (D2) is first alphanumerically to the case conditions name (D8). If this wasnt the case (such as the names WT and MUT) then you would need to relevel the factors such that WT is first (google this if it would ever effect you).

#Create the DESeq data object
dds<-DESeqDataSetFromMatrix(countData=counts,colData=design, design= ~ condition)
#here the object dds is formed which contains the count data and the experimental design specified by the matrix we just made

#Calculate differential gene expression between D2 and D8
dds <- DESeq(dds)
#In R you can overwrite objects even if the data is used as the input in the functions whose results will overwrite it

#Save the results of the expression test to a new object
res<-results(dds)

#Filter the results for only those which are significant (optional)
sig <- subset(res, baseMean>20 & abs(log2FoldChange) > log2(1.5) & padj < 0.05 )
#This line of code will subset the res matrix for genes meeting the minimum expression level (basemean 20), with an absolute + or - fold change over 1.5 and a significant p value under 0.05.

#look at the results
sig
#If followed correctly, you should have 6319 rows remaining as significantly differentially expressed. 
#From here you can save the results of this to a file to look at further in excel or you can continue in R.
write.table(sig,"DESeq2results.tsv")

#Look at the top differentially expressed genes 
sig[order(-sig$log2FoldChange),]
#Here we use [] on a matrix to subset rows based on the criteria inside, where we reverse (-) sort on the log2fold change such that the top 5 upresulted are shown first before the bottom 5 downregulated:
#Top 5 - Col3a1, Hapln1, Postn, Tnx, Tnnt2
#Bottom 5 - Cyct, Bhlhe23, Pramel1, D030068K23Rik, Sycp1

#Look at the expression of the genes we highlighted at the start as in figure S4
res[c("Gata4","Fgf5","T"),]
#Here we subset on rownames by supplying a list within []. When doing this you have to remember to include the comma between the rownames you want and the columns, even if you arent specifying columns.

```

## Step 9 - Creating BW files
``` bash
mkdir ~/exampleRnaSeq/analysis/bigwigs
bamCoverage -b ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep1_Aligned.out.sort.bam --outFileName ~/exampleRnaSeq/analysis/bigwigs/KO_EB_D2_rep1.bw --numberOfProcessors 40 --normalizeUsing CPM
bamCoverage -b ~/exampleRnaSeq/analysis/alignments/KO_EB_D2_rep2_Aligned.out.sort.bam --outFileName ~/exampleRnaSeq/analysis/bigwigs/KO_EB_D2_rep2.bw --numberOfProcessors 40 --normalizeUsing CPM
bamCoverage -b ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep1_Aligned.out.sort.bam --outFileName ~/exampleRnaSeq/analysis/bigwigs/KO_EB_D8_rep1.bw --numberOfProcessors 40 --normalizeUsing CPM
bamCoverage -b ~/exampleRnaSeq/analysis/alignments/KO_EB_D8_rep2_Aligned.out.sort.bam --outFileName ~/exampleRnaSeq/analysis/bigwigs/KO_EB_D8_rep2.bw --numberOfProcessors 40 --normalizeUsing CPM
# Here we use the --normalizeUsing CPM flag to library size normalise the RNA-Seq tracks.

#cleanup by removing the alignment files you likely no longer need 
rm ~/exampleRnaSeq/analysis/alignments -r
#-r  tells the program to work recursively so that the directory can be removed

```
