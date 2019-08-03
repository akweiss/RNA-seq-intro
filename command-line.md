# Pipeline and Overview

The following framework is provided in the article [Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.nature.com/articles/nprot.2016.095):
<p align="center">
  <img width="300" height="500" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/RNA-seq%20Pipe.png">
</p>

This effectively illustrates the pipeline we will be implementing for our project. Generally speaking, the overall goal of this section is to assemble our data with the Ballgown package, as Ballgown acts as an intermediary between the command line software and R (where we will do our analysis and visualizations). Therefore, our objective here is to pass our initial data through the pipeline laid out above.

The data we are given initially are millions of reads of mRNA that have no non-coding introns; the HISAT2 software will take these reads and align them to the genome. Essentially, this means the software works to piece together where in the structure of the genome our collection of reads belongs (or "aligns"). We then will pass these alignments into StringTie, which will assemble our reads into transcripts. Because some of the assembled transcripts may only be partially described by our reads, we need pass our assemblies through StringTie's merge function in order to create a set of transcripts that is consistent across all samples. StringTie also provides estimates for expression levels for each gene and isoform both before and after the merging process, as well as information about relative abundance - both of these concepts I will discuss in more depth in the R section. For now, we can pass this information into gffcompare for useful statistics on our transcripts. StringTie can then compile coverage tables, which are descriptions of our reads or merged transcripts that align to or "cover" our known reference genome. Finally, these tables will be fed into Ballgown in R, which will take all of the information on transcripts and abundances and determine which genes and transcripts are differentially expressed between two experimental groups.

# Computations in Command Line

For starters, we need to download all of the appropriate software. This project uses HISAT2, StringTie, and SAMTools - we also need to install Ballgown through the R console. These were all fairly straightforward installations with clear instructions on their respective websites and READMEs. After unpacking and running the makefiles, we then need to copy the executables into our project directory. Then we can download the data (also supplied in the article), move it to our project directory, and get started on the first step.

The first objective in our pipeline is to read our RNA sequences into HISAT2 in order to generate our alignments. To reiterate, our goal here is map these millions of short reads of mRNA that we are given to our reference genome, where we will then attempt to pinpoint the position of where these short sequences belong within that genome. After navigating to our project directory, we can run the following commands for each of our twelve samples to do this:

``` sh
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188104_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188104_chrX_2.fastq.gz -S ERR188104_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188234_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188234_chrX_2.fastq.gz -S ERR188234_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188245_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188245_chrX_2.fastq.gz -S ERR188245_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188257_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188257_chrX_2.fastq.gz -S ERR188257_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188273_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188273_chrX_2.fastq.gz -S ERR188273_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188337_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188337_chrX_2.fastq.gz -S ERR188337_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188383_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188383_chrX_2.fastq.gz -S ERR188383_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188401_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188401_chrX_2.fastq.gz -S ERR188401_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188454_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR188454_chrX_2.fastq.gz -S ERR188454_chrX.sam
$ ./hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR204916_chrX_1.fastq.gz 
-2 chrX_data/samples/ERR204916_chrX_2.fastq.gz -S ERR204916_chrX.sam
```

After some trial and error (and a lot of unintentional typos at first), each command eventually yielded something like this in the terminal.

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/hisat2-terminal.png">
</p>

Nice! If interested, you can read about concordant and discordant pairing [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#concordant-pairs-match-pair-expectations-discordant-pairs-dont); essentially, concordant pairs behave in a way that is expected given our reference genome, and discordant pairs do not.

From our HISAT2 commands, we now have a collection of .SAM files. We need to convert these into .BAM to prep them for being run through StringTie. We can do that by running the following for each sample:

``` sh
$ samtools sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam
$ samtools sort -@ 8 -o ERR188104_chrX.bam ERR188104_chrX.sam
$ samtools sort -@ 8 -o ERR188234_chrX.bam ERR188234_chrX.sam
$ samtools sort -@ 8 -o ERR188245_chrX.bam ERR188245_chrX.sam
$ samtools sort -@ 8 -o ERR188257_chrX.bam ERR188257_chrX.sam
$ samtools sort -@ 8 -o ERR188273_chrX.bam ERR188273_chrX.sam
$ samtools sort -@ 8 -o ERR188337_chrX.bam ERR188337_chrX.sam
$ samtools sort -@ 8 -o ERR188383_chrX.bam ERR188383_chrX.sam
$ samtools sort -@ 8 -o ERR188401_chrX.bam ERR188401_chrX.sam
$ samtools sort -@ 8 -o ERR188428_chrX.bam ERR188428_chrX.sam
$ samtools sort -@ 8 -o ERR188454_chrX.bam ERR188454_chrX.sam
$ samtools sort -@ 8 -o ERR204916_chrX.bam ERR204916_chrX.sam
```

Now that we have our read alignments, let's look back at our pipeline and general outline. Our next step is transcript assembly. This is where we generate the full mRNA transcripts from our initial reads of only coding exons. To assemble our transcripts, we run the following StringTie commands:

``` sh
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188044_chrX.gtf –I ERR188044 ERR188044_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188104_chrX.gtf –I ERR188104 ERR188104_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188234_chrX.gtf –I ERR188234 ERR188234_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188245_chrX.gtf –I ERR188245 ERR188245_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188257_chrX.gtf –I ERR188257 ERR188257_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188273_chrX.gtf –I ERR188273 ERR188273_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188337_chrX.gtf –I ERR188337 ERR188337_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188383_chrX.gtf –I ERR188383 ERR188383_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188401_chrX.gtf –I ERR188401 ERR188401_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188428_chrX.gtf –I ERR188428 ERR188428_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188454_chrX.gtf –I ERR188454 ERR188454_chrX.bam
$ ./stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR204916_chrX.gtf –I ERR204916 ERR204916_chrX.bam
```

Now we can merge all of our samples with a single command. As a reminder, we run this merging command in an effort to account for transcriptions that were only partially described in the above step. Merging gives us a consistent set of transcripts for our eventual analysis.

``` sh
$ ./stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt
```

From here, our framework suggests that we could run gffcompare to generate transcription statistics. We can accomplish this by running:

``` sh
$ ./gffcompare -r chrX_data/genes/chrX.gtf -G -o merged stringtie_merged.gtf
```

Which yields:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/gffcompare-terminal.png">
</p>

We can follow this with the command:

``` sh
$ cat merged.stats
```

Which prints the following useful statistics for us.

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/gffcompare-terminal-2.png">
</p>

Sensitivity is the proportion of genes that are correctly constructed, which translates to a low false negative rate. Precision is the proportion of output that overlaps the annotation. High sensitivity and precision are crucial for our transcriptome assembler; if large quantities of our initial reads are not aligned, then the aligner will have significant difficulty reconstructing genes - particularly low abundance genes. In our case, for our merged transcriptions, the results indicate that we clearly have a low false negative rate, but the precision tells us that many of our transcripts are not within the list of StringTie's known transcripts. This means they are either false positives or *de novo* (new) transcripts.

Our last step is to create coverage tables, which we will pass into Ballgown in R (described in the following section). These coverage tables consolidate the description of our merged transcripts, comparing them to our reference genome to determine how much of that genome they cover.

``` sh
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188104/ERR188104_chrX.gtf ERR188104_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188234/ERR188234_chrX.gtf ERR188234_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188245/ERR188245_chrX.gtf ERR188245_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188257/ERR188257_chrX.gtf ERR188257_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188273/ERR188273_chrX.gtf ERR188273_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188337/ERR188337_chrX.gtf ERR188337_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188383/ERR188383_chrX.gtf ERR188383_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188401/ERR188401_chrX.gtf ERR188401_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188428/ERR188428_chrX.gtf ERR188428_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188454/ERR188454_chrX.gtf ERR188454_chrX.bam
$ ./stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR204916/ERR204916_chrX.gtf ERR204916_chrX.bam
```

With that, we're done with the command line portion of the project and can move on to our [analysis in R!](https://github.com/akweiss/RNA-seq-intro/blob/master/R.md)

# References
[EBI: Read mapping or alignment](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/read-mapping-or)

[RNA-seqlopedia](https://rnaseq.uoregon.edu/)

[Imperial College London: Read alignment](https://www.imperial.ac.uk/bioinformatics-data-science-group/resources/software/next-generation-sequencing-ngs-software/read-alignment/)

[Getting started with HISAT, StringTie, and Ballgown](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/)
