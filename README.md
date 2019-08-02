# Introductory RNA-seq Using Command Line Programs and R

This repo contains my work done on a basic RNA sequencing project, following the example provided in Pertea M, Kim, Pertea G, Leek & Salzberg's [Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.nature.com/articles/nprot.2016.095). I knew virtually nothing about RNA-seq techniques going into this project, so a considerable amount of my time was spent reading and doing research to understand the concepts being presented. My motivation for wanting to do this was purely out of my own curiosity; computational biology and biological data analysis are fields that I'm deeply interested in, and I wanted to leverage my programming knowledge to give myself a better understanding of both.

In the future I would love to do a similar, slightly more ambitious project like this to further expand my understanding. There were so many tools and terms that came up throughout the process that were completely unfamiliar to me, but fortunately I felt confident enough in my programming abilities to carry me through when I felt stumped on the science. Overall, I'm happy with the results and simultaneously motivated to learn more.

## Some Context for Understanding Gene Expression Analysis

Starting out on this project, I took several hours to essentially surf Wikipedia and read a few academic articles to get familiar with some of the concepts being covered - mainly to justify to myself that I had some grasp of what was going on. From that, I've done my best to condense that down into what I think the most relevant information for this project is: namely, transcription and gene expression.

Transcription is a process where RNA polymerase generates complementary nucleotides to some template strand of DNA (with Uracil being substituted for Thymine). The RNA molecule only produces these nucleotides on the 3' end - or the "tail end" - adding one nucleotide at a time complementary to the template 3' -> 5' DNA strand. Thus, the 5' -> 3' RNA strand is identical to the coding 5' -> 3' DNA strand; for example, a coding DNA strand of "TAG" is transcribed as "ATC" on its complementary template strand, to which RNA would bind and write "UAG". Eventually, proteins bind to an AAUAAA sequence on the RNA, which signals where the pre-mRNA should be cut approximately 20 nucleotides down. After this, two caps are put on the RNA strand - a single G at the 5' end and hundreds of As at the 3' end - and intervening sequences (noncoding introns) are removed from the strand. Thus, what remains are the coding sequences, or the exons, which are joined together after the introns are purged. 

In the case of mRNA, these transcriptions carry coding for the synthesis of proteins. In the translation process, amino acids are chained together by triplets of nucleotides called codons to generate proteins, which then fold into their predefined three-dimensional structures. These two processes in concert make up the basis of the science of this project, because they are the building blocks for gene expression. While it is relatively common knowledge that genes are molecules with functions coded by nucleotides in DNA or RNA, gene expression is essentially the implicit amalgamation of that fact. The proteins synthesized are responsible for performing various functions throughout the cell or tissue, and depending on their prevalence they can generate varying manifestations.

Gene expression analysis is defined on [NCBI](https://www.ncbi.nlm.nih.gov/probe/docs/applexpression/) as "the determination of the pattern of genes expressed at the level of genetic transcription, under specific circumstances or in a specific cell." There are various approaches to measuring gene expression - the most obvious being observing the protein produced - but another approach is to measure mRNA and use these measurements to deduce expression levels. Gene annotation is the process of interpreting what certain coding regions of a given gene do, which is crucial for the deduction of expression levels. Finally, expression profiling is both of these processes on a large scale - instead of inferring singular proteins, we are interested in the measurement of thousands of genes in concert. Doing so gives us a broader understanding of the function or behavior of a cell.

This project involves RNA-seq, which is widely used to process huge quantities of data; it was designed with the intention to analyze all RNA molecules of a given cell and facilitate a profile of expression. However, the goal of this project is to analyze an abbreviated subset of the human X chromosome; the reason for this abbreviation is largely  due to processing constraints, as files for this sort of analysis are - as previously mentioned - often incredibly large. Because the X chromosome is relatively gene-rich, however, the subset should suffice for our purposes.

## Modules

I've separated my work into two separate markdown files here: the first describes the processes and commands performed in the terminal, and the second expands on the work done in R to explore the data and create the visualizations. In-depth explanations and a step-by-step approach are included in both.

* [Command Line](https://github.com/akweiss/RNA-seq-intro/blob/master/command-line.md)
* [R](https://github.com/akweiss/RNA-seq-intro/blob/master/R.md)

## Resources
Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL. [Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.nature.com/articles/nprot.2016.095).

Kim D, Langmead B and Salzberg SL. [HISAT: a fast spliced aligner with low memory requirements](https://www.nature.com/articles/nmeth.3317). *[Nature Methods 2015](https://www.nature.com/nmeth/).*

[StringTie: Transcript assembly and quantification for RNA-seq](https://ccb.jhu.edu/software/stringtie/)

[Ballgown: Flexible, isoform-level differential expression analysis](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html)
