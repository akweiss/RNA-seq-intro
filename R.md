## Analysis in R

Now that we have our Ballgown data from our command line programs, we can shift gears and move into R to do our analysis and generate some visualizations. Generally speaking, in differential gene expression analysis we are interested in performing statistical analyses in order to detect quantitative differences in expression levels between two experimental groups. In our case, we are interested in investigating the difference between the human X chromosome expression in males and females.

To get started, let's import the appropriate libraries into R.

``` sh
> library(ballgown)
> library(RSkittleBrewer)
> library(genefilter)
> library(dplyr)
> library(devtools)
```

Next we need to read our phenotype data into R, using the .csv also provided by the article outlining this project. This .csv contains our sample IDs in the rows, and our variables (sex and population) in the columns. Sex can take the value of male or female, and population can take the value of GBR (Great Britain) or YRI (Yoruba in Ibadan, Nigeria).

``` sh
> pheno_data = read.csv("~/bin/chrX_data/geuvadis_phenodata.csv")
```

Now we can create our Ballgown object in R using the data compiled in the terminal (assigned to the dataDir parameter) and our phenotype data imported in the previous step (assigned to the pData parameter). From this object, we will also create bg_chrX_filt, which removes low abundance genes. As the name implies, low abundance genes are genes that match to the genome with lower frequency - therefore the bg_chrX_filt variable is identical to bg_chrX, with the exception of having the rarely-occuring genes weeded out. The "rowVars(texpr(bg_chrX)) > 1" translates to only including genes that appear in our data more than once. We will use the bg_chrX_filt object to perform the majority of our analysis.

``` sh
> bg_chrX = ballgown(dataDir = "~/bin/ballgown", samplePattern = "ERR", pData = pheno_data)
> bg_chrX_filt = subset(bg_chrX, "rowVars(texpr(bg_chrX)) > 1", genomesubset = TRUE)
```

Using the stattest function built into Ballgown, we can test for statistically significant differences between our two experimental groups. That is, in results_transcripts we look for transcripts (feature = "transcript") that are differentially expressed between males and females, while simultaneously accounting for variation due to population (adjustvars = c("population")). Similarly, in results_genes we are looking to identify genes (feature = "gene") that show statistically significant differences between sexes while accounting for variation due to population.

``` sh
> results_transcripts = stattest(bg_chrX_filt, feature = "transcript", covariate = "sex", 
adjustvars = c("population"), getFC = TRUE, meas = "FPKM")
> results_genes = stattest(bg_chrX_filt, feature = "gene", covariate = "sex", 
adjustvars = c("population"), getFC = TRUE, meas = "FPKM")
```

Next, we will append the gene names and gene IDs to our results_transcripts dataframe.

``` sh
> results_transcripts = data.frame(geneNames = ballgown::geneNames(bg_chrX_filt), 
geneIDs = ballgown::geneIDs(bg_chrX_filt), results_transcripts)
```

We can then re-sort our results in order of their p-values, from smallest to largest.

``` sh
> results_transcripts = arrange(results_transcripts, pval)
> results_genes = arrange(results_genes, pval)
```

Now we can print our results with q-values < 0.05. Note that q-values are simply adjusted p-values that take into account false discovery rate; a q-value of 0.05 means we are willing to accept that 5% of our tests found to be statistically significant are actually false positives.

``` sh
> subset(results_transcripts, results_transcripts$qval < 0.05)
> subset(results_genes, results_genes$qval < 0.05)
```
These commands give us the following output:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/RNA-seq-qval-1.png">
</p>

From this, we now know that our X chromosome has eleven different transcripts that are differentially expressed between sexes, one of which corresponding to an isoform of a known gene (XIST). Additionally, we see that our X chromosome also has eight differentially expressed genes.

## Data Visualizations

First off, we can change our colors to help aid our visualizations:

``` sh
> colorscheme = c('hotpink', 'dodgerblue', 'darkorange', 'limegreen', 'yellow')
> palette(colorscheme)
```

Let's start with a plot of gene abundances across samples, measured as FPKM values. FPKM is "fragments per kilobase of exon model per million reads mapped" - it is an estimation of gene expression based on our data, calculated by counting the number of reads mapped to each gene sequence. 

``` sh
> fpkm = texpr(bg_chrX, meas = "fpkm")
> fpkm = log2(fpkm + 1)
> boxplot(fpkm, col = as.numeric(pheno_data$sex), las = 2, ylab = 'log2(FPKM + 1)')
```

From this, we get the following visualization with males being represented by blue and females represented by pink:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FPKM-samples-dist.png">
</p>

Here we can see that our medians (the bold bars within each box for each individual sample) are all nearly zero. Hmm. Intuitively, this seems to indicate to me that we have a LOT of low abundance genes kicking around in there. Let's try the same thing again, calculating FPKM with bg_chrX_filt this time.

``` sh
> fpkm = texpr(bg_chrX, meas = "fpkm")
> fpkm = log2(fpkm + 1)
> boxplot(fpkm, col = as.numeric(pheno_data$sex), las = 2, ylab = 'log2(FPKM + 1)')
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FPKM_filt_samples.png">
</p>

Now we see our medians all vary between 0 and ~3. Interesting! Implicitly, this tells us something about our transcription assembly process: the sensitivity was very good?

Next, let's investigate which genes missed our cutoff of q < 0.05. We can do this by looking at the first 20 rows of our results_transcript dataframe.

```sh
> head(results_transcripts, 20)
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/results_20.png">
</p>

From this, we can see that the gene ATP6AP2 narrowly missed our cutoff. Let's investigate this gene. From the readout, we know its ID is 1042 - after a little investigation, we can find that in our bg_chrX_filt object, ATP6AP2 is indexed at position 797. Using this, we can generate some visualizations.

First off, boxplot to sdlkfjsdklfj illustrate our differential expression.

``` sh
> plot(fpkm[797,] ~ pheno_data$sex, border = c(1,2), main = paste(ballgown::geneNames(bg_chrX_filt)[797],' : ', 
ballgown::transcriptNames(bg_chrX_filt)[797]), pch = 19, xlab = "Sex", ylab = "log2(FPKM + 1)")
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/ATP6AP2-diff.png">
</p>

Expand

Noting that from our above transcripts readout, we know the geneID of ATP6AP2 is MSTRG.240. Thus, we can expand further by plotting the average expression levels for all transcripts of ATP6AP2 between males and females with the following:

``` sh
> plotMeans(MSTRG.240', bg_chrX_filt, groupvar = "sex", legend = FALSE)
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/ATP6AP2-means.png">
</p>

Analysis here

We can perform similar analyses for the genes PNPLA4 and FMR1.

Using identical commands as above (but with updated position within the bg_chrX_filt object and appropriate geneID), we get the following results for PNPLA4:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/PNPLA4-diff.png">
</p>

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/PNPLA4-means.png">
</p>

Talk about transcript level vs gene level; .63 is PNP's gene ID, which DOES appear on the readout for results_genes.

Similarly, for FMR1 we have:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FMR1-diff.png">
</p>

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FMR1-means.png">
</p>

FMR1 matches more on transcript level. Does not even appear on the top 30 for gene comparison, which we can somewhat verify from the means plot.

Lastly, let's look at a gene that is known to be expressed differentially and *did* show up on our list with the q-value < 0.05 cutoff: XIST. According to our reference article for this project, XIST is known to be more highly expressed in females than males. We can verify this quickly by producing a boxplot similar to those done above.

IMAGE HERE.

If we recall our readout for transcript_results, we can see that XIST has the ID of 2394:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/RNA-seq-qval-2.png">
</p>

Thus, let's look at XIST in one of our female samples. From our .csv, we know ERR188234 is female - therefore:

``` sh
> plotTranscripts(ballgown::geneIDs(bg_chrX)[2394], bg_chrX, 
main = c('Gene XIST in Sample ERR188234'), sample = c('ERR188234'))
```

Gives us this plot:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/XIST_ERR188234_Vis.png">
</p>

Additionally, we can generate a similar plot of the expression of the XIST gene in our other samples by simply changing the "sample = c('ERR188234')" parameter (and changing the title of our plot for the purposes of clarity):

``` sh
> plotTranscripts(ballgown::geneIDs(bg_chrX)[2394], bg_chrX, 
main = c('Gene XIST in Sample ERR188428'), sample = c('ERR188428'))
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/XIST_ERR188428_Vis.png">
</p>

Analysis

# Final Thoughts

wrap it up homie

# Resources
[EBI: FPKM](https://www.ebi.ac.uk/training/online/glossary/fpkm)
