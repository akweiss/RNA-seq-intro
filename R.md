## Setting Up for Analysis in R Using Ballgown

Now that we have our Ballgown data from our command line programs, we can shift gears and move into R to do our analysis and generate some visualizations. Generally speaking, in differential gene expression analysis we are interested in performing statistical analyses in order to detect quantitative differences in expression levels between two experimental groups. In our case, we are interested in investigating the difference between the human X chromosome gene expression in males and females.

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

Now we can create our Ballgown object in R using the data compiled in the terminal (assigned to the dataDir parameter) and our phenotype data imported in the previous step (assigned to the pData parameter). From this object, we will also create bg_chrX_filt, which removes low abundance genes. As the name implies, low abundance genes are genes that match to our genome with lower frequency - therefore the bg_chrX_filt variable is identical to bg_chrX, with the exception of having the rarely-occuring genes weeded out. The "rowVars(texpr(bg_chrX)) > 1" translates to only including genes that appear in our data more than once. We will use the bg_chrX_filt object to perform the majority of our analysis.

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

## Data Exploration and Visualizations

First off, we can change our colors to help aid our visualizations:

``` sh
> colorscheme = c('hotpink', 'dodgerblue', 'darkorange', 'limegreen', 'yellow')
> palette(colorscheme)
```

Let's start with a plot of gene abundances across samples, measured as FPKM values. FPKM is "fragments per kilobase of exon model per million reads mapped" - it is an estimation of gene expression based on our data, calculated by counting the number of reads mapped to each gene sequence. We apply the transform of log2(fpkm + 1) because log(0) is undefined.

``` sh
> fpkm = texpr(bg_chrX, meas = "fpkm")
> fpkm = log2(fpkm + 1)
> boxplot(fpkm, col = as.numeric(pheno_data$sex), las = 2, ylab = 'log2(FPKM + 1)')
```

From this, we get the following visualization of our twelve samples with males being represented by blue and females represented by pink:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FPKM-samples-dist.png">
</p>

Here we can see that our medians (the bold bars within each box for each individual sample) are all nearly zero. Hmm. Let's try the same thing again, calculating FPKM with bg_chrX_filt this time.

``` sh
> fpkm = texpr(bg_chrX, meas = "fpkm")
> fpkm = log2(fpkm + 1)
> boxplot(fpkm, col = as.numeric(pheno_data$sex), las = 2, ylab = 'log2(FPKM + 1)')
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FPKM_filt_samples.png">
</p>

Now we see our medians all vary between 0 and ~3. Interesting! After some thought, this does make sense; reads that only appear once in our data can only be mapped once - if they map at all - and log2(1) = 0. Thus, this seems to imply that our bg_chrX data has a lot of genes with very low or zero abundance in our samples, while the bg_chrX_filt data removes these to give a better feel for the proportion of genes with higher abundance.

Next, let's investigate which transcripts missed our cutoff of q < 0.05. We can do this by looking at the first 20 rows of our results_transcript dataframe.

```sh
> head(results_transcripts, 20)
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/results_20.png">
</p>

From this, we can see that the gene ATP6AP2 narrowly missed our cutoff. Let's investigate this gene. From the readout, we know its transcript ID is 1042 - after a little investigation, we can find that in our bg_chrX_filt object, ATP6AP2 is indexed at position 797. Using this, we can generate some visualizations.

First off, let's use a boxplot to show FPKM distribution of the individual ATP6AP2 transcript NM_005765 across all samples.

``` sh
> plot(fpkm[797,] ~ pheno_data$sex, border = c(1,2), main = paste(ballgown::geneNames(bg_chrX_filt)[797],' : ', 
ballgown::transcriptNames(bg_chrX_filt)[797]), pch = 19, xlab = "Sex", ylab = "log2(FPKM + 1)")
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/ATP6AP2-diff.png">
</p>

While this might seem like a reasonable difference between males and females, note the scale of the y-axis: in actuality, it appears that the difference between medians of NM_005765 is approximately 0.5. In this way, the boxplot can appear a bit deceiving - we should consider manually rescaling the y-axis to better illustrate this FPKM distribution.

Going back to our transcripts_results readout above, we know the gene ID of ATP6AP2 is MSTRG.240. Thus, we can expand further by plotting the average expression levels for all transcripts of ATP6AP2 between males and females with the following:

``` sh
> plotMeans(MSTRG.240', bg_chrX_filt, groupvar = "sex", legend = FALSE)
```

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/ATP6AP2-means.png">
</p>

From this visualization, we now know that in our data the gene ATP6AP2 has five distinct isoforms. Isoforms are proteins that are functionally similar, but are not identical on their encoding. We can see the differential expression of ATP6AP2 in males and females with this side-by-side comparison, where "hotter" colors indicate higher expression. In this case, for instance, the fourth isoform seems to be the most highly expressed.

We can perform similar analyses for the genes PNPLA4 and FMR1.

Using identical commands as above (but with updated position within the bg_chrX_filt object and appropriate geneID), we get the following results for PNPLA4:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/PNPLA4-diff.png">
</p>

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/PNPLA4-means.png">
</p>

The case of PNPLA4 provides a good example of the difference between transcript level and gene level differential expression. Note that PNPLA4's gene ID, MSTRP.63, *is* listed in our statistically significant gene_results table, but is *not* in the transcript_results table (though it does appear in the top 20). Essentially, this means that the overall gene of PNPLA4 was found to be differentially expressed with statistical significance, but its transcripts were not (EDIT). From inspecting our boxplot of NM_001142389, this makes some sense; there is a tremendous amount of variation in females, to the extent that for this transcript nearly all of the data for the males is encompassed in the space just between the first quartile and the median for the females. This can make it difficult to conclude if there really is a significant difference. On the gene level, however, PNPLA4 appears to be more highly expressed in four out of the five distinct isoforms.

It's also worth noting that the gene ID for ATP6AP2, MSTRP.240, just *barely* missed the cutoff of q < 0.05 on the gene level. In my opinion, this - paired with also narrowly missing the cutoff on the transcript level - enables an argument to be made for our data showing ATP6AP2 to be differentially expressed.

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/gene_results_10.png">
</p>

Continuing on and following a similar approach, for FMR1 we get:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FMR1-diff.png">
</p>

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/FMR1-means.png">
</p>

FMR1 appears to be more statistically significant on transcript level. In fact, if we query our table of gene_results, it does not appear in even the top 30. To some degree, we can verify this through our plots - in the side-by-side comparison between males and females on the gene level, there is no clear perceptual distinction or "winner" between the two. But on the transcript level, we have two fairly compact boxplots where the median of the females approaches the maximum of the males. However, of the three genes we've examined the evidence for FMR1 is certainly the least compelling, and I don't think our data reasonably delineates FMR1 as a candidate for differential expression the same way that it does for PNPLA4 and ATP6AP2.

Lastly, let's look at a gene that is known to be expressed differentially and *did* show up on our list with a q-value < 0.05 cutoff: XIST. According to our reference article for this project, XIST is known to be more highly expressed in females than males. We can verify this quickly by producing a boxplot similar to those done above.

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/XIST_diff.png">
</p>

Yeah, that's... pretty convincing.

If we recall our readout for transcript_results, we can see that XIST has the ID of 2394:

<p align="center">
  <img width="700" src="https://github.com/akweiss/RNA-seq-intro/blob/master/images/RNA-seq-qval-2.png">
</p>

Thus, let's look at XIST in one of our female samples. From our first visualization (and our provided .csv), we know ERR188234 is female - therefore:

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

From these plots, we see that in our samples XIST has thirteen distinct isoforms. Of these, it appears that the eleventh isoform is the most highly expressed - which we can verify if we iterate through all of the female samples we have. I have included these plots in the [images](https://github.com/akweiss/RNA-seq-intro/tree/master/images) folder, for anyone interested in studying them in more detail. 

## Final Thoughts

From this project and our data, we can clearly conclude that both the NR_001564 transcript of the XIST gene and the XIST gene itself are differentially expressed between males and females. Additionally, PNPLA4 was found to be differentially expressed on the gene level, while ATP6AP2 just barely missed the cutoff in both categories. I think both of these genes have strong evidence of differential expression given our data, despite not making the cutoff of q < 0.05. It's worth noting that this cutoff is entirely arbitrary; if we had set it to q < 0.08, for instance, ATP6AP2 would be considered differentially expressed on both levels.

Overall, I'm happy with my analysis because my results were slightly different from the results provided in the tutorial. It was fun to explore my own data and make my own conclusions, knowing that if I ran the same command line protocol again I would get a different outcome and subsequently different data to analyze. Most of all, I feel like this analysis forced me to learn and understand the subject on a deeper level, which is largely why I found it so enjoyable.

## Resources
[EBI: FPKM](https://www.ebi.ac.uk/training/online/glossary/fpkm)
