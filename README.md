# gmea

## Background

Gene Methylation Enrichment Analysis is a novel yet simple approach to detect epigenetic dysregulation.
Infinium arrays are widely used for profiling DNA methylation differences due to ageing, lifestyle, development and disease [1].
However the data provided at each position is rather noisy and therefore if we consider the contribution of all probes belonging to a gene, the signal to noise ratio will be higher [2].

In this project we will use methods normally reserved for gene set enrichment analysis [3,4] to detect genes that have probes “enriched” in the hyper or hypomethylated direction.
The idea is that we can aggregate information from all the probes for each gene we will get a clearer picture of whether a gene is differentially methylated or not.
Specifically, a Wilcox test on the limma t-statistics for each gene.
A self-contained test can run a one sample test to assess the probability that the mean is zero.
A competitive test can compare the values for one gene compared to all other genes.

## Research plan

If we are to write a journal article on this work we will need the following:

* Ability to detect differential gene methylation where existing methods cannot.
Existing methods include looking at individual probes, DMR finding, etc. 
To do this we will apply existing and GMEA method to an older dataset (eg BPROOF[5]) and identify novel trends.


* To show that GMEA is robust, split the BPROOF data in half to show the results are consistent.

* Examine the association between differential gene expression and methylation.
Does GMEA give better prediction of gene expression changes as compared to existing approaches?
We'll need datasets with matching infinium array and expression data; a brief search found a couple candidates [6.7]

* Optimise the code for speed.
Currently it is slow and could benefit from some optimisation.
We could write some C helper code to do the heavy lifting, life FGSEA[8].
Alternatively, we could write a new wilcox test function which strips off a lot of the unnecessary stuff.

* Provide it as an R/bioconductor package.
The epigenetics community will probably find this useful, especially with identifying slight trends in otherwise null data.

## References

1. Robinson MD, Kahraman A, Law CW, Lindsay H, Nowicka M, Weber LM, Zhou X. Statistical methods for detecting differentially methylated loci and regions. Front Genet. 2014 Sep 16;5:324. doi: 10.3389/fgene.2014.00324. PMID: 25278959; PMCID: PMC4165320.

2. Abraham G, Kowalczyk A, Loi S, Haviv I, Zobel J. Prediction of breast cancer prognosis using gene set statistics provides signature stability and biological context. BMC Bioinformatics. 2010 May 25;11:277. doi: 10.1186/1471-2105-11-277. PMID: 20500821; PMCID: PMC2895626.

3. Irizarry RA, Wang C, Zhou Y, Speed TP. Gene set enrichment analysis made simple. Stat Methods Med Res. 2009 Dec;18(6):565-75. doi: 10.1177/0962280209351908. Erratum in: Stat Methods Med Res. 2011 Oct;20(5):571. PMID: 20048385; PMCID: PMC3134237.

4. Wallace C. Gene set enrichment analysis using Wilcoxon tests [Internet]. Microsoft.com. 2012 [cited 2022 May 10]. Available from: https://mran.microsoft.com/snapshot/2015-04-25/web/packages/wgsea/vignettes/wgsea.pdf

5. Jung AY, Smulders Y, Verhoef P, Kok FJ, Blom H, Kok RM, Kampman E, Durga J. No effect of folic acid supplementation on global DNA methylation in men and women with moderately elevated homocysteine. PLoS One. 2011;6(9):e24976. doi: 10.1371/journal.pone.0024976. Epub 2011 Sep 23. PMID: 21966393; PMCID: PMC3179474.

6. Gevaert O. MethylMix: an R package for identifying DNA methylation-driven genes. Bioinformatics. 2015 Jun 1;31(11):1839-41. doi: 10.1093/bioinformatics/btv020. Epub 2015 Jan 20. PMID: 25609794; PMCID: PMC4443673.

7. Jiao Y, Widschwendter M, Teschendorff AE. A systems-level integrative framework for genome-wide DNA methylation and gene expression data identifies differential gene expression modules under epigenetic control. Bioinformatics. 2014 Aug 15;30(16):2360-6. doi: 10.1093/bioinformatics/btu316. Epub 2014 May 2. PMID: 24794928.

8. Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. Fast gene set enrichment analysis. bioRxiv 060012; doi: https://doi.org/10.1101/060012
