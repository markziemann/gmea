# gmea

## Background



Infinium arrays are widely used for profiling DNA methylation differences due to ageing, lifestyle,
development and disease [1].
Methods have been devised to conduct over-representation based functional enrichment analysis, yet
there are no widely recognised approaches to applying functional class scoring (FCS) methods like GSEA
which are thought to have better sensitivity.
Conducting FCS analysis of infinium arrays is complicated due to the presence of multiple probes for each
gene.
Existing methods also do not explicitly consider hyper- and hypo-methylated gene pathways separately, which
is odd given direction of methylation changes is likely to have a bearing upon the direction of gene
regulation, which is key in understanding disease processes.

We consider three different approaches towards aggregating probe differential methylation to gene set
measurements.

1. limma-aggregate-enrichment (LA): In this approach, probe-based differential results from limma are
aggregated to gene level by taking the median t-statistic value for each gene.
This value can then be subject to downstream enrichment tests.

2. aggregate-limma-enrichment (AL): In this approach, the methylation values themselves for all probes are
aggregated (median) prior to differential analysis with limma.
The downstream differential methylation results are gene based rather than probe based.
These results can be subject to downstream enrichment tests.

3. aggregate-aggregate-limma (AA): In this approach, the methylation values are aggregated as in (2), but
these values then undergo another aggregation step to summarise gene set methylation.
These gene set methylation values then undergo differential methylation analysis using limma.

The purpose of this work is to determine which of these approaches is "best" and to compare the results
obtained with existing approaches/packages.

## Research plan

If we are to write a journal article on this work we will need the following:

1. Show that it is robust.
The three approaches outlines above will be tested with a large dataset that describes normal and cancer
samples from 37 patients.
From this dataset 18 patients will be subset randomly into two arms and analysed separately.
This will allow us to examine the proportion of common, uncommon and discordant pathways.
A high proportion of common pathways is good as it indicates high reproducibility.
Uncommon pathways are not good because they suggest lack of reproducibility.
Discordant pathways are very bad as they indicate obtaining contradictory results.
This experiment will indicate which of the approaches is best.
Make a venn diagram showing overlap between the results with the full 37 patients.

2. Show that it is sensitive.
Using the same cancer dataset, show that sensitivity is better with GMEA as compared to other approaches.
Simply downsample patient numbers and look at the recall.
Use GMEA, ebgsea and gsameth for comparison.

3. Examine the association between differential gene expression and methylation.
Does GMEA give better prediction of gene expression changes as compared to existing approaches?
[6.7]

TODO:

* Optimise the code for speed.
Currently it is still a bit slow and could benefit from more optimisation.
We could write some C helper code to do the heavy lifting, life FGSEA[8].

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
