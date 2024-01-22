# Shiny app for differential pathway methylation analysis

This app will conduct differential pathway methylation analysis using the LAM
approach described in the "manuscript" area of this repository.

## Input data requirements

For this to work, you will need a limma output file of a differential
methylation array analysis.
There are some criteria:

* The methylation array used is HM450K or EPIC v1.
Others are not supported at the moment.

* The file needs to be in TSV format.

* Shiny doesn't recognise normal R row names, so the probe identifier will need
to have its own column.
It needs to be the first column and needs its own row name.

* No spaces in the column names.

* The limma object needs one column with the name "t" which is what LAM uses
to score genes.

## Computational requirements

Computer with Docker installed.

CPU with 4 or more CPU threads and 32 GB RAM recommended.

A computer where you can launch a webservice, and access it via a custom port.

## Instructions

First pull the image.

```
docker pull mziemann/gmea
```

Run a container

```
docker run -e PASSWORD=bioc -p 8787:8787 mziemann/gmea
```

Now use a web browser to access the Rstudio instance by visiting localhost:8787.

Username is rstudio and passwd is bioc.

Next make your way to the terminal and type these commands.
It will copy the codes to an area where user rstudio can make modifications,
and receive recent updates.

```
cp /gmea ~ 
cd gmea
git pull
```

The shiny app is found in the "app" folder, so you can use the files pane (bottom
right) to open it up and then click "Run app".

This will pop up a new window where you can use the app.

Give your analysis a name and select the array type.
Just EPIC and HM450k are currently supported.

Select the gene sets you want to query.
Reactome is a good choice.
Take note that the GO gene set analysis will take a few minutes longer due to the
size of that dataset.

Select prioritisation scheme.
Significance is good when there are relatively few statistically significant sets.
Enrichment score is good when there are a large number of statistically significant
sets.
Enrichment score is a surrogate measure of effect size.

Upload your file and hit the "Analyze" button.

The analysis might take a few minutes.
You will know it is completed when the top pathway results appear, upon which you
can download the full results table and the HTML report.
The results table download will be nearly instant, but the report might take
a minute or so to generate.
Be patient and just hit the button once.

If you are unsure whether the app is still running an analysis, you can check the
CPU usage.
It is using >50% of one thread then it is probably still running something.
If it gets stuck after a few analyses, then just hit the "Reload app" button to
get a fresh app window.

## Gene information used here

The array annotations are quite old, and include lots of gene symbols that have
bee made defunct.
To solve this issue, I have used HGNChelper to update the gene symbols (our gene
info was downloaded on 21-Dec-2023).
Reactome, KEGG and GO gene sets were downloaded from MSigDB on 22-Jan-2024.
Here are the names:

* c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt

* c2.cp.reactome.v2023.2.Hs.symbols.gmt

* c5.all.v2023.2.Hs.symbols.gmt

## Session information

R version 4.3.2 (2023-10-31)

mitch_1.15.1

## Getting help

Post an issue to this repository if you have a question or having problems with the
app.
