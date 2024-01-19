# GMEA

## Summary

Infinium arrays are still widely used for profiling DNA methylation differences due to ageing, lifestyle,
development and disease [1].
Methods have been devised to conduct over-representation based functional enrichment analysis, yet
there are no widely recognised approaches to applying functional class scoring (FCS) methods like GSEA
which are thought to have better sensitivity.
Existing FCS methods do not explicitly consider hyper- and hypo-methylated gene pathways separately, which
is odd given direction of methylation changes is likely to have a bearing upon the direction of gene
regulation, which is key in understanding disease processes.
Conducting FCS analysis of infinium arrays is complicated due to the presence of multiple probes for each
gene.
In this work we examined various methods for performing FCS analysis of infinium methylation data using
simulated methylation profiles.
The method identified with the highest accuracy involved arithmetic mean aggregation of probe t-statistics
to gene level, followed by pathway enrichment with the mitch package that conducts an ANOVA-on-ranks test.
In our paper, we describe the simulation work conducted, and apply this method to various public datasets
including cancer, ageing, in vitro fertilisation, and a large scale EWAS involving 18k volunteers.
This work introduces a new pathway enrichment technique that will provide new epigenetic insights from
Infinium methylation array data.

About the name: we originally called this project Gene Methylation Enrichment Analysis (GMEA), but
currently don't think it is suitable.

## Contents

Important contents:

* figs: this contains all the data analysis scripts.

* manuscript: this contains an early version of the draft manuscript.
A live compiled version is found online [here](https://ziemann-lab.net/public/gmea/manuscript.html).

* example workflow: this contains a reproducible example workflow so you can see how the method works.
A live compiled version is found [here](https://ziemann-lab.net/public/gmea/example_workflow.html).

* app/trial: this is an attempt to create a shiny app (currently not working).

* Dockerfile: this was used to create a docker image which was used to execute all of the data analysis for the project.
This makes our work perfectly reproducible, and portable across different systems.
You can download the image from [Dockerhub](https://hub.docker.com/repository/docker/mziemann/gmea/general).

## Reproduction

You will need a Linux computer with 32 GB RAM and 20 GB of free space and docker installed.
You will also need to be in the docker group so you can use it without needing sudo password.
The docker container has Rstudio and all the packages we'll need to run this example.

In a terminal run the following command to start a container.

```
docker run -e PASSWORD=bioc -p 8787:8787 mziemann/gmea
```

Now visit the Rstudio webserver with your browser.
If you ran this on our local PC, then the address will be localhost:8787.
If you ran the command on another server, you will need to find its ip address and visit that (eg:123.456.789.42:8787).
You will be greeted with the Rstudio login the username is rstudio and password is bioc.

Next, copy the /gmea folder to the rstudio user home directory using the R console.

```
file.copy("/gmea",".",recursive=TRUE)
```

You will see the folder appear in the files pane.

In the terminal, execute this command to update the codes.

```
cd gmea && git pull
```

From there, enter the 'figs' directory where the scripts are.

```
cd figs
```

In this folder, the script called 'main.Rmd' will conduct all the analyses in the correct order,
and generate all the figures.
Be warned, that the analysis required to generate figures 1 and 2 takes approximately 4 weeks on a
server with 32 threads and 128 GB RAM.
(Alternatively, SLURM scripts are provided to speed this step up with a HPC cluster.)

```
Rscript -e 'rmarkdown::render("main.Rmd")'
```

## Contributions

If you have questions, or are experiencing problems applying this method or reproducing our work, please raise an issue on this repository.
