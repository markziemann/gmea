## Run the example yourself.

You will need a Linux computer with 32 GB RAM and 20 GB of free space and docker installed.
You will also need to be in the docker group so you can use it without needing sudo password.
The docker container has Rstudio and all the packages we'll need to run this example.

In a terminal run the following command to start a container.

```
docker run -e PASSWORD=bioc -p 8787:8787 mziemann/gmea
```

Now visit the rstudio webserver with your browser.
If you ran this on our local PC, then the address will be `localhost:8787`.
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

In the files pane, click on the gmea folder to reveal the contents. 
Click on "example_workflow" and you will see there is the "example_workflow.Rmd" script
as well as a gmt file of gene ontologies and the DMPs.csv.gz file with the methylation
data.

Click on the "example_workflow.Rmd" script and the contents will appear in the text editor
window.
Click on the "Knit" button to execute it.
It might take a few minutes.
When it's complete you will also see that a html report called "example_mitchreport.html"
and "example_mitchplots.pdf" have been generated as part of this workflow.
It contains some data and charts to help interpret the enrichment results.

In order to customise it for your own data, replace the input file "DMPs.csv.gz" with your own data file.
Ensure you are using the appropriate gene table eg: EPIC array or 450K.
I have examples of each in the example_workflow.Rmd script.
