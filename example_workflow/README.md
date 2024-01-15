## Run the example yourself.

You will need a Linux computer with 32 GB RAM and 20 GB of free space and docker installed.
You will also need to be in the docker group so you can use it without needing sudo password.
Here we provide two approaches.
First one is using the Rstudio IDE and the other using the command line.

## Using Rstudio

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

In the files pane, click on the gmea folder to reveal the contents. Click on "example_workflow" and the
"example_workflow.Rmd" file to bring up the Rmarkdown script.

Click on the "Knit" button to execute it.
It might take a few minutes.

In order to customise it for your own data, replace the input file "DMPs.csv.gz" with your own data file.
Ensure you are using the appropriate gene table eg: EPIC array or 450K. I have examples of each in the example_workflow.Rmd script.

## Using the command line

In a terminal run the following to start a container and get a bash command prompt.

```
docker run -it mziemann/gmea /bin/bash
```

You are located in the /gmea main folder.
The contents may not have been updated for a while, so pull to get the latest code.

```
git pull
```

Make your way to the "example_workflow" folder and list the contents.

```
cd example_workflow/
ls
```

In the folder you will see the following: 

* "example_workflow.Rmd" script

* "DMPs.csv.gz" differential methylation data from limma

* "c5.go.v2023.2.Hs.symbols.gmt" Gene ontology gene sets

As the packages are already installed, we should be able to run the codes.
Use the `R` command to open a prompt, then type the following command.

```
rmarkdown::render("example_workflow.Rmd")
```

Which will generate a html output for you.
To visualise the html on your PC, exit R (`q()`)and the container (`exit`), and type the following command:

```
docker cp `docker ps -alq`:/gmea/example_workflow/example_workflow.html .
```

Then open it up in your favourite browser.

```
firefox example_workflow.html
```

You will see that mitch has generated some charts and a report for you.
You can fetch them using the following:

```
docker cp `docker ps -alq`:/gmea/example_workflow/example_mitchreport.html .
docker cp `docker ps -alq`:/gmea/example_workflow/example_mitchcharts.pdf .
```

In order to run this with your own data, you can replace the DMPs.csv.gz file of limma differential methylation
with your own data.

