## Run the example yourself.

You will need a Linux computer with 32 GB RAM and 20 GB of free space and docker installed.
You will also need to be in the docker group so you can use it without needing sudo password.
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

