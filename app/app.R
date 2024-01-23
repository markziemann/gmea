library("shiny")
library("vroom")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("mitch")

options(shiny.maxRequestSize = 300 * 1024^2) #300MB
# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)

arraytype <- c("EPIC","450k")

genesettype <- c("Reactome","KEGG MEDICUS","GO")

prioritisationtype <- c("Significance","Enrichment score")

ui <- fluidPage(
  titlePanel(
    "GMEA"
  ),
  sidebarPanel(
    "Welcome to the GMEA analysis server. Optionally give your analysis a name, select your array platform, your preferred gene set annotation set, upload your limma data set and hit the analyse button. Please be patient as analysis might take a few minutes. Keep in mind that only TSV is supported at this time. Please ensure that the first column is the probe ID, and that there is one column with the heading 't'. This app doesn't recognise default R row names, so please check that the sample t values shown match your data. After selecting the file to analyze, hit the 'Analyse!' button and once the top results are shown go ahead and download the full table and the HTML report.",
    textInput("dataname", "Dataset name:"),
    radioButtons("arraytype", "Array platform:", arraytype),
    radioButtons("genesettype", "Gene set database:", genesettype),
    radioButtons("prioritisationtype", "Prioritisation:", prioritisationtype),
    fileInput("upload", multiple = FALSE,label="Upload limma data:"),
    downloadButton("downloadData", label = "Download full results table"),
    downloadButton("downloadData2", label = "Download HTML report"),
    actionButton("analyse", "Analyse!",icon = icon("gears"))
  ),
  mainPanel(
    "This area will populate with results once analysis is complete.",
    textOutput("jobdata"),
    "file information",
    tableOutput("fileinfo"),
    "t values",
    tableOutput("tvals"),
    "Top differentially methylated pathways",
    dataTableOutput("mtable"),
    
  )
)

server <- function(input, output, session) {

  mydataname <- eventReactive(input$analyse, {
    as.character(input$dataname)
  })
    
  myarraytype <- eventReactive(input$analyse, {
    as.character(input$arraytype)
  })

  mygenesettype <- eventReactive(input$analyse, {
    as.character(input$genesettype)
  })
  
  myprioritisationtype <- eventReactive(input$analyse, {
    as.character(input$prioritisationtype)
  })  

  fileinfo <- eventReactive(input$analyse, {
    input$upload
  })
 

  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
 genetable <- reactive({
   if (myarraytype()=="EPIC") {
     gt <- readRDS("/home/app/epic.rds")
   } else {
     gt <- readRDS("/home/app/hm450k.rds")
   }
   gt
 })
    
  genesets <- reactive({
    if((mygenesettype())  == "Reactome") {
      gs <- gmt_import("/home/app/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
    }
    
    if((mygenesettype())  == "KEGG MEDICUS") {
      gs <- gmt_import("/home/app/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
    }
    
    if((mygenesettype())  == "GO") {
      gs <- gmt_import("/home/app/c5.all.v2023.2.Hs.symbols.gmt")
    }
    gs
  })
  
  prioritisation <- reactive({
    if((myprioritisationtype())  == "Significance") {
      prioritisation = "significance"
    } else {
      prioritisation = "effect"
    }
    prioritisation
  })
  
  mres <- reactive({
    m <- data()
    m <- cbind(m[,1],m[,"t"])
    rownames(m) <- m[,1]
    m[,1]=NULL
    m2 <- mitch_import(x=m,DEtype="prescored",geneTable=genetable())
    mres <- mitch_calc(x=m2,genesets=genesets(),minsetsize=5,cores=1, priority=prioritisation())    
  })
  
  mtable <- reactive({
        myres <- mres()
        myres$enrichment_result
  })
  
  output$jobdata <- renderText(paste(mydataname(),myarraytype(),mygenesettype()))
  
  output$fileinfo <- renderTable(fileinfo())
  
  output$tvals <- renderTable({
    req(input$upload)
    m <- data()
    m <- cbind(m[,1],m[,"t"])
    rownames(m) <- m[,1]
    m[,1]=NULL
    head(m)
  })
  
  output$mtable <- renderDataTable(
    mtable(),options = list(pageLength = 10, info = FALSE)
  )

  output$downloadData <- downloadHandler(
    filename = function() {
      paste(mydataname(), '.tsv', sep='')
    },
    content = function(file) {
      req(mtable())
      write.table(mtable(),file,sep="\t",row.names = FALSE)
    }
  )

  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste(mydataname(), '.html.zip', sep='')
    },
    content = function(file) {
      req(mtable())
      tmpdir <- tempdir()
      setwd(tempdir())
      mitch_report(res=mres(),outfile = "mitchreport.html",overwrite = TRUE)
      zip(zipfile=file, files="mitchreport.html")
    },
    contentType = "application/zip"
  )    

}

shinyApp(ui, server)
