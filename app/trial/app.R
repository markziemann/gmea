library("shiny")
library("vroom")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("mitch")

options(shiny.maxRequestSize = 300 * 1024^2) #300MB

arraytype <- c("EPIC","450k")

genesettype <- c("Reactome","KEGG","GO")

ui <- fluidPage(
  titlePanel(
    "GMEA"
  ),
  sidebarPanel(
    "Welcome to the GMEA analysis server. Optionally give your analysis a name, select your array platform, your preferred gene set annotation set, upload your limma data set and hit the analyse button. Please be patient as analysis might take a few minutes. Keep in mind that only TSV is supported at this time. Please ensure that the first column is the probe ID, and that there is one column with the heading 't'. This app doesn't default R row names, so please check that the sample t values shown match your data.",
    textInput("dataname", "Dataset name:"),
    radioButtons("arraytype", "Array platform:", arraytype),
    radioButtons("genesettype", "Gene set database:", genesettype),
    fileInput("upload", multiple = FALSE,label="Upload limma data:"),
    downloadLink("downloadData", "Download full table of results"),
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
    dataTableOutput("mtable")
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
  
  gt2 <- reactive({
    req(input$upload)
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    myann <- data.frame(anno[,c("UCSC_RefGene_Name","UCSC_RefGene_Group","Islands_Name","Relation_to_Island")])
    gp <- myann[,"UCSC_RefGene_Name",drop=FALSE]
    gp2 <- strsplit(gp$UCSC_RefGene_Name,";")
    names(gp2) <- rownames(gp)
    gp2 <- lapply(gp2,unique)
    gt2 <- stack(gp2)
    colnames(gt2) <- c("gene","probe")
    gt2$probe <- as.character(gt2$probe)
    gt2
  })
  
  genesets <- reactive({
    if((mygenesettype())  == "Reactome") {
      gmt_import("c2.cp.reactome.v2023.1.Hs.symbols.gmt")
    }
  })
  
  mres <- reactive({
    m <- data()
    m <- cbind(m[,1],m[,"t"])
    rownames(m) <- m[,1]
    m[,1]=NULL
    m2 <- mitch_import(x=m,DEtype="prescored",geneTable=gt2())
    mres <- mitch_calc(x=m2,genesets=genesets(),minsetsize=5,cores=1, priority="effect")    
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
      #write.csv(mtcars, file)
    }
  )
}

shinyApp(ui, server)
