library("shiny")
library("vroom")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("mitch")

options(shiny.maxRequestSize = 100 * 1024^2) #100MB

arraytype <- c("450k", "EPIC")
genesettype <- c("Reactome","KEGG","GO")

ui <- fluidPage(
  titlePanel(
    "GMEA"
  ),
  sidebarPanel(
    textInput("dataname", "Dataset name:"),
    radioButtons("arraytype", "Array platform:", arraytype),
    radioButtons("genesettype", "Gene set database:", genesettype),
    fileInput("upload", multiple = FALSE,label="Upload limma data:"),
    actionButton("analyse", "Analyse!",icon = icon("gears"))
  ),
  mainPanel(
    textOutput("jobdata"),
    tableOutput("fileinfo"),
    tableOutput("head"),
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
 
  output$jobdata <- renderText(paste(mydataname(),myarraytype(),mygenesettype()))
  output$fileinfo <- renderTable(fileinfo())

  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  output$head <- renderTable({
    head(data(), 10)
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
    if((mygenesettype)  == "Reactome") {
      gmt_import("c2.cp.reactome.v2023.1.Hs.symbols.gmt")
    }
  })
  
  mtable <- reactive({
    m2 <- mitch_import(x=data(),DEtype="limma",geneTable=gt2())
    mres <- mitch_calc(x=m2,genesets=genesets(),minsetsize=5,cores=16, priority="effect")
    mres$enrichment_result
  })
  output$mtable <- renderDataTable(mtable())
}

shinyApp(ui, server)