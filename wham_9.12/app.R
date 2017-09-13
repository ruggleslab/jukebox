##R Shiny gfam app

## update 8/28 correlation plots have been fixed, now we need to adjust input to reflect
## new python scripts and univ input!!
#!


library(shiny)
library(ggplot2)
library(plyr)
library(rowr)
library(plotly)
library(grid)
library(base)
library(shinythemes)

## ggplot legend extract
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


### corr plot func

get_plot_output_list <- function(max_plots, input_n, mat_list, sym_list, groupNames) {
  # Insert plot output objects the list
  plot_output_list <- lapply(1:input_n, function(i) {
    plotname <- paste("plot", i, sep="")
    savename <- paste("save", i, sep="")
    plot_output_object <- plotOutput(plotname, height = 280, width = 250)
    plot_output_object <- renderPlot({
      plot_mat = mat_list[[i]]
      plot_mat[is.na(plot_mat)] <- 0
      mat_title = c(paste(groupNames[i], "Correlation Values", sep=" "))
      sym_mat = sym_list[[i]]
      heatmap.2(plot_mat,
                cellnote = sym_mat, notecol="black",
                main=mat_title,
                density.info="none",
                #key = TRUE,
                #keysize = 1.0,
                breaks = seq(-1, 1, by = 0.02),
                col=c(colorRampPalette(c("blue", "white", "red"))(n=100)),
                dendrogram="both",
                margins=c(5,5),
                cexRow=1,
                cexCol=1.2,
                trace=c("none"),
                na.color="gray60"
      )
    })
  })

  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
}

download_plot_output_list <- function(max_plots, input_n, mat_list, sym_list, groupNames) {
  # Insert plot output objects the list
  plot_output_list <- lapply(1:input_n, function(i) {
    savename <- paste("save", i, sep="")
    plot_mat = mat_list[[i]]
    plot_mat[is.na(plot_mat)] <- 0
    mat_title = c(paste(groupNames[i], "Correlation Values", sep=" "))
    sym_mat = sym_list[[i]]
    heatmap.2(plot_mat,
              cellnote = sym_mat, notecol="black",
              main=mat_title,
              density.info="none",
              #key = TRUE,
              #keysize = 1.0,
              breaks = seq(-1, 1, by = 0.02),
              col=c(colorRampPalette(c("blue", "white", "red"))(n=100)),
              dendrogram="both",
              margins=c(5,5),
              cexRow=1,
              cexCol=1.2,
              trace=c("none"),
              na.color="gray60"
    )
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
}


max_plots <- 10
##### gif loader function

loadingLogo <- function(href, src, loadingsrc, height = NULL, width = NULL, alt = NULL) {
  tagList(
    tags$head(
      tags$script(
        "setInterval(function(){
        if ($('html').attr('class')=='shiny-busy') {
        $('div.busy').show();
        $('div.notbusy').hide();
        } else {
        $('div.busy').hide();
        $('div.notbusy').show();
        }
},100)")
  ),
  tags$a(href=href,
         div(class = "busy",
             img(src=loadingsrc,height = height, width = width, alt = alt)),
         div(class = 'notbusy',
             img(src = src, height = height, width = width, alt = alt))
         )
  )
}


## Initiate multiplot function ##

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#allows user to choose samples
sampleUploadUI <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        column(8,
             #makes checkbox of experiment choices
             #uiOutput(ns("group_pre")),
             uiOutput(ns("exptSamples"))
             )
        )
      )
    )
}

sampleUpload <- function(input,output,session, acc_nums){
  #user clicks to group samples
  output$exptSamples <- renderUI({
    req(acc_nums())
    ns <- session$ns
    selectInput(ns("exptallSamples"), 'Group', colnames(acc_nums()), multiple=TRUE, selectize=TRUE)
  })
  # output$group_pre <- renderUI({
  #   req(acc_nums())
  #   ns <- session$ns
  #   textInput(ns("group_pre_samps"), "Group Prefix")
  # })
  
  # output$sampName <- renderUI({
  #   req(acc_nums())
  #   ns <- session$ns
  #   textInput(ns("sampName"), label = 'Group Name')
  # })
  
  #makes a list of grouped samples that can be used to reorder
  #original matrix
  listcond <- reactive({
    cond2_t = input$exptallSamples
    transcription_conditions = list(cond2_t)
    return(transcription_conditions)
  })
  return(listcond)
}

## multiplot initiated!
#system('python acc_collapse.py')
#system('python species_collapse.py')

#acc_data <- read.delim("Acc_collapsed.tsv", header = TRUE)
#acc_data2 = subset(acc_data, Gene.Family != " NO_NAME")
#acc_data3 = subset(acc_data2, Species.Number > 10)
#spec_data <- read.delim("gfam_collapsed.tsv", header=TRUE)
#spec_data2 = subset(spec_data, Gene.Family != " NO_NAME")

ui <- navbarPage(title = "Workflow Hub for Automated Metagenomic Exploration",
                 theme = shinytheme("superhero"),
                       tabPanel("Home",
                                  fluidRow(column(8, 
                                                  fluidRow(column(12, h3("Methods Workflow"),
                                                                  mainPanel(img(src='methodsKVR.png', height = 842, width = 650))))),
                                           column(4,
                                                  fluidRow(column(12, h3("Resources"))),
                                                  fluidRow(column(12, h3("External Information")))))
                                ),
                       tabPanel("Upload",
                                fluidPage(
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput('file1', 'Choose TSV File',
                                                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                      tags$hr(),
                                      #checkboxInput('header', 'Header', TRUE),
                                      radioButtons('input_type', "Input Type",
                                                   choiceNames = c("Wham", "HUMANn2"),
                                                   choiceValues = c("Wham", "HUMANn2")),
                                      radioButtons('sep', 'Separator',
                                                   c(Tab='\t'),
                                                   '\t'),
                                      checkboxInput("testme", "Try a Sample Dataset!", 
                                                    value = FALSE)
                                    ),
                                    mainPanel(
                                      dataTableOutput('contents')
                                    )
                                  ))
                                
                       ),
                        tabPanel("Groups",
                                 fluidPage(
                                   titlePanel("Assign each sample to a Group"),
                                   sidebarLayout(
                                     sidebarPanel(
                                       numericInput("numInputs", "Select Number of Groups", 1, min = 1, max = 10),
                                       uiOutput("group_pre")
                                     ),
                                    mainPanel(
                                   fluidRow(
                                     h4("Group Selection"),
                                     textOutput("tutorialGroup"),
                                     tags$head(tags$style(
                                       "#tutorialGroup{color: red; font-size: 18px}")),
                                     textOutput("group_warning"),
                                     tags$head(tags$style(
                                       "#group_warning{color: red; font-size: 18px}")),
                                     # place to hold dynamic inputs
                                     uiOutput("inputGroup"))
                                   )
                                 )
                        )),
                        tabPanel("Explore Taxa",
                                 mainPanel(
                                   fluidRow(numericInput("taxaDims", "How many taxa levels are present in your data?", 1, min = 1, max = 2),
                                            textInput("taxaSep", "If more than one level, what character seperates your taxa levels? \n Permitted characters include { . , _ }")),
                                   fluidRow(uiOutput("TaxaDimExp")),
                                   fluidRow(plotlyOutput("species_explore")),
                                   fluidRow(downloadButton("species_download", "Download Plot"), 
                                            downloadButton("species_legend_download", "Download Legend")),
                                   width = 12)
                                 ),
                        tabPanel("Explore Genes",
                                 textOutput("instructions"),
                                 fluidRow(plotlyOutput("gene_explore")),
                                 fluidRow(downloadButton("gene_explore_download", "Download Plot"))
                                 ),
                        tabPanel("Expression",
                                 fluidRow(column(4,
                                                 selectizeInput('acc_list', choices=NULL,
                                                         label = h3("Begin by selecting pathways of interest"),
                                                         multiple = TRUE)),
                                          column(10,
                                                 checkboxInput('xy_switch', label = 'group by Gene Family?')
                                 )),
                                 hr(),
                                 fluidRow(plotOutput("plot1")),
                                 fluidRow(downloadButton("expression_download", "Download Plot")),
                                 fluidRow(column(6, align= "center", tableOutput("exp_table"))),
                                 fluidRow(downloadButton("expression_table_download", "Download Table"))
                                 ),
                        tabPanel("Taxa",
                                 fluidRow(column(12, uiOutput("spec_plot"))),
                                 fluidRow(downloadButton("expression_taxa_download", "Download Plots"))
                        ),
                        tabPanel("Correlation", selectizeInput('sig_select', choices=NULL,
                                                               label = h3("Begin by selecting two pathways of interest"),
                                                               multiple = TRUE),
                                 fluidPage(column(4, verbatimTextOutput("sig_message"))),
                                 fluidPage(column(10, plotOutput("corr_plot"))),
                                 fluidRow(downloadButton("corr_download", "Download Plot")),
                                 fluidPage(column(10, uiOutput("group_corrs"))),
                                 fluidPage(uiOutput("group_download")),
                                 fluidPage(column(10, tableOutput("corr_labels"))),
                                 fluidPage(downloadButton("corr_table_download", "Download Labels"))),
                        tabPanel(title = loadingLogo("http://google.com", "wham_logo_trans.png",
                                                            'wham_grey_inf.gif', height = 135, width = 280)),
                 tags$head(tags$style('.navbar {
                            font-size: 18px}', '.navbar-brand {font-size:32px}'))
                                                            
)

### control size of input file
options(shiny.maxRequestSize=500*1024^2)

`%then%` <- shiny:::`%OR%`

server <- function(input, output, session) {
  file_test <- reactive({
    if (input$testme){
      args = "univ_input_top.tsv"
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
      args = c(input$file1)
    }
    command ="python"
    acc_script="univ_input_gfam.py"
    
    # Add path to script as first arg
    accArgs = c(acc_script, args)
    
    result <- tryCatch(
      {system2(command, args=accArgs, stdout=TRUE)},
      warning=function(cond) {
        message(paste("Invalid File Type"))
        message(cond)
        return("PYTHON_ERROR")}
    )
    result
  })
  
  output$contents <- renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects and uploads a 
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
    # columns. The 'datapath' column will contain the local filenames where the 
    # data can be found.
    file_warning = file_test()
    if (file_warning == "PYTHON_ERROR"){
      message = c(rep("Invalid File Input", 50))
      message_mat = matrix(message, nrow = 10, ncol = 5)
      full_file = data.frame(message_mat)
    }
    else if (input$testme) {
      full_file <- read.csv("univ_input_top.tsv", header=TRUE, sep=input$sep)
      #full_file
    }
    else {
      inFile <- input$file1
      if (is.null(inFile)) {
        return(NULL)}
      full_file = read.csv(inFile$datapath, header=TRUE, sep=input$sep)
      #full_file
    }
    full_file
  })
  sample_ids <- reactive({
    if (input$testme) {
      full_file <- read.csv("univ_input_top.tsv", header=TRUE, sep=input$sep)
      FullFile = data.frame(full_file)
      samples = c(colnames(FullFile))
      samp_num = length(samples)
      samps = samples[4:samp_num]
      samps
      }
    else {
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)
      full_file = read.csv(inFile$datapath, header=TRUE, sep=input$sep)
      FullFile = data.frame(full_file)
      samples = c(colnames(FullFile))
      samp_num = length(samples)
      samps = samples[4:samp_num]
      samps
    }
  })
  
  acc_select <- reactive({
    if (input$testme){
      args = "univ_input_top.tsv"
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
      args = c(input$file1)
    }
    command ="python"
    acc_script="univ_input_gfam.py"
    
    # Add path to script as first arg
    accArgs = c(acc_script, args)
    
    result <- tryCatch(
      {system2(command, args=accArgs, stdout=TRUE)},
      warning=function(cond) {
        message(paste("Invalid File Type"))
        message(cond)
        return("PYTHON_ERROR")}
    )
    validate(
      need(result != "PYTHON_ERROR", "Invalid File Type")
    )
    acc_data = system2(command, args=accArgs, stdout=TRUE)
    acc_dat = acc_data[-1]
    acc_da = acc_dat[-1]
    accs  = read.table(text = acc_da, sep = '\t')
    #samples = paste("sample", 1:82, sep='')
    colnames(accs) = c("Acc", "Gene.Family", "SpeciesNumber", "Species", sample_ids())
    accs$Gene.Family = paste(accs$Gene.Family, " ", sep = "")
    accs$Gene.Family
  })
  observe({
    if (input$testme) {
      updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
      }
    else {
      if(is.null(input$file1)){}
      else {
        updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
      }
    }
  })
  acc_full <- reactive({
    if (input$testme){
      args = "univ_input_top.tsv"
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
      args = c(input$file1)
    }
    command ="python"
    acc_script="univ_input_gfam.py"


    # Add path to script as first arg
    accArgs = c(acc_script, args)
    
    result <- tryCatch(
      {system2(command, args=accArgs, stdout=TRUE)},
      warning=function(cond) {
        message(paste("Invalid File Type"))
        message(cond)
        return("PYTHON_ERROR")}
    )
    validate(
      need(result != "PYTHON_ERROR", "Invalid File Type")
    )
    
    acc_data = system2(command, args=accArgs, stdout=TRUE)
    acc_dat = acc_data[-1]
    acc_da = acc_dat[-1]
    #print(sample_ids())
    accs  = read.table(text = acc_da, sep = '\t')
    #samples = paste("sample", 1:82, sep='')
    colnames(accs) = c("Acc", "Gene.Family", "SpeciesNumber", "Species", sample_ids())
    accs
  })
  
  acc_nums <- reactive({
    if (input$testme){
      args = "univ_input_top.tsv"
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
      args = c(input$file1)
    }
    command ="python"
    acc_script="univ_input_gfam.py"
    
    # Add path to script as first arg
    accArgs = c(acc_script, args)
    
    acc_data = system2(command, args=accArgs, stdout=TRUE)
    acc_dat = acc_data[-1]
    acc_da = acc_dat[-1]
    #print(sample_ids())
    accs  = read.table(text = acc_da, sep = '\t')
    #samples = paste("sample", 1:82, sep='')
    colnames(accs) = c("Acc", "Gene.Family", "SpeciesNumber", "Species", sample_ids())
    accs2 = accs[,5:ncol(accs)]
    rownames(accs2) = accs$Gene.Family
    accs3 = as.matrix(accs2)
    accs3
  })
  
  ###input grouped samples based on number of groups, output reordered matrix###  
  results <- c()
  makeReactiveBinding('results')
  
  # create number of  columns based on input number of groups
  observeEvent(input$numInputs, {
    output$inputGroup = renderUI({
      if (input$testme){
        updateNumericInput(session, 'numInputs', value = 4)
        #updateTextInput(session, 'group_pre', value = "Antibiotic")
        group1 = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
                   "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
                   "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
                   "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
        group2 = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
                   "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
                   "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
                   "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
        group3 = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
                   "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
                   "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
                   "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
        group4 = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
                   "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
                   "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
                   "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")
        fluidRow(selectInput("test_groups", label = "",
                             choices = group1, multiple = TRUE, selected = group1),
                 selectInput("test_groups", label = "", choices = group2, 
                             multiple = TRUE, selected = group2),
                 selectInput("test_groups", label = "", choices = group3, 
                             multiple = TRUE, selected = group3),
                 selectInput("test_groups", label = "", choices = group4, 
                             multiple = TRUE, selected = group4))
      } else {
        validate(
          need(input$file1, 'Please provide a file in the Upload Tab')) 
        input_list <- lapply(1:input$numInputs, function(i) {
          inputName <- paste0("input", i)
          sampleUploadUI(inputName)
          })
      }
    })
    
    results <<- lapply(1:input$numInputs, function(i) {
      inputName <- paste0("input", i)
      callModule(sampleUpload, inputName, acc_nums)
    })
  })
  
  output$group_pre = renderUI({
    lapply(1:input$numInputs, function(i) {
      inputName <- paste0("group", i)
      textInput(inputName, label = "Group Prefix")
    })
  })
  
  group_names <- reactive({
    if (input$testme){
      groups <- c("Arm", "Vagina", "Saliva", "Stool")
    }
    else{
      groups <- sapply(1:input$numInputs, function(i){
        input[[paste0("group", i)]][1]})
    }
    return(groups)
  })

  #group results
  grouped_samps <- reactive({
    if (input$testme) {
      group1 = c("SRR532024_Abundance.RPKs", "SRR532015_Abundance.RPKs", "SRR532006_Abundance.RPKs",
                 "SRR532040_Abundance.RPKs", "SRR532504_Abundance.RPKs", "SRR532507_Abundance.RPKs",
                 "SRR638753_Abundance.RPKs", "SRR640357_Abundance.RPKs", "SRR640340_Abundance.RPKs",
                 "SRR640452_Abundance.RPKs", "SRR640499_Abundance.RPKs", "SRR545546_Abundance.RPKs")
      group2 = c("SRR062353_Abundance.RPKs", "SRR062357_Abundance.RPKs", "SRR062301_Abundance.RPKs",
                 "SRR062276_Abundance.RPKs", "SRR1804686_Abundance.RPKs", "SRR1804628_Abundance.RPKs",
                 "SRR514191_Abundance.RPKs", "SRR514180_Abundance.RPKs", "SRR513168_Abundance.RPKs",
                 "SRR514231_Abundance.RPKs", "SRR513448_Abundance.RPKs")
      group3 = c("SRR062435_Abundance.RPKs", "SRR062441_Abundance.RPKs", "SRR062389_Abundance.RPKs",
                 "SRR062413_Abundance.RPKs", "SRR062402_Abundance.RPKs", "SRR062396_Abundance.RPKs",
                 "SRR346673_Abundance.RPKs", "SRR346681_Abundance.RPKs", "SRR062371_Abundance.RPKs",
                 "SRR062372_Abundance.RPKs", "SRR062462_Abundance.RPKs", "SRR062415_Abundance.RPKs")
      group4 = c("SRR528423_Abundance.RPKs", "SRR528353_Abundance.RPKs", "SRR528300_Abundance.RPKs",
                 "SRR528261_Abundance.RPKs", "SRR528183_Abundance.RPKs", "SRR528155_Abundance.RPKs",
                 "SRR532178_Abundance.RPKs", "SRR532183_Abundance.RPKs", "SRR532190_Abundance.RPKs",
                 "SRR532191_Abundance.RPKs", "SRR533152_Abundance.RPKs", "SRR533153_Abundance.RPKs")
      g1 = list(group1)
      g2 = list(group2)
      g3 = list(group3)
      g4 = list(group4)
      glist = list(g1,g2,g3,g4)
      glist
    } else {
      lapply(1:input$numInputs, function(i) {
        results[[i]]()})
    }
  }) 

  output$group_warning <- renderText({
    validate(
      need(grouped_samps(), 'Select groups!')
    )
    group_list = c(unlist(grouped_samps()))
    library(base)
    duplicates = c(duplicated(group_list))
    if ('TRUE' %in% duplicates) {
      message = c("Warning: A sample has been assigned to more than one group!")
    } else {
      message = c(" ")
    }
    message
  })
  
  output$tutorialGroup <- renderText({
    if (input$testme){
      message = c("Groups have been named and assigned!\n Please continue to the next Tab")
      message
    }
  })
  
  # group_names <- reactive({
  #   if (input$testme) {
  #     polynames <- sapply(1:input$numInputs, function(i){
  #       c(paste('Antibiotic', 'group', i, sep='_'))})
  #   }else{
  #     polynames <- sapply(1:input$numInputs, function(i){
  #       c(paste(input$group_pre, 'group', i, sep='_'))})
  #   }
  #   return(polynames)
  # })
  
  group_dims <- reactive({
    req(grouped_samps())
    tl <- sapply(1:input$numInputs, function(i){
      sapply(grouped_samps()[[i]], length)
    })
    sample_num <- c(tl)
    return(sample_num)
  }) 
  
  reorder_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = acc_nums()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    return(exp2)
  })

  gene_explore_plot = reactive({
    if (input$testme) {

    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    orig_df = acc_full()
    #"Acc", "Gene.Family", "SpeciesNumber", "Species"
    
    plot_df = reorder_mat()
    plot_df$Acc = orig_df$Acc
    plot_df$Gene.Family = orig_df$Gene.Family
    plot_df$SpeciesNumber = orig_df$SpeciesNumber
    plot_df$Species = orig_df$Species
    
    #gfam_list = input$acc_list
    #heyo = plot_df[grep(paste(gfam_list, collapse="|"), acc_select()), ]
    
    heyo = plot_df
    
    col = (ncol(heyo))-4
    row = nrow(heyo)
    
    RelExp = data.frame(heyo[1:(col)])
    RelExp2 = t(RelExp)
    df_RelExp = data.frame(RelExp2)
    
    datas = c()
    for (i in 1:row){
      data=as.vector(t(df_RelExp[i]))
      datas = c(datas,data)
    }
    #class(datas) = "numeric"
    
    groupings = group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles2 = rep(group_titles, row)
    
    Accs = as.vector(t(heyo[1+col]))
    Acc = rep(Accs, each = col)
    gfams = as.vector(t(heyo[2 + col]))
    Gfam = rep(gfams, each = col)
    sample_num = paste(1:(col))
    samp=strtoi(sample_num)
    isamp = rep(samp,row)
    g_data = data.frame(datas, Acc, Gfam, isamp, group_titles2)
    colnames(g_data) = c("Relative_Expression","Acc", "Gene_Family", "Sample_num", "Group")
    #g_data without zeros
    
    g_data2 = subset(g_data, Relative_Expression > 0, select=c(Relative_Expression,Acc,Gene_Family,Sample_num, Group))
    
    uniq_gfam_num = length(unique(g_data2$Gene_Family))
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #
    
    if (uniq_gfam_num > 74) {
      col_vector_edit = sample(col_vector, uniq_gfam_num, TRUE)
    } else {
      col_vector_edit = col_vector
    }
    
    
    exp_plot = ggplot(g_data2, aes(x = Group, y = Relative_Expression, fill = Gene_Family)) +
      #geom_bar(position = 'fill', stat = 'identity') +
      stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
      scale_fill_manual(values = col_vector_edit) +
      theme(legend.position = "none")
    
    # myText = paste("Group = ", g_data3$groups, '\n', "AverageExp = ", g_data3$averageExp,
    #                '\n', "Gene Family = ", g_data3$Gfam)
    # pp=plotly_build(exp_plot)   
    # style(pp, text=myText, hoverinfo = "text")
    # print(str(pp$data[[1]]$hoverinfo))
    # plotly_build(pp)
    # ggplotly(tooltip = c('x+y'))
  })
  
  output$gene_explore <- renderPlotly({
    gene_explore_plot()
  })
  
  output$gene_explore_download <- downloadHandler(
    filename = function() { paste("gene_explore", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = gene_explore_plot() +
               theme(axis.title.y = element_text(size = 22),
                     axis.title.x = element_text(size = 22),
                     axis.text.y = element_text(size = 22),
                     axis.text.x = element_text(size = 22)),
             device = 'png', 
             width = 30, height = 24, units = "cm")
    }
  )
  
  output$instructions <- renderText({
    # validate(
    #   need(input$file1, "") %then%
    #   need(unlist(grouped_samps()), "")
    #   )
    instruction = c("Hover over bars to view Gene Family Name")
    instruction
  })
  
  exp_plot <- reactive({
    if (input$testme) {
      # updateSelectizeInput(session,'acc_list', choices = acc_select(), label = "",
      #                      selected = c(" Acetyltransferase ", " GCN5_related_N_acetyltransferase "))
      # validate(
      #   need(input$acc_list, "Select at least one Gene Family!"))
      }
    else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab")) #%then%
        #need(input$acc_list, 'Select at least one Gene Family!'))
      }
    validate(
      need(input$acc_list, 'Select at least one Gene Family!')
    )
    
    orig_df = acc_full()
    #"Acc", "Gene.Family", "SpeciesNumber", "Species"
    plot_df = reorder_mat()
    plot_df$Acc = orig_df$Acc
    plot_df$Gene.Family = orig_df$Gene.Family
    plot_df$SpeciesNumber = orig_df$SpeciesNumber
    plot_df$Species = orig_df$Species
    
    gfam_list = input$acc_list
    heyo = plot_df[grep(paste(gfam_list, collapse="|"), acc_select()), ]

    col = (ncol(heyo))-4
    row = nrow(heyo)
    
    RelExp = data.frame(heyo[1:(col)])
    RelExp2 = t(RelExp)
    df_RelExp = data.frame(RelExp2)
    
    datas = c()
    for (i in 1:row){
      data=as.vector(t(df_RelExp[i]))
      datas = c(datas,data)
    }
    #class(datas) = "numeric"
    
    groupings = group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles2 = rep(group_titles, row)
    
    Accs = as.vector(t(heyo[1+col]))
    Acc = rep(Accs, each = col)
    gfams = as.vector(t(heyo[2 + col]))
    Gfam = rep(gfams, each = col)
    sample_num = paste(1:(col))
    samp=strtoi(sample_num)
    isamp = rep(samp,row)
    g_data = data.frame(datas, Acc, Gfam, isamp, group_titles2)
    colnames(g_data) = c("Exp","Acc", "Gfam", "Sample_num", "groups")
    #g_data without zeros
    
  
    #g_data2 = subset(g_data, Exp > 0, select=c(Exp,Acc,Gfam,Sample_num, groups))
    g_data2 = subset(g_data, select=c(Exp,Acc,Gfam,Sample_num, groups))
    g_data2
  })
  expression_plot <- reactive({
    g_data2 = data.frame(exp_plot())
    class(g_data2$Exp) = "numeric"
    gmax = max(g_data2$Exp)
    gmin = min(g_data2$Exp)
    
    #print(g_data2$Exp)
    
    med = median(gmax,gmin)
    #seq(gmin,gmax,1e-5)
    break_points = seq(gmin,gmax,1e-4)
    #print(g_data2)
    #pdf("Expression_Levels_box.pdf")
    g_data2$Gfam2 <- strtrim(g_data2$Gfam, 25)
    
    if (input$xy_switch){
      plot_exp = ggplot(g_data2, aes(x = Gfam2, y = Exp, fill = groups)) +
        geom_boxplot(outlier.shape=3) +
        ggtitle("Logarthimic Gene Expression") +
        #scale_colour_gradientn(colours = rainbow(length(Gfam))) +
        #scale_colour_brewer(palette="rainbow") +
        scale_y_log10() +
        xlab("Gene Family") +
        ylab("Log Relative Expression") +
        #theme(plot.margin = unit(1, "cm")) +
        #theme(plot.margin = unit(c(1, 1, 1, 2), "cm")) +
        #geom_jitter(color = 'red') +
        guides(color=FALSE, fill = guide_legend(title = "Group")) +
        theme(#axis.text.x = element_blank(), 
          #plot.margin = unit(c(.1, .1, .1, .1), "cm"),
          plot.title = element_text(hjust = 0, size = 22),
          axis.title.y = element_text(size = 22),
          axis.title.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          axis.text.x = element_text(size = 15, angle = 15, hjust=1),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18))
      plot_exp
    } else {
      plot_exp = ggplot(g_data2, aes(x = groups, y = Exp, fill = Gfam2)) +
        geom_boxplot(outlier.shape=3) +
        ggtitle("Logarthimic Gene Expression Across Relevant Samples") +
        scale_colour_gradientn(colours = rainbow(length(g_data2$Gfam))) +
        #scale_colour_brewer(palette="rainbow") +
        #scale_y_log10(limits = c(gmin,gmax)) +
        scale_y_log10() +
        xlab("Groups") +
        ylab("Log Relative Expression") +
        #theme(plot.margin = unit(1, "cm")) +
        #theme(plot.margin = unit(c(1, 1, 1, 2), "cm")) +
        #geom_jitter(color = 'red') +
        guides(color=FALSE, fill = guide_legend(title = "Gene Family")) +
        theme(#axis.text.x = element_blank(), 
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          plot.title = element_text(hjust = 0, size = 22),
          axis.title.y = element_text(size = 22),
          axis.title.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18))
      plot_exp
    }
  })
  
  output$plot1 <- renderPlot({
    expression_plot()
  })
  
  output$expression_download <- downloadHandler(
    filename = function() { paste("gene_family_expression", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = expression_plot(), device = 'png', 
             width = 60, height = 24, units = "cm")
    }
  )
  
  expression_table <- reactive({
    #req(exp_plot())
    g_data2 = data.frame(exp_plot())
    g_data2 = g_data2[order(g_data2$groups),]
    
    Gfam_uniq = unique(g_data2$Gfam)
    group_uniq = as.character(unique(g_data2$groups))
    #print(group_uniq)
    if (length(group_uniq > 1)) {
      hio=c()
      #dim = c()
      for (i in Gfam_uniq){
        acc_split = g_data2[grep(i, g_data2$Gfam),]
        acc_split[c("Exp")][is.na(acc_split[c("Exp")])] <- 0
        hi = pairwise.t.test(acc_split$Exp, acc_split$groups, p.adjust = 'BH')
        hiP <- hi$p.value
        if (length(group_uniq > 3)) {
          cutoff <- length(group_uniq) -1
          hiP <- hiP[1:cutoff-1,-1]
          hiP[upper.tri(hiP)] = t(hiP)[upper.tri(hiP)]
          hiP <- suppressWarnings(cbind(hi$p.value[cutoff,], hiP))
          hiP <- rbind(hiP, hi$p.value[cutoff,])
          rownames(hiP) <- rownames(hi$p.value)
          colnames(hiP) <- colnames(hi$p.value)
          hiP <- data.frame(hiP)
        }
        hio = rbind.fill(hio, hiP)
        #dim = c(dim, unlist(dimnames(hi$p.value)[1]))
      }
      #print(hio)
      hi3 = data.frame(rep(Gfam_uniq, each = length(group_uniq)-1))
      colnames(hi3) = "Gfam"
      # if (input$testme){
      #   hi3$group_comparisons = c(rep(group_uniq[1], length(Gfam_uniq)))
      # } else {
      #   hi3$group_comparisons = c(rep(group_uniq[-1], length(Gfam_uniq)))
      # }
      hi3$group_comparisons = c(rep(group_uniq[-1], length(Gfam_uniq)))
      hi3 = cbind.fill(hi3, hio)
      hi3
    }
  })
  
  output$exp_table = renderTable({
    expression_table()
  }, digits = 5)
  
  output$expression_table_download <- downloadHandler(
    filename = function() { paste("gene_family_expression_table", '.txt', sep='') },
    content = function(file) {
      exp = expression_table()
      exp2 = data.frame(exp)
      write.table(exp2, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  gfam_full <- reactive({
    if (input$testme){
      args = "univ_input_top.tsv"
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
      args = c(input$file1)
    }
    command ="python"
    gfam_script="univ_input_species.py"
    
    
    # Add path to script as first arg
    gfamArgs = c(gfam_script, args)
    
    gfam_data = system2(command, args=gfamArgs, stdout=TRUE)
    gfam_data
    gfam_dat = gfam_data[-1]
    gfam_da = gfam_dat[-1]
    gfams  = read.table(text = gfam_da, sep = '\t')
    #samples = paste("sample", 1:82, sep='')
    #colnames(gfams) = c("Acc", "Gene.Family", "Species", "Full", samples)
    colnames(gfams) = c("Acc", "Gene.Family", "SpeciesNumber", "Species", sample_ids())
    gfams
  })
  gfam_select <- reactive({
    if (input$testme){
      args = "univ_input_top.tsv"
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
      args = c(input$file1)
    }
    command ="python"
    gfam_script="univ_input_species.py"
    
    # Add path to script as first arg
    gfamArgs = c(gfam_script, args)
    
    gfam_data = system2(command, args=gfamArgs, stdout=TRUE)
    gfam_data
    gfam_dat = gfam_data[-1]
    gfam_da = gfam_dat[-1]
    gfams  = read.table(text = gfam_da, sep = '\t')
    #samples = paste("sample", 1:82, sep='')
    colnames(gfams) = c("Acc", "Gene.Family", "Species", "Full", sample_ids())
    gfams$Gene.Family = paste(gfams$Gene.Family, " ", sep = "")
    gfams$Gene.Family
  })
  
  reorder_spec_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = gfam_full()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    return(exp2)
  })
  
  taxa_explore_plot <- reactive({
    validate(
      need(input$taxaSep, "Please provide a Character Delimiter")
    )
    if (input$testme) {
      
    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    gfam_df = data.frame(gfam_full())
    
    spec_plot_df = reorder_spec_mat()
    spec_plot_df$Acc = gfam_df$Acc
    spec_plot_df$Gene.Family = gfam_df$Gene.Family
    spec_plot_df$SpeciesNumber = gfam_df$SpeciesNumber
    spec_plot_df$Species = gfam_df$Species
    
    #specs = subset(spec_plot_df, SpeciesNumber != "unclassified")
    specs = spec_plot_df
    
    spec_col = (ncol(specs))-4
    spec_row = nrow(specs)
    
    spec_RelExp = data.frame(specs[1:(spec_col)])
    spec_RelExp2 = t(spec_RelExp)
    df_spec_RelExp = data.frame(spec_RelExp2)
    
    spec_datas = c()
    for (i in 1:spec_row) {
      spec_data=as.vector(t(df_spec_RelExp[i]))
      spec_datas = c(spec_datas,spec_data)
    }
    
    groupings = group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles2 = rep(group_titles, spec_row)
    
    spec_Accs = as.vector(t(specs[1+spec_col]))
    spec_Acc = rep(spec_Accs, each = spec_col)
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)
    
    spec_specs = as.vector(t(specs[3+spec_col]))
    spec_spec = rep(spec_specs, each = spec_col)
    
    spec_gfams0 = as.vector(t(specs[2+spec_col]))
    #spec_gfams1 = gsub("Probable_","",spec_gfams0)
    #spec_gfams2 = tolower(spec_gfams1)
    
    spec_gfam = rep(spec_gfams0, each = spec_col)
    
    spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Relative_Expression","Acc", "Gfam", "Sample_num", "Group", "Species")
    
    #print(spec_g_data$Species)
    taxa_sub_pattern = paste0("\\", input$taxaSep, ".*")
    spec_g_data$Genus = unlist(gsub(taxa_sub_pattern,"",as.character(spec_g_data$Species))) ## what does this do?
    uniq_genus_num = length(unique(spec_g_data$Genus))
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #
    
    exp_plot = ggplot(spec_g_data, aes(x = Group, y = Relative_Expression, fill = Genus)) +
      #geom_bar(position = 'fill', stat = 'identity')
      stat_summary(fun.y = "mean", geom = "bar", position = "fill")
    if (uniq_genus_num > 74) {
      exp_plot = exp_plot +
        scale_fill_manual(values = sample(col_vector, uniq_genus_num, TRUE))
    } else {
      exp_plot = exp_plot +
        scale_fill_manual(values = col_vector)
    }
    exp_plot
    # exp2 = ggplotly(exp_plot) #make the plot a ggplot object then print it and save?
    # exp2
  })
  
  output$taxa_explore <- renderPlotly({
    taxa_explore_plot() + theme(legend.text = element_text(size = 10))
  })
  
  taxa_explore_legend <- reactive({
    legend1 <- g_legend(taxa_explore_plot())
    legend1
  })
  
  output$taxa_download <- downloadHandler(
    filename = function() { paste("taxa_explore", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = taxa_explore_plot() + 
               theme(legend.position = 'none',
                     axis.title.y = element_text(size = 22),
                     axis.title.x = element_text(size = 22),
                     axis.text.y = element_text(size = 22),
                     axis.text.x = element_text(size = 22)),
             device = 'png', 
             width = 30, height = 24, units = "cm")
    }
  )
  
  output$taxa_legend_download <- downloadHandler(
    filename = function() { paste("taxa_explore_legend", '.png', sep='') },
    content = function(file) {
      png(file, width = 60, height = 15, units ='cm', res = 300)
      par(mfrow = c(1,1))
      grid.draw(taxa_explore_legend())
      dev.off()
    }
  )
  
  species_explore_plot <- reactive({
    if (input$testme) {
      
    } else {
      validate(
        need(input$file1, "") %then%
        need(unlist(grouped_samps()), "")) 
    }
    
    gfam_df = data.frame(gfam_full())
    
    spec_plot_df = reorder_spec_mat()
    spec_plot_df$Acc = gfam_df$Acc
    spec_plot_df$Gene.Family = gfam_df$Gene.Family
    spec_plot_df$SpeciesNumber = gfam_df$SpeciesNumber
    spec_plot_df$Species = gfam_df$Species
    
    #specs = subset(spec_plot_df, SpeciesNumber != "unclassified")
    specs = spec_plot_df
    
    spec_col = (ncol(specs))-4
    spec_row = nrow(specs)
    
    spec_RelExp = data.frame(specs[1:(spec_col)])
    spec_RelExp2 = t(spec_RelExp)
    df_spec_RelExp = data.frame(spec_RelExp2)
    
    spec_datas = c()
    for (i in 1:spec_row) {
      spec_data=as.vector(t(df_spec_RelExp[i]))
      spec_datas = c(spec_datas,spec_data)
    }
    
    groupings = group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles2 = rep(group_titles, spec_row)
    
    spec_Accs = as.vector(t(specs[1+spec_col]))
    spec_Acc = rep(spec_Accs, each = spec_col)
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)
    
    spec_specs = as.vector(t(specs[3+spec_col]))
    spec_spec = rep(spec_specs, each = spec_col)
    
    spec_gfams0 = as.vector(t(specs[2+spec_col]))
    #spec_gfams1 = gsub("Probable_","",spec_gfams0)
    #spec_gfams2 = tolower(spec_gfams1)
    
    spec_gfam = rep(spec_gfams0, each = spec_col)
    
    spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Relative_Expression","Acc", "Gfam", "Sample_num", "Group", "Species")
    
    #print(spec_g_data$Species)
    spec_g_data$Genus = unlist(gsub("\\..*","",as.character(spec_g_data$Species)))
    spec_g_data$Species = unlist(gsub(".*-","",as.character(spec_g_data$Species)))
    
    uniq_spec_num = length(unique(spec_g_data$Species))
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #
    
    
    exp_plot = ggplot(spec_g_data, aes(x = Group, y = Relative_Expression, fill = Species)) +
      #geom_bar(position = 'fill', stat = 'identity')
      stat_summary(fun.y = "mean", geom = "bar", position = "fill")
    if (uniq_spec_num > 74) {
      exp_plot = exp_plot +
        scale_fill_manual(values = sample(col_vector, uniq_spec_num, TRUE))
    } else {
      exp_plot = exp_plot +
        scale_fill_manual(values = col_vector)
    }

    exp_plot
  })
  
  output$species_explore <- renderPlotly({
    species_explore_plot() + theme(legend.text = element_text(size = 10))
  })
  
  species_explore_legend <- reactive({
    legend1 <- g_legend(species_explore_plot())
    legend1
  })
  
  output$species_download <- downloadHandler(
    filename = function() { paste("species_explore", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = species_explore_plot() + 
               theme(legend.position = "none",
                     axis.title.y = element_text(size = 22),
                     axis.title.x = element_text(size = 22),
                     axis.text.y = element_text(size = 22),
                     axis.text.x = element_text(size = 22)), 
             device = 'png', 
             width = 30, height = 24, units = "cm")
    }
  )
  
  output$species_legend_download <- downloadHandler(
    filename = function() { paste("species_explore_legend", '.png', sep='') },
    content = function(file) {
      png(file, width = 60, height = 15, units ='cm', res = 300)
      grid.draw(species_explore_legend())
      dev.off()
    }
  )
  
  output$TaxaDimExp <- renderUI({
    if (input$taxaDims > 1){
      mainPanel(fluidRow(plotlyOutput("taxa_explore")),
                fluidRow(downloadButton("taxa_download", "Download Plot"),
                         downloadButton("taxa_legend_download", "Download Legend")),
                width = 12)
    }
  })
  
  output$spec <- renderPlot({
    validate(
      need(input$acc_list, 'Please Select at least one Gene Family in the Expression Tab')
    )

    spec_list = input$acc_list
    acc_column = as.vector(gfam_select())
    gfam_df = data.frame(gfam_full())
    
    spec_plot_df = reorder_spec_mat()
    spec_plot_df$Acc = gfam_df$Acc
    spec_plot_df$Gene.Family = gfam_df$Gene.Family
    spec_plot_df$SpeciesNumber = gfam_df$SpeciesNumber
    spec_plot_df$Species = gfam_df$Species
    
    specs = spec_plot_df[grep(paste(spec_list, collapse="|"), acc_column), ]
    
    spec_col = (ncol(specs))-4
    spec_row = nrow(specs)

    spec_RelExp = data.frame(specs[1:(spec_col)])
    spec_RelExp2 = t(spec_RelExp)
    df_spec_RelExp = data.frame(spec_RelExp2)

    spec_datas = c()
    for (i in 1:spec_row) {
      spec_data=as.vector(t(df_spec_RelExp[i]))
      spec_datas = c(spec_datas,spec_data)
    }
    
    groupings = group_names()
    grouping_nums = group_dims()

    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles2 = rep(group_titles, spec_row)

    spec_Accs = as.vector(t(specs[1+spec_col]))
    spec_Acc = rep(spec_Accs, each = spec_col)
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)

    spec_specs = as.vector(t(specs[3+spec_col]))
    spec_spec = rep(spec_specs, each = spec_col)

    spec_gfams0 = as.vector(t(specs[2+spec_col]))
    #spec_gfams1 = gsub("Probable_","",spec_gfams0)
    #spec_gfams2 = tolower(spec_gfams1)

    spec_gfam = rep(spec_gfams0, each = spec_col)

    spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Exp","Acc", "Gfam", "Sample_num", "Groups", "Species")
    
    #spec_g_data2 = subset(spec_g_data, Exp > 0, select=c(Exp,Acc,Gfam,Sample_num,Species))
    #head(spec_g_data2)
    #tail(spec_g_data2)
    
    #spec_g_data$Species = unlist(gsub(".*-","",as.character(spec_g_data$Species)))

    gfams_og = spec_g_data$Gfam
    gfam_og = unique(gfams_og)
    uniq_gfam = length(gfam_og)
    spec_plot_num = uniq_gfam

    library(RColorBrewer)
    n <- 60
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #

    plot_list = list()
    for (i in 1:uniq_gfam) {
      #g3 = subset(spec_g_data2, Acc == Accs[i])
      g3 = subset(spec_g_data, Gfam == gfam_og[i], drop = TRUE)
      g3 = data.frame(g3)
      uniq_spec = unique(g3$Species)
      p = ggplot(g3, aes(x=Groups, y=Exp, fill = Species)) +
        geom_bar(position = "fill", stat='identity') +
        #stat_summary(aes(g3$Groups), fun.y = "mean", geom = "bar",
        #             position = 'stack') +
        ggtitle(paste("Species Contribution for \nExpression Levels in", gfam_og[i], sep =' ')) +
        ylab("Relative Expression") +
        xlab("Group") +
        #scale_x_continuous(breaks=seq(0,spec_col,5)) +
        theme(plot.title = element_text(size = 22))
      if (length(uniq_spec) > 30) {
        p = p +
          theme(legend.position="right",
                legend.key.width = unit(0.1, "cm"),
                legend.key.height = unit(0.1, "cm"),
                legend.text = element_text(size=7),
                legend.title = element_text(size=10),
                axis.title.y = element_text(size = 22),
                axis.title.x = element_text(size = 22),
                axis.text.y = element_text(size = 22),
                axis.text.x = element_text(size = 22)) +
          guides(fill=guide_legend(nrow=30)) +
          scale_color_manual(values=sample(col_vector, 294, replace=FALSE))
        if (length(uniq_spec) > 100) {
          p = p +
            theme(legend.position="right",
                  legend.key.width = unit(0.08, "cm"),
                  legend.key.height = unit(0.08, "cm"),
                  legend.text = element_text(size=5),
                  legend.title = element_text(size=10),
                  axis.title.y = element_text(size = 22),
                  axis.title.x = element_text(size = 22),
                  axis.text.y = element_text(size = 22),
                  axis.text.x = element_text(size = 22)) +
            guides(fill=guide_legend(nrow=60)) +
            scale_color_manual(values=sample(col_vector, 294, replace=FALSE))
        }
      }
      else {
        p = p + 
          theme(legend.position="right",
                #legend.key.width = unit(0.1, "cm"),
                #legend.key.height = unit(0.1, "cm"),
                legend.text = element_text(size=10),
                legend.title = element_text(size=11),
                axis.title.y = element_text(size = 22),
                axis.title.x = element_text(size = 22),
                axis.text.y = element_text(size = 22),
                axis.text.x = element_text(size = 22)) +
          guides(fill=guide_legend(nrow=10)) +
          scale_color_manual(values=sample(col_vector, 294, replace=FALSE))
      }
      plot_list[[i]] = p
      p = NULL
    }
    multiplot(plotlist = plot_list, cols = 1)
  })
  
  spec_height = reactive({
    validate(
      need(input$acc_list, '')
    )
    
    spec_list = input$acc_list
    acc_column = as.vector(gfam_select())
    gfam_df = data.frame(gfam_full())
    specs = gfam_df[grep(paste(spec_list, collapse="|"), acc_column), ]

    spec_col = (ncol(specs))-4
    spec_row = nrow(specs)
    
    spec_RelExp = data.frame(specs[5:(spec_col+4)])
    spec_RelExp2 = t(spec_RelExp)
    df_spec_RelExp = data.frame(spec_RelExp2)
    
    spec_datas = c()
    for (i in 1:spec_row)
    {spec_data=as.vector(t(df_spec_RelExp[i]))
    spec_datas = c(spec_datas,spec_data)
    }
    
    spec_Accs = as.vector(t(specs[1]))
    spec_Acc = rep(spec_Accs, each = spec_col)
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)
    
    spec_specs = as.vector(t(specs[3]))
    spec_spec = rep(spec_specs, each = spec_col)
    
    spec_gfams0 = as.vector(t(specs[2]))
    #spec_gfams1 = gsub("Probable_","",spec_gfams0)
    #spec_gfams2 = tolower(spec_gfams1)
    
    spec_gfam = rep(spec_gfams0, each = spec_col)
    
    spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, spec_isamp, spec_spec)
    colnames(spec_g_data) = c("Exp","Acc", "Gfam", "Sample_num", "Species")
    #head(spec_g_data)
    
    #spec_g_data2 = subset(spec_g_data, Exp > 0, select=c(Exp,Acc,Gfam,Sample_num,Species))
    #head(spec_g_data2)
    #tail(spec_g_data2)
    
    gfams_og = spec_g_data$Gfam
    gfam_og = unique(gfams_og)
    uniq_gfam = length(gfam_og)
    uniq_gfam*300
  })
    
  PlotHeight = reactive(
    return(spec_plot_num*300)
  )
  output$spec_plot <- renderUI({
    if (input$testme) {
      
    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab") %then%
        need(input$acc_list, 'Please select at least one Gene Family in the Expression Tab')
      )
    }
    plotOutput("spec", height = spec_height())
  })
  
  output$expression_taxa_download <- downloadHandler(
    filename = function() { paste("gene_family_taxa", '.png', sep='') },
    content = function(file) {
      png(file, width = 40, height = 12*(spec_height()/300), units ='cm', res = 300)
      spec_list = input$acc_list
      acc_column = as.vector(gfam_select())
      gfam_df = data.frame(gfam_full())
      
      spec_plot_df = reorder_spec_mat()
      spec_plot_df$Acc = gfam_df$Acc
      spec_plot_df$Gene.Family = gfam_df$Gene.Family
      spec_plot_df$SpeciesNumber = gfam_df$SpeciesNumber
      spec_plot_df$Species = gfam_df$Species
      
      specs = spec_plot_df[grep(paste(spec_list, collapse="|"), acc_column), ]
      
      spec_col = (ncol(specs))-4
      spec_row = nrow(specs)
      
      spec_RelExp = data.frame(specs[1:(spec_col)])
      spec_RelExp2 = t(spec_RelExp)
      df_spec_RelExp = data.frame(spec_RelExp2)
      
      spec_datas = c()
      for (i in 1:spec_row) {
        spec_data=as.vector(t(df_spec_RelExp[i]))
        spec_datas = c(spec_datas,spec_data)
      }
      
      groupings = group_names()
      grouping_nums = group_dims()
      
      group_titles = c()
      for (i in 1:length(groupings)) {
        title = groupings[i]
        reps = grouping_nums[i]
        title_rep = rep(title, reps)
        group_titles = c(group_titles, title_rep)
      }
      
      group_titles2 = rep(group_titles, spec_row)
      
      spec_Accs = as.vector(t(specs[1+spec_col]))
      spec_Acc = rep(spec_Accs, each = spec_col)
      spec_sample_num = paste(1:(spec_col))
      spec_samp=strtoi(spec_sample_num)
      spec_isamp = rep(spec_samp,spec_row)
      
      spec_specs = as.vector(t(specs[3+spec_col]))
      spec_spec = rep(spec_specs, each = spec_col)
      
      spec_gfams0 = as.vector(t(specs[2+spec_col]))
      #spec_gfams1 = gsub("Probable_","",spec_gfams0)
      #spec_gfams2 = tolower(spec_gfams1)
      
      spec_gfam = rep(spec_gfams0, each = spec_col)
      
      spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
                               spec_isamp, group_titles2, spec_spec)
      colnames(spec_g_data) = c("Exp","Acc", "Gfam", "Sample_num", "Groups", "Species")
      
      #spec_g_data2 = subset(spec_g_data, Exp > 0, select=c(Exp,Acc,Gfam,Sample_num,Species))
      #head(spec_g_data2)
      #tail(spec_g_data2)
      
      spec_g_data$Species = unlist(gsub(".*-","",as.character(spec_g_data$Species)))
      
      gfams_og = spec_g_data$Gfam
      gfam_og = unique(gfams_og)
      uniq_gfam = length(gfam_og)
      spec_plot_num = uniq_gfam
      
      library(RColorBrewer)
      n <- 60
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #
      
      plot_list = list()
      for (i in 1:uniq_gfam) {
        #g3 = subset(spec_g_data2, Acc == Accs[i])
        g3 = subset(spec_g_data, Gfam == gfam_og[i], drop = TRUE)
        g3 = data.frame(g3)
        uniq_spec = unique(g3$Species)
        p = ggplot(g3, aes(x=Groups, y=Exp, fill = Species)) +
          geom_bar(position = "fill", stat='identity') +
          #stat_summary(aes(g3$Groups), fun.y = "mean", geom = "bar",
          #             position = 'stack') +
          ggtitle(paste("Species Contribution for \nExpression Levels in", gfam_og[i], sep =' ')) +
          ylab("Relative Expression") +
          xlab("Group") +
          #scale_x_continuous(breaks=seq(0,spec_col,5)) +
          theme(plot.title = element_text(size = 22))
        if (length(uniq_spec) > 30) {
          p = p +
            theme(legend.position="right",
                  legend.key.width = unit(0.1, "cm"),
                  legend.key.height = unit(0.1, "cm"),
                  legend.text = element_text(size=7),
                  legend.title = element_text(size=10),
                  axis.title.y = element_text(size = 22),
                  axis.title.x = element_text(size = 22),
                  axis.text.y = element_text(size = 22),
                  axis.text.x = element_text(size = 22)) +
            guides(fill=guide_legend(nrow=30)) +
            scale_color_manual(values=sample(col_vector, 294, replace=FALSE))
          if (length(uniq_spec) > 100) {
            p = p +
              theme(legend.position="right",
                    legend.key.width = unit(0.08, "cm"),
                    legend.key.height = unit(0.08, "cm"),
                    legend.text = element_text(size=5),
                    legend.title = element_text(size=10),
                    axis.title.y = element_text(size = 22),
                    axis.title.x = element_text(size = 22),
                    axis.text.y = element_text(size = 22),
                    axis.text.x = element_text(size = 22)) +
              guides(fill=guide_legend(nrow=60)) +
              scale_color_manual(values=sample(col_vector, 294, replace=FALSE))
          }
        }
        else {
          p = p + 
            theme(legend.position="right",
                  #legend.key.width = unit(0.1, "cm"),
                  #legend.key.height = unit(0.1, "cm"),
                  legend.text = element_text(size=10),
                  legend.title = element_text(size=11),
                  axis.title.y = element_text(size = 22),
                  axis.title.x = element_text(size = 22),
                  axis.text.y = element_text(size = 22),
                  axis.text.x = element_text(size = 22)) +
            guides(fill=guide_legend(nrow=10)) +
            scale_color_manual(values=sample(col_vector, 294, replace=FALSE))
        }
        plot_list[[i]] = p
        p = NULL
      }
      multiplot(plotlist = plot_list, cols = 1)
      dev.off()
    }
  )
  
  sig_tab <- reactive({
    gfam_DF = acc_full()
    samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=50,]
    samp_paths$Gene.Family = paste(samp_paths$Gene.Family, " ", sep="")
    samp_paths$Gene.Family
  })
  
  observe({
    if (input$testme) {
      updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
    }
    else {
      if (is.null(input$file1)) {}
      else{
        updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
      }
    }
  })
  

  actual_corr_plot <- reactive({
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    corr_list = input$sig_select

    gfam_DF1 = acc_full()
    col_num = ncol(gfam_DF1)
    row_num = nrow(gfam_DF1)
    reorder_samps = gfam_DF1[,5:col_num]
    reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
    header_DF = gfam_DF1[,1:4]
    
    gfam_DF = cbind(header_DF, reorder_samps)
    
    samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=50,]
    samp_paths$Gene.Family = paste(samp_paths$Gene.Family, " ", sep="")
    
    # groupings = group_names()
    # grouping_nums = group_dims()
    
    # group_titles = c()
    # for (i in 1:length(groupings)) {
    #   title = groupings[i]
    #   reps = grouping_nums[i]
    #   title_rep = rep(title, reps)
    #   group_titles = c(group_titles, title_rep)
    # }
    
    
    heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene.Family), ]
    
    col = (ncol(samp_paths))-4
    row = nrow(heyo)
    
    gfam = heyo$Gene.Family
    heyo_small = data.frame(gfam, heyo[5:(col+4)])
    rownames(heyo_small) = heyo$Gene.Family
    heyo_small = heyo_small[,-1]
    
    heyo_side = data.frame(t(heyo_small))
    colnames(heyo_side) = rownames(heyo_small)
    
    library(psych)
    hah2 = corr.test(heyo_side, method = "spearman")
    
    
    corr_mat = as.matrix(hah2$r)
    rownames(corr_mat) = paste(1:nrow(corr_mat))
    colnames(corr_mat) = paste(1:ncol(corr_mat))
    
    orig_p = hah2$p
    orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
    new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
    new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))
    
    new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
    diag(new_p_mat) = rep(1,nrow(new_p_mat))
    #new_p_dat = data.frame(new_p_mat)
    
    p_list = c(new_p_mat)
    p_dims = dim(new_p_mat)
    
    sym_list = c()
    for (i in p_list){
      if (is.na(i)){
        sym = ""
      }
      else {
        if (i > 0){
          sym = ""
          if (i < 0.05){
            sym = "*"
            if (i<0.01){
              sym = "**"
              if(i<0.001){
                sym = "***"
              }
            }
          }
        }
      }
      sym_list = c(sym_list, sym)
    }
    
    sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])

    library(gplots)
    my_palette <-colorRampPalette(c("blue", "white", "red"))(n=100)
    corr_mat[is.na(corr_mat)] <- 0
    v <- heatmap.2(corr_mat,
                   cellnote = sym_mat, notecol="black",
                   main="Correlation Values across \nAll Samples", #heatmap title
                   density.info="none",
                   #key = TRUE,
                   #keysize = 1.0,
                   breaks = seq(-1, 1, by = 0.02),
                   col=my_palette,
                   dendrogram="both",
                   margins=c(5,5),
                   cexRow=1,
                   cexCol=1.2,
                   trace=c("none"),
                   na.color="gray60"
    )
    v
  })
  
  output$corr_plot <- renderPlot({
    actual_corr_plot()
  })
  
  group_corr_plist <- reactive({
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    corr_list2 = input$sig_select

    groupings = group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles_uniq = unique(group_titles)
    
    corr_mat_list=list()
    for (j in 1:input$numInputs){

      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,5:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1:4]

      gfam_DF = cbind(header_DF, reorder_samps)

      samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=50,]
      samp_paths$Gene.Family = paste(samp_paths$Gene.Family, " ", sep="")

      heyo = samp_paths[grep(paste(corr_list2, collapse = '|'), samp_paths$Gene.Family), ]

      col = (ncol(samp_paths))-4
      row = nrow(heyo)

      gfam = heyo$Gene.Family
      heyo_small = data.frame(gfam, heyo[5:(col+4)])
      rownames(heyo_small) = heyo$Gene.Family
      heyo_small = heyo_small[,-1]

      heyo_side = data.frame(t(heyo_small))
      colnames(heyo_side) = rownames(heyo_small)
      
      heyo_side$Group = group_titles
      heyo_side_real = subset(heyo_side, Group == group_titles_uniq[j], select = -c(Group))
      

      library(psych)
      hah2 = corr.test(heyo_side_real, method = "spearman")

      corr_mat = as.matrix(hah2$r)
      rownames(corr_mat) = paste(1:nrow(corr_mat))
      colnames(corr_mat) = paste(1:ncol(corr_mat))

      corr_mat_list[[j]] <- corr_mat
    }
    corr_mat_list
  })

  group_sym_plist <- reactive({
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    corr_list2 = input$sig_select

    groupings = group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    group_titles_uniq = unique(group_titles)
    
    sym_mat_list=list()
    for (j in 1:input$numInputs){
      
      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,5:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1:4]
      
      gfam_DF = cbind(header_DF, reorder_samps)
      
      samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=50,]
      samp_paths$Gene.Family = paste(samp_paths$Gene.Family, " ", sep="")
      
      heyo = samp_paths[grep(paste(corr_list2, collapse = '|'), samp_paths$Gene.Family), ]
      
      col = (ncol(samp_paths))-4
      row = nrow(heyo)
      
      gfam = heyo$Gene.Family
      heyo_small = data.frame(gfam, heyo[5:(col+4)])
      rownames(heyo_small) = heyo$Gene.Family
      heyo_small = heyo_small[,-1]
      
      heyo_side = data.frame(t(heyo_small))
      colnames(heyo_side) = rownames(heyo_small)
      
      heyo_side$Group = group_titles
      heyo_side_real = subset(heyo_side, Group == group_titles_uniq[j], select = -c(Group))
      
      
      library(psych)
      hah2 = corr.test(heyo_side_real, method = "spearman")
      
      
      corr_mat = as.matrix(hah2$r)
      rownames(corr_mat) = paste(1:nrow(corr_mat))
      colnames(corr_mat) = paste(1:ncol(corr_mat))
      
      orig_p = hah2$p
      orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
      new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
      new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))
      
      new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
      diag(new_p_mat) = rep(1,nrow(new_p_mat))
      #new_p_dat = data.frame(new_p_mat)
      
      p_list = c(new_p_mat)
      p_dims = dim(new_p_mat)
      
      sym_list = c()
      for (i in p_list){
        if (is.na(i)){
          sym = ""
        }
        else {
          if (i > 0){
            sym = ""
            if (i < 0.05){
              sym = "*"
              if (i<0.01){
                sym = "**"
                if(i<0.001){
                  sym = "***"
                }
              }
            }
          }
        }
        sym_list = c(sym_list, sym)
      }
      
      sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
      sym_mat_list[[j]] <- sym_mat
    }
    sym_mat_list
  })
  
  # output$group_corrs <- renderPlot({
  #   grid.newpage()
  #
  #   lay <- grid.layout(nrow = 2, ncol=1)
  #   pushViewport(viewport(layout = lay))
  #   grid.draw(editGrob(group_corr_plist()[[1]], vp=viewport(layout.pos.row = 1,
  #                                          layout.pos.col = 1, clip=TRUE)))
  #   grid.draw(editGrob(group_corr_plist()[[2]], vp=viewport(layout.pos.row = 2,
  #                                          layout.pos.col = 1, clip=TRUE)))
  #})

  # output$group_corrs <- renderUI({
  #   plot_output_list <- lapply(1:input$numInputs, function(i) {
  #     plotname <- paste("plot", i, sep="")
  #     plotOutput(plotname, height = 280, width = 250)
  #   })
  #   # Convert the list to a tagList - this is necessary for the list of items
  #   # to display properly.
  #   do.call(tagList, plot_output_list)
  # })

  output$sig_message <- renderText({
    paste("* = p < 0.05", "** = p < 0.01", "*** = p < 0.001", sep='\n')
  })
  
  corr_label_table <- reactive({
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "") %then%
        need(unlist(grouped_samps()), ""))
    }
    
    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    
    Gene_Families = (input$sig_select)
    nums = length(Gene_Families)
    Label = paste(1:nums)
    Label_Matrix = data.frame(Gene_Families, Label)
    Label_Matrix
  })
  
  output$corr_labels <- renderTable({
    corr_label_table()
  })
  
  output$corr_table_download <- downloadHandler(
    filename = function() { paste("gene_family_expression_table", '.txt', sep='') },
    content = function(file) {
      corr_tab = corr_label_table()
      corr_tab2 = data.frame(corr_tab)
      write.table(corr_tab2, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  observe({
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "") %then%
        need(unlist(grouped_samps()), ""))
    }
    
    validate(
      need(length(input$sig_select) > 1, '')
    )
    if (input$numInputs > 1) {
    checker = input$sig_select
    output$group_corrs <- renderUI({ get_plot_output_list(max_plots, 
                                                          input$numInputs, 
                                                          group_corr_plist(),
                                                          group_sym_plist(),
                                                          group_names())
      })
    }
  })
  
  output$corr_download <- downloadHandler(
    filename = function() { paste("gene_family_correlation", '.png', sep='') },
    content = function(file) {
      png(file, width = 25, height = 20, units ='cm', res = 300)
      corr_list = input$sig_select
      
      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,5:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1:4]
      
      gfam_DF = cbind(header_DF, reorder_samps)
      
      samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=50,]
      samp_paths$Gene.Family = paste(samp_paths$Gene.Family, " ", sep="")
      
      # groupings = group_names()
      # grouping_nums = group_dims()
      
      # group_titles = c()
      # for (i in 1:length(groupings)) {
      #   title = groupings[i]
      #   reps = grouping_nums[i]
      #   title_rep = rep(title, reps)
      #   group_titles = c(group_titles, title_rep)
      # }
      
      
      heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene.Family), ]
      
      col = (ncol(samp_paths))-4
      row = nrow(heyo)
      
      gfam = heyo$Gene.Family
      heyo_small = data.frame(gfam, heyo[5:(col+4)])
      rownames(heyo_small) = heyo$Gene.Family
      heyo_small = heyo_small[,-1]
      
      heyo_side = data.frame(t(heyo_small))
      colnames(heyo_side) = rownames(heyo_small)
      
      library(psych)
      hah2 = corr.test(heyo_side, method = "spearman")
      
      
      corr_mat = as.matrix(hah2$r)
      rownames(corr_mat) = paste(1:nrow(corr_mat))
      colnames(corr_mat) = paste(1:ncol(corr_mat))
      
      orig_p = hah2$p
      orig_p[col(orig_p) == row(orig_p) | upper.tri(orig_p)] <- NA
      new_ps = p.adjust(orig_p, method = "BH", n = length(orig_p))
      new_p_mat = matrix(new_ps, ncol = ncol(orig_p), nrow = nrow(orig_p))
      
      new_p_mat[upper.tri(new_p_mat)] = t(new_p_mat)[upper.tri(new_p_mat)]
      diag(new_p_mat) = rep(1,nrow(new_p_mat))
      #new_p_dat = data.frame(new_p_mat)
      
      p_list = c(new_p_mat)
      p_dims = dim(new_p_mat)
      
      sym_list = c()
      for (i in p_list){
        if (is.na(i)){
          sym = ""
        }
        else {
          if (i > 0){
            sym = ""
            if (i < 0.05){
              sym = "*"
              if (i<0.01){
                sym = "**"
                if(i<0.001){
                  sym = "***"
                }
              }
            }
          }
        }
        sym_list = c(sym_list, sym)
      }
      
      sym_mat = matrix(sym_list, ncol = p_dims[1], nrow = p_dims[2])
      
      library(gplots)
      my_palette <-colorRampPalette(c("blue", "white", "red"))(n=100)
      corr_mat[is.na(corr_mat)] <- 0
      v <- heatmap.2(corr_mat,
                     cellnote = sym_mat, notecol="black",
                     main="Correlation Values across \nAll Samples", #heatmap title
                     density.info="none",
                     #key = TRUE,
                     #keysize = 1.0,
                     breaks = seq(-1, 1, by = 0.02),
                     col=my_palette,
                     dendrogram="both",
                     margins=c(5,5),
                     cexRow=1,
                     cexCol=1.2,
                     trace=c("none"),
                     na.color="gray60"
      )
      v
      dev.off()
    }
  )
  

  output$group_download <- renderUI({
    lapply(1:input$numInputs, function(i) {
      display_name = group_names()
      downloadButton(paste0("downloadData", i), paste("Download", display_name[i], sep=" "))
    })
  })
  
  observe({
    lapply(1:input$numInputs, function(i) {
      output[[paste0("downloadData", i)]] <- downloadHandler(
        filename = function() { paste(group_names()[i], "_correlation", '.png', sep='') },
        content = function(file) {
          png(file, width = 25, height = 20, units ='cm', res = 300)
          download_plot_output_list(max_plots,
                                    1,
                                    group_corr_plist()[i],
                                    group_sym_plist()[i],
                                    group_names()[i])
          dev.off()
        }
      )
    })
  })
  
  
}

shinyApp(ui, server)
