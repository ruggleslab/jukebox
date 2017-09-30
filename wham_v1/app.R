##R Shiny gfam app

## update 9.23 taxa search with legends

library(shiny)
library(ggplot2)
library(plyr)
library(rowr)
library(plotly)
library(grid)
library(base)
library(rPython)
library(shinythemes)
library(data.table)
library(stringr)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #

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

get_spec_plot_list <- function(max_plots, input_n, exprDF, supplied_col_vector) {
  # Insert plot output objects the list
  plot_output_list <- lapply(1:input_n, function(i) {
    plotname <- paste("plot", i, sep="")
    savename <- paste("save", i, sep="")
    plot_output_object <- plotlyOutput(plotname, height = 280, width = 250)
    plot_output_object <- renderPlotly({
      spec_df <- exprDF
      gfam_uniq <- unique(spec_df$Gfam)
      spec <- subset(spec_df, Gfam == gfam_uniq[i])
      
      spec$Species <- as.character(spec$Species)
      spec$Species <- str_trunc(spec$Species, 50, "center")

      p = ggplot(spec, aes(x=Groups, y=Exp, fill = Species)) +
        geom_bar(position = "fill", stat='identity') +
        ggtitle(paste("Species Contribution for", gfam_uniq[i], sep =' ')) +
        ylab("Relative Expression") +
        xlab("Group") +
        theme(plot.title = element_text(size = 10),
              legend.position="right",
              #legend.key.width = unit(0.1, "cm"),
              #legend.key.height = unit(0.1, "cm"),
              legend.text = element_text(size=10),
              legend.title = element_text(size=11),
              axis.title.y = element_text(size = 11),
              axis.title.x = element_text(size = 11),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10)) +
        scale_fill_manual(values = supplied_col_vector)
    })
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
}

download_spec_plot_list <- function(max_plots, i, exprDF, supplied_col_vector) {
  spec <- exprDF
  spec$Species <- as.character(spec$Species)
  spec$Species <- str_trunc(spec$Species, 50, "center")
  p = ggplot(spec, aes(x=Groups, y=Exp, fill = Species)) +
    geom_bar(position = "fill", stat='identity') +
    ggtitle(paste("Species Contribution for \n", spec$Gfam[1], sep =' ')) +
    ylab("Relative Expression") +
    xlab("Group") +
    guides(fill = guide_legend(nrow = 40))+
    theme(plot.title = element_text(size = 10),
          legend.position="right",
          #legend.key.width = unit(0.1, "cm"),
          #legend.key.height = unit(0.1, "cm"),
          legend.text = element_text(size=10),
          legend.title = element_text(size=11),
          axis.title.y = element_text(size = 11),
          axis.title.x = element_text(size = 11),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10)) +
    scale_fill_manual(values = supplied_col_vector)
  # })
  # 
  # do.call(tagList, plot_output_list) # needed to display properly.
  # 
  # return(plot_output_list)
  p
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
                                  fluidRow(column(7, 
                                                  fluidRow(column(12, h3("Methods Workflow"),
                                                                  mainPanel(img(src='WHAM_workflow.png', height = 842, width = 650))))),
                                           column(5,
                                                  fluidRow(column(12, h3("Resources")),
                                                           mainPanel(textOutput("resource_text1"),
                                                           uiOutput("samp_url"),
                                                           htmlOutput("resource_text2"),
                                                           uiOutput("hmp_url"))),
                                                  fluidRow(column(12, h3("External Information")),
                                                           mainPanel(textOutput("ext_text"),
                                                                     uiOutput("ruggles_url"),
                                                                     uiOutput("github"),
                                                                     uiOutput("github_url"),
                                                                     uiOutput("citation")))))
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
                                                   choiceNames = c("Wham"),
                                                   choiceValues = c("Wham")),
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
                                            textInput("taxaSep", "If more than one level, what character seperates your taxa levels? Permitted characters include { . , _ }")),
                                   fluidRow(uiOutput("TaxaDimExp")),
                                   fluidRow(plotlyOutput("species_explore")),
                                   fluidRow(downloadButton("species_download", "Download Plot"), 
                                            downloadButton("species_legend_download", "Download Legend")),
                                   width = 12)
                                 ),
                        tabPanel("Explore Genes",
                                 textOutput("instructions"),
                                 selectizeInput("excluder", choices=NULL,
                                                label = "Select genes to exclude",
                                                multiple = TRUE),
                                 fluidRow(plotlyOutput("gene_explore")),
                                 fluidRow(downloadButton("gene_explore_download", "Download Plot"))
                                 ),
                        tabPanel("Gene Search",
                                 fluidRow(column(4,
                                                 selectizeInput('acc_list', choices=NULL,
                                                         label = h3("Begin by selecting gene families of interest"),
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
                        tabPanel("Taxa Search",
                                 h3("Taxa contribution will appear here based on selections in the Gene Search tab"),
                                 fluidRow(column(12, uiOutput("spec_plot"))),
                                 fluidPage(uiOutput("search_taxa_download"))
                        ),
                        tabPanel("Correlation", selectizeInput('sig_select', choices=NULL,
                                                               label = h3("Begin by selecting two gene families of interest"),
                                                               multiple = TRUE),
                                 fluidPage(column(4, verbatimTextOutput("sig_message"))),
                                 fluidPage(column(10, plotOutput("corr_plot"))),
                                 fluidRow(downloadButton("corr_download", "Download Plot")),
                                 fluidPage(column(10, uiOutput("group_corrs"))),
                                 fluidPage(uiOutput("group_download")),
                                 fluidPage(column(10, tableOutput("corr_labels"))),
                                 fluidPage(downloadButton("corr_table_download", "Download Labels"))),
                        tabPanel(title = loadingLogo("https://www.youtube.com/watch?v=pIgZ7gMze7A", "wham_logo_trans.png",
                                                            'wham_grey_inf.gif', height = 135, width = 280)),
                 tags$head(tags$style('.navbar {
                            font-size: 18px}', '.navbar-brand {font-size:32px}')),
                 tags$head(tags$style("*{ font-family: Helvetica; }"))
                                                            
)

### control size of input file
options(shiny.maxRequestSize=500*1024^2)

`%then%` <- shiny:::`%OR%`

server <- function(input, output, session) {
  # file_test <- reactive({
  #   if (input$testme){
  #     #full_file = fread("univ_input.tsv", header=TRUE, sep="\t")
  #   }
  #   else {
  #     validate(
  #       need(input$file1, 'Please upload a file!')
  #     )
  #     args = c(input$file1$datapath)
  #   }
  #   result
  # })
  
  output$resource_text1 <- renderText({
    message <- c("For a sample WHAM input file please visit the link below.")
    message
  })

  url1 <- a("Sample Input File", href="https://github.com/ruggleslab/jukebox/blob/master/wham_v1/univ_input_top.tsv")
  output$samp_url <- renderUI({
    tagList("", url1)
  })
  
  
  output$resource_text2 <- renderText({
    paste("<br>", "Sample Data was derived from the HMP (Human Microbiome Project Consortium (2012) Structure, function and diversity of the healthy human microbiome. Nature, 486, 207–214.)")
  })
  
  url2 <- a("HMP HomePage", href="https://hmpdacc.org/hmp/")
  output$hmp_url <- renderUI({
    tagList("", url2)
  })
  
  output$ext_text <- renderText({
    message <- c("WHAM! is an open-source project developed in the Ruggles Lab at NYU Langone Medical Center.")
  })
  
  url3 <- a("Ruggles Lab HomePage", href="http://www.ruggleslab.org/home.html")
  output$ruggles_url <- renderUI({
    tagList("", url3)
  })
  
  output$github <- renderText({
    paste("<br>", "Source Code can be found on the Ruggles Lab Github")
  })
  
  url4 <- a("Ruggles Lab Github", href="https://github.com/ruggleslab/jukebox")
  output$github_url <- renderUI({
    tagList("", url4)
  })
  
  output$citation <- renderText({
    paste("<br>", "For additional information of using WHAM! or to cite WHAM! in your work, please refer to the following paper:",
          "*Citation Info*")
  })
  
  
  full_file <- reactive({
    # file_warning = file_test()
    # if (file_warning == "PYTHON_ERROR"){
    #   message = c(rep("Invalid File Input", 50))
    #   message_mat = matrix(message, nrow = 10, ncol = 5)
    #   full_file = data.frame(message_mat)
    # }
    if (input$testme) {
      full_file <- fread("univ_input.tsv", header=TRUE, sep=input$sep)
    }
    else {
      inFile <- input$file1
      if (is.null(inFile)) {
        return(NULL)}
      full_file = fread(inFile$datapath, header=TRUE, sep=input$sep)
      correct_cols <- c("Acc", "Gene_Family", "Species")
      validate(
        need(colnames(full_file)[1:3]==correct_cols, paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab."))
      )
    }
    full_file
  })
  
  output$contents <- renderDataTable({
    full_file <- full_file()
    full_file[1:10000]
  })
  
  sample_ids <- reactive({
    if (input$testme) {
      full_file <- fread("univ_input.tsv", header=TRUE, sep=input$sep)
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
      full_file = fread(inFile$datapath, header=TRUE, sep=input$sep)
      FullFile = data.frame(full_file)
      samples = c(colnames(FullFile))
      samp_num = length(samples)
      samps = samples[4:samp_num]
      samps
    }
  })
  
  acc_full <- reactive({
    full_file <- full_file()
    DT <- data.table(full_file[,4:50])
    DT$Gene_Family = full_file$Gene_Family
    DT2 <- DT[, lapply(.SD,sum), by = "Gene_Family"]
    #samples = paste("sample", 1:82, sep='')
    #colnames(accs) = c("Acc", "Gene.Family", "SpeciesNumber", "Species", sample_ids())
    DT2
  })
  
  acc_select <- reactive({
    accs <- acc_full()
    accs$Gene_Family = paste(" ", accs$Gene_Family, " ", sep = "")
    accs$Gene_Family
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
  
  observe({
    if (input$testme) {
      updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
    }
    else {
      if(is.null(input$file1)){}
      else {
        updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
      }
    }
  })
  
  acc_nums <- reactive({
    # validate(
    #   need(acc_full())
    # )
    accs <- acc_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Gene_Family
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
        group1 = c("SRR532024_Abundance-RPKs", "SRR532015_Abundance-RPKs", "SRR532006_Abundance-RPKs",
                   "SRR532040_Abundance-RPKs", "SRR532504_Abundance-RPKs", "SRR532507_Abundance-RPKs",
                   "SRR638753_Abundance-RPKs", "SRR640357_Abundance-RPKs", "SRR640340_Abundance-RPKs",
                   "SRR640452_Abundance-RPKs", "SRR640499_Abundance-RPKs", "SRR545546_Abundance-RPKs")
        group2 = c("SRR062353_Abundance-RPKs", "SRR062357_Abundance-RPKs", "SRR062301_Abundance-RPKs",
                   "SRR062276_Abundance-RPKs", "SRR1804686_Abundance-RPKs", "SRR1804628_Abundance-RPKs",
                   "SRR514191_Abundance-RPKs", "SRR514180_Abundance-RPKs", "SRR513168_Abundance-RPKs",
                   "SRR514231_Abundance-RPKs", "SRR513448_Abundance-RPKs")
        group3 = c("SRR062435_Abundance-RPKs", "SRR062441_Abundance-RPKs", "SRR062389_Abundance-RPKs",
                   "SRR062413_Abundance-RPKs", "SRR062402_Abundance-RPKs", "SRR062396_Abundance-RPKs",
                   "SRR346673_Abundance-RPKs", "SRR346681_Abundance-RPKs", "SRR062371_Abundance-RPKs",
                   "SRR062372_Abundance-RPKs", "SRR062462_Abundance-RPKs", "SRR062415_Abundance-RPKs")
        group4 = c("SRR528423_Abundance-RPKs", "SRR528353_Abundance-RPKs", "SRR528300_Abundance-RPKs",
                   "SRR528261_Abundance-RPKs", "SRR528183_Abundance-RPKs", "SRR528155_Abundance-RPKs",
                   "SRR532178_Abundance-RPKs", "SRR532183_Abundance-RPKs", "SRR532190_Abundance-RPKs",
                   "SRR532191_Abundance-RPKs", "SRR533152_Abundance-RPKs", "SRR533153_Abundance-RPKs")
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
      group1 = c("SRR532024_Abundance-RPKs", "SRR532015_Abundance-RPKs", "SRR532006_Abundance-RPKs",
                 "SRR532040_Abundance-RPKs", "SRR532504_Abundance-RPKs", "SRR532507_Abundance-RPKs",
                 "SRR638753_Abundance-RPKs", "SRR640357_Abundance-RPKs", "SRR640340_Abundance-RPKs",
                 "SRR640452_Abundance-RPKs", "SRR640499_Abundance-RPKs", "SRR545546_Abundance-RPKs")
      group2 = c("SRR062353_Abundance-RPKs", "SRR062357_Abundance-RPKs", "SRR062301_Abundance-RPKs",
                 "SRR062276_Abundance-RPKs", "SRR1804686_Abundance-RPKs", "SRR1804628_Abundance-RPKs",
                 "SRR514191_Abundance-RPKs", "SRR514180_Abundance-RPKs", "SRR513168_Abundance-RPKs",
                 "SRR514231_Abundance-RPKs", "SRR513448_Abundance-RPKs")
      group3 = c("SRR062435_Abundance-RPKs", "SRR062441_Abundance-RPKs", "SRR062389_Abundance-RPKs",
                 "SRR062413_Abundance-RPKs", "SRR062402_Abundance-RPKs", "SRR062396_Abundance-RPKs",
                 "SRR346673_Abundance-RPKs", "SRR346681_Abundance-RPKs", "SRR062371_Abundance-RPKs",
                 "SRR062372_Abundance-RPKs", "SRR062462_Abundance-RPKs", "SRR062415_Abundance-RPKs")
      group4 = c("SRR528423_Abundance-RPKs", "SRR528353_Abundance-RPKs", "SRR528300_Abundance-RPKs",
                 "SRR528261_Abundance-RPKs", "SRR528183_Abundance-RPKs", "SRR528155_Abundance-RPKs",
                 "SRR532178_Abundance-RPKs", "SRR532183_Abundance-RPKs", "SRR532190_Abundance-RPKs",
                 "SRR532191_Abundance-RPKs", "SRR533152_Abundance-RPKs", "SRR533153_Abundance-RPKs")
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
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    orig_df = acc_full()
    #"Acc", "Gene.Family", "SpeciesNumber", "Species"
    
    plot_df = reorder_mat()
    plot_df$Gene_Family = orig_df$Gene_Family
    
    plot_df$Sum = rowSums(plot_df[,1:ncol(plot_df)-1])
    plot_df = plot_df[order(-plot_df$Sum),]
    
    if (nrow(plot_df) < 500){
      out_bound <- nrow(plot_df)
    }
    else {
      out_bound = 500
    }
    plot_df2 <- plot_df[1:out_bound,]
    #gfam_list = input$acc_list
    #heyo = plot_df[grep(paste(gfam_list, collapse="|"), acc_select()), ]
    
    heyo = plot_df2[,1:ncol(plot_df2)-1]

    col = (ncol(heyo))-1
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
    
    gfams = as.vector(t(heyo[1 + col]))
    Gfam = rep(gfams, each = col)
    Gfam = paste0(" ",Gfam," ")
    sample_num = paste(1:(col))
    samp=strtoi(sample_num)
    isamp = rep(samp,row)
    g_data = data.frame(datas, Gfam, isamp, group_titles2)
    colnames(g_data) = c("Relative_Expression", "Gene_Family", "Sample_num", "Group")
    #g_data without zeros
    
    g_data2 = subset(g_data, Relative_Expression > 0, select=c(Relative_Expression,Gene_Family,Sample_num, Group))
    
    if (length(input$excluder)>0){
      exclude_list = input$excluder
      g_data2 = g_data2[grep(paste(exclude_list, collapse="|"), g_data2$Gene_Family, invert=TRUE), ]
    }
    
    #gfam_list = input$acc_list
    #heyo = plot_df[grep(paste(gfam_list, collapse="|"), acc_select()), ]
    
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
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
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
    plot_df$Gene_Family = paste(" ", orig_df$Gene_Family, " ", sep="")
    
    gfam_list = input$acc_list
    heyo = plot_df[grep(paste(gfam_list, collapse="|"), acc_select()), ]

    col = (ncol(heyo))-1
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
    
    gfams = as.vector(t(heyo[1 + col]))
    Gfam = rep(gfams, each = col)
    #Gfam = paste0(" ", Gfam, " ")
    sample_num = paste(1:(col))
    samp=strtoi(sample_num)
    isamp = rep(samp,row)
    g_data = data.frame(datas, Gfam, isamp, group_titles2)
    colnames(g_data) = c("Exp", "Gfam", "Sample_num", "groups")
    #g_data without zeros
    
  
    #g_data2 = subset(g_data, Exp > 0, select=c(Exp,Acc,Gfam,Sample_num, groups))
    g_data2 = subset(g_data, select=c(Exp,Gfam,Sample_num, groups))
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
    g_data2$Gfam2 <- gsub(" ", "", g_data2$Gfam)
    g_data2$Gfam2 <- str_trunc(g_data2$Gfam2, 35, "center")
    #g_data2$Gfam2 <- strtrim(g_data2$Gfam, 25)
    
    if (input$xy_switch){
      plot_exp = ggplot(g_data2, aes(x = Gfam2, y = Exp, fill = groups)) +
        geom_boxplot(outlier.shape=3) +
        ggtitle("Logarthimic Gene Expression") +
        scale_fill_manual(values = col_vector) +
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
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 22),
          axis.text.x = element_text(size = 15, angle = 15, hjust=1),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18))
      plot_exp
    } else {
      plot_exp = ggplot(g_data2, aes(x = groups, y = Exp, fill = Gfam2)) +
        geom_boxplot(outlier.shape=3) +
        ggtitle("Logarthimic Gene Expression Across Relevant Samples") +
        scale_fill_manual(values = col_vector) +
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
    if (length(group_uniq) > 1) {
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
      #hi3
    }
    else if (length(group_uniq)==1){
      hi3 <- data.frame(Seletion=1:length(Gfam_uniq), Gene_Family=Gfam_uniq)
      hi3
    }
    hi3
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
    }
    else {
      validate(
        need(input$file1, 'Please upload a file!')
      )
    }
    gfams <- full_file()
    #samples = paste("sample", 1:82, sep='')
    #colnames(gfams) = c("Acc", "Gene.Family", "Species", "Full", samples)
    #colnames(gfams) = c("Acc", "Gene_Family", "Taxa", sample_ids())
    gfams
  })

  gfam_select <- reactive({
    validate(
      need(gfam_full(), '')
      )
    gfams <- gfam_full()
    gfams$Gene_Family = paste(" ", gfams$Gene_Family, " ", sep = "")
    gfams$Gene_Family
  })

  reorder_spec_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = gfam_full()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    return(exp2)
  })
  
  spec_genus_data <- reactive({
    validate(
      need(input$taxaSep, "Please provide a Character Delimiter")
    )
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    gfam_df = data.frame(gfam_full())
    
    spec_plot_df = reorder_spec_mat()
    
    DT <- data.table(spec_plot_df)
    DT$Species = gfam_df$Species
    DT2 <- DT[, lapply(.SD,sum), by = "Species"]
    
    specs = DT2
  
    spec_col = (ncol(specs))-1
    spec_row = nrow(specs)
    
    spec_RelExp = data.frame(specs[,-1])
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
    
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)
    
    spec_spec <- rep(specs$Species, each = spec_col)
    
    spec_g_data = data.frame(spec_datas, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Relative_Expression", "Sample_num", "Group", "Species")
    
    #print(spec_g_data$Species)
    taxa_sub_pattern = paste0("\\", input$taxaSep, ".*")
    spec_g_data$Genus = unlist(gsub(taxa_sub_pattern,"",as.character(spec_g_data$Species))) 
    
    spec_g_data$Genus <- str_trunc(spec_g_data$Genus, 70, "center")
    spec_g_data
    
  })
  
  taxa_explore_plot <- reactive({
    spec_g_data <- spec_genus_data()
    uniq_genus_num = length(unique(spec_g_data$Genus))
    
    exp_plot = ggplot(spec_g_data, aes(x = Group, y = Relative_Expression, fill = Genus)) +
      #geom_bar(position = 'fill', stat = 'identity')
      stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
      scale_fill_manual(values = sample(col_vector, uniq_genus_num, TRUE))
    exp_plot
    # exp2 = ggplotly(exp_plot) #make the plot a ggplot object then print it and save?
    # exp2
  })
  
  output$taxa_explore <- renderPlotly({
    taxa_explore_plot() + theme(legend.text = element_text(size = 10))
  })
  
  taxa_explore_legend <- reactive({
    legend1 <- g_legend(taxa_explore_plot() + guides(fill=guide_legend(ncol=3)))
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
      genera <- length(unique(spec_genus_data()$Genus))
      png(file, width = 45, height = 0.25*genera, units ='cm', res = 300)
      par(mfrow = c(1,1))
      grid.draw(taxa_explore_legend())
      dev.off()
    }
  )
  
  spec_species_data <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab")) 
    }
    
    gfam_df = data.frame(gfam_full())
    
    spec_plot_df = reorder_spec_mat()
    
    DT <- data.table(spec_plot_df)
    DT$Species = gfam_df$Species
    DT2 <- DT[, lapply(.SD,sum), by = "Species"]
    
    specs = DT2
    
    spec_col = (ncol(specs))-1
    spec_row = nrow(specs)
    
    spec_RelExp = data.frame(specs[,-1])
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
    
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)
    
    spec_spec <- rep(specs$Species, each = spec_col)
    
    spec_g_data = data.frame(spec_datas, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Relative_Expression", "Sample_num", "Group", "Species")
    
    #print(spec_g_data$Species)
    spec_g_data$Genus = unlist(gsub("\\..*","",as.character(spec_g_data$Species)))
    spec_g_data$Species = unlist(gsub(".*-","",as.character(spec_g_data$Species)))
    
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #
    
    spec_g_data$Species <- str_trunc(spec_g_data$Species, 70, "center")
    spec_g_data
    
  })
  
  species_explore_plot <- reactive({
    spec_g_data <- spec_species_data()
    uniq_spec_num = length(unique(spec_g_data$Species))
    
    exp_plot = ggplot(spec_g_data, aes(x = Group, y = Relative_Expression, fill = Species)) +
      #geom_bar(position = 'fill', stat = 'identity')
      stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
      scale_fill_manual(values = sample(col_vector, uniq_spec_num, TRUE))
    exp_plot
  })
  
  output$species_explore <- renderPlotly({
    species_explore_plot() + theme(legend.text = element_text(size = 8))
  })
  
  species_explore_legend <- reactive({
    legend1 <- g_legend(species_explore_plot()+guides(fill=guide_legend(ncol=3)))
    legend1
  })
  
  output$species_download <- downloadHandler(
    filename = function() { paste("species_explore", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_species_data()$Species))
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
      species = length(unique(spec_species_data()$Species))
      png(file, width = 45, height = 0.25*species, units ='cm', res = 300)
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
  
  spec <- reactive({
    validate(
      need(input$acc_list, 'Please Select at least one Gene Family in the Expression Tab')
    )

    spec_list = input$acc_list
    acc_column = as.vector(gfam_select())
    gfam_df = data.frame(gfam_full())
    
    spec_plot_df = reorder_spec_mat()
    spec_plot_df$Acc = gfam_df$Acc
    spec_plot_df$Gene_Family = gfam_df$Gene_Family
    spec_plot_df$Species = gfam_df$Species
    
    specs = spec_plot_df[grep(paste(spec_list, collapse="|"), acc_column), ]
    
    spec_col = (ncol(specs))-3
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
    spec_gfam = paste0(" ", spec_gfam, " ")

    spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Exp","Acc", "Gfam", "Sample_num", "Groups", "Species")
    spec_g_data
  })
  
  spec_colors <- reactive({
    spec_df <- spec()
    most_specs <- max(table(spec_df$Gfam))
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" ##manually change color here if you don't like what you see in the figs! #
    
    cols <- sample(col_vector, most_specs, replace=TRUE)
    cols
  })
  
  observe({
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
          need(unlist(grouped_samps()), "Please select groups in the Group Tab"))
    }
    output$spec_plot <- renderUI({
      validate(
        need(length(input$acc_list) > 0, 'Please select at least one gene family on the previous tab')
      )
      get_spec_plot_list(max_plots, 
                         length(input$acc_list), 
                         spec(),
                         spec_colors())
      })
  })
  
  output$search_taxa_download <- renderUI({
    lapply(1:length(input$acc_list), function(i) {
      display_name = input$acc_list
      downloadButton(paste0("downloadTaxa", i), paste("Download", display_name[i], sep=" "))
    })
  })
  
  observe({
    lapply(1:length(input$acc_list), function(i) {
      output[[paste0("downloadTaxa", i)]] <- downloadHandler(
        filename = function() { paste(input$acc_list[i], "_taxa", '.png', sep='') },
        content = function(file) {
          spec <- spec()
          spec_curr <- subset(spec, spec$Gfam==input$acc_list[i])
          specs_nums <- length(unique(spec_curr$Species))
          png(file, width = 25+0.3*specs_nums, height = 25, units ='cm', res = 300)
          plotter<-download_spec_plot_list(max_plots, 
                                  1,
                                  subset(spec, spec$Gfam==input$acc_list[i]),
                                  spec_colors())
          print(plotter)
          dev.off()
        }
      )
    })
  })
    
  sig_tab <- reactive({
    gfam_DF = acc_full()
    samp_paths = gfam_DF
    #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
    samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
    samp_paths$Gene_Family
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
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    }
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
    reorder_samps = gfam_DF1[,-1]
    reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
    header_DF = gfam_DF1[,1]
    
    gfam_DF = cbind(header_DF, reorder_samps)
    samp_paths = gfam_DF
    #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
    samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
    
    
    heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene_Family), ]
    
    col = (ncol(samp_paths))-1
    row = nrow(heyo)
    
    heyo_small <- heyo
    rownames(heyo_small) = heyo$Gene_Family
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
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]

      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
      samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")

      heyo = samp_paths[grep(paste(corr_list2, collapse = '|'), samp_paths$Gene_Family), ]

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small<-heyo
      rownames(heyo_small) = heyo$Gene_Family
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
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]
      
      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
      samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
      
      heyo = samp_paths[grep(paste(corr_list2, collapse = '|'), samp_paths$Gene_Family), ]
      
      col = (ncol(samp_paths))-1
      row = nrow(heyo)
      
      heyo_small <- heyo
      rownames(heyo_small) = heyo$Gene_Family
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
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "")
      )
    }
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
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]
      
      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
      samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
      
      # groupings = group_names()
      # grouping_nums = group_dims()
      
      # group_titles = c()
      # for (i in 1:length(groupings)) {
      #   title = groupings[i]
      #   reps = grouping_nums[i]
      #   title_rep = rep(title, reps)
      #   group_titles = c(group_titles, title_rep)
      # }
      
      
      heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene_Family), ]
      
      col = (ncol(samp_paths))-1
      row = nrow(heyo)
      
      heyo_small<-heyo
      rownames(heyo_small) = heyo$Gene_Family
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
