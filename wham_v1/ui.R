library(shiny)
library(shinythemes)

ui <- navbarPage(title = "Workflow Hub for Automated Metagenomic Exploration",
                 theme = shinytheme("superhero"),
                       tabPanel("Home",
                                  fluidRow(column(7, 
                                                  wellPanel(h3("Methods Workflow"),
                                                            img(src='Figure_1.png', height = 647.7, width = 500))),
                                           column(5,
                                                  mainPanel(h3("Resources"),
                                                            textOutput("resource_text1"),
                                                            uiOutput("samp_url"),
                                                            htmlOutput("resource_text2"),
                                                            uiOutput("hmp_url"),
                                                            htmlOutput("humann2wham_text"),
                                                            uiOutput("humann2wham_url"),
                                                            h3("External Information"),
                                                            textOutput("ext_text"),
                                                            uiOutput("ruggles_url"),
                                                            uiOutput("github"),
                                                            uiOutput("github_url"),
                                                            uiOutput("citation"))))
                                ),
                       tabPanel("Upload",
                                fluidPage(
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput('file1', 'Choose TSV File (Max=200MB)',
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
                                      numericInput("filter_level", "Enter Filter Level",
                                                   min = 0, max = 1, value = 0.9),
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
                                   fluidRow(numericInput("taxaDims", 
                                                         "How many taxa levels 
                                                         are present in your 
                                                         data?", 1, min = 1, max = 2, width = '35%'),
                                            textInput("taxaSep", 
                                                      "If more than one level, 
                                                      what character seperates 
                                                      your taxa levels? 
                                                      Permitted characters include 
                                                      { . , _ }", width = '35%')),
                                   fluidRow(verbatimTextOutput("ex_delimiter")),
                                   fluidRow(numericInput("upper_limit",
                                                         "Taxa upper limit",
                                                         1, min = 0, max = 1, width = '35%'),
                                            numericInput("lower_limit",
                                                         "Taxa lower limit",
                                                         0, min=0, max=1, width = '35%')),
                                   fluidRow(uiOutput("TaxaDimExp")),
                                   fluidRow(plotlyOutput("species_explore")),
                                   fluidRow(downloadButton("species_download", "Download Plot"), 
                                            downloadButton("species_legend_download", "Download Legend"),
                                            downloadButton("species_raw_data", "Download Table")),
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
                                 fluidRow(column(12,plotOutput("plot1", height = '600px'))),
                                 fluidRow(downloadButton("expression_download", "Download Plot")),
                                 h5("Pairwise T-test results are performed here if at least 2 groups are selected"),
                                 fluidRow(column(8,uiOutput("exp_heat"))),
                                 fluidRow(uiOutput("exp_heat_download")),
                                 fluidRow(column(6, align= "center", tableOutput("exp_table"))),
                                 fluidRow(downloadButton("expression_table_download", "Download Table"))
                                 ),
                        tabPanel("Taxa Search",
                                 h5("Taxa contribution will appear here based on selections in the Gene Search tab"),
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