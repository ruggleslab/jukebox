#when publsihing bioconductor source error
#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")

### control size of input file
options(shiny.maxRequestSize=200*1024^2)

`%then%` <- shiny:::`%OR%`

server <- function(input, output, session) {
  
  #### Home Page & Resource Information
  output$resource_text1 <- renderText({
    message <- c("For a sample WHAM input file please visit the link below.")
    message
  })

  url1 <- a("Sample Input File", href="https://github.com/ruggleslab/jukebox/blob/master/wham_v1/sample_input.tsv.zip", target="_blank")
  output$samp_url <- renderUI({
    tagList("", url1)
  })
  
  
  output$resource_text2 <- renderText({
    paste("<br>", "Sample Data was derived from the HMP (Human Microbiome Project Consortium (2012) Structure, function and diversity of the healthy human microbiome. Nature, 486, 207–214.)")
  })
  
  url2 <- a("HMP HomePage", href="https://hmpdacc.org/hmp/", target="_blank")
  output$hmp_url <- renderUI({
    tagList("", url2)
  })
  
  output$humann2wham_text <- renderText({
    paste("<br>", "HUMANn2 users can readily convert HUMANn2 gene tables 
          to a wham-friendly input using our collection of python conversion scripts,
          available on the Ruggles Lab Github.")
  })
  
  url3 <- a("HUMANn2 File Converter", href="https://github.com/ruggleslab/jukebox/tree/master/wham_v1/file_conversion_scripts", target="_blank")
  output$humann2wham_url <- renderUI({
    tagList("", url3)
  })
  
   
  output$ext_text <- renderText({
    message <- c("WHAM! is an open-source project developed in the Ruggles Lab at NYU Langone Medical Center.")
  })
  
  url4 <- a("Ruggles Lab HomePage", href="http://www.ruggleslab.org/home.html", target="_blank")
  output$ruggles_url <- renderUI({
    tagList("", url4)
  })
  
  output$github <- renderText({
    paste("<br>", "Source Code can be found on the Ruggles Lab Github")
  })
  
  url5 <- a("Ruggles Lab Github", href="https://github.com/ruggleslab/jukebox", target="_blank")
  output$github_url <- renderUI({
    tagList("", url5)
  })
  
  output$citation <- renderText({
    paste("<br>", "For additional information of using WHAM! or to cite WHAM! in your work, please refer to the following paper:",
          "https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4870-z")
  })
  
  observe({
    if (input$input_type == "Biobakery"){
      output$file_selector <- renderUI({
        fileInput('file1', 'Choose TSV File (Max=200MB)',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
      })
    }
    if (input$input_type == "EBI"){
      output$file_selector <- renderUI({
        fluidRow(
          fileInput('features1', 'Choose Feature File (Max=200MB)',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
          fileInput('taxa1', 'Choose Taxa File (Max=200MB)',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
        )
      })
    }
  })
  
  
  ## Load in the data or try a sample dataset
  full_file_feature <- reactive({
    if (input$testme) {
      start_col = 4
      #full_file <- fread("sample_input2.tsv", header=TRUE, sep=input$sep)
      full_file_feature <- fread("sample_input.tsv", header=TRUE, sep='\t')
      #full_file_feature <- fread("antibiotic_count_est_u90_wham.tsv", header=TRUE, sep='\t') ##temp test
      m_count <<- c("")
    }
    else {
      if (input$input_type == "Biobakery"){
        start_col = 4
        inFile <- input$file1
        if (is.null(inFile)) {
          return(NULL)}
        full_file <- try(
          {fread(inFile$datapath, header=TRUE, sep='\t')})
        
        correct_cols <- c("Acc", "Feature", "Taxa")
        colnames(full_file)[1:3] = correct_cols
        validate(
          need(class(full_file)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
            need(all(unlist(lapply(full_file[,1:3], class))=="character"), paste0("The input file provided is not an appropriate format. Please be sure the first three columns are characters.")) %then%
            need(length(colnames(full_file))>4, "Please provide more than 1 sample")
        )
        
        nums <- data.matrix(full_file[,start_col:ncol(full_file)])
        rownames(nums) <- rownames(full_file)
        
        cols <- colnames(nums)
        checker <- apply(nums, 1, is.integer)
        if (!all(checker==T)){
          num_check <- mean(colSums(nums))
          if(num_check < 100000){
            nums <- nums*100000
            nums <- round(nums,0)
            nums <- t(apply(nums, 1, as.integer))
            colnames(nums) <- cols
            m_count <<- c("Input not a count matrix, values converted to integers and scaled to 1million where applicable")
          } else {
            nums <- round(nums,0)
            nums <- t(apply(nums, 1, as.integer))
            colnames(nums) <- cols
            }
        } else {m_count <<- c("")}
        checker <- apply(nums, 1, is.integer)
        validate(
          need(all(checker==T), "It appears the provided input is not a count matrix!")
        )
        full_file <- cbind(full_file[,1:3], nums)
        full_file <- subset(full_file, Feature != "NO_NAME")
        full_file$Feature <- gsub("[^[:alnum:]']", "_", full_file$Feature)
        full_file_feature <- full_file
      }
      if (input$input_type == "EBI"){
        start_col = 3
        inFile_taxa <- input$taxa1
        inFile_feature <- input$features1
        
        if (is.null(inFile_feature)) {
          return(NULL)}
        full_file_feature <- try(
          {fread(inFile_feature$datapath, header=TRUE, sep='\t')})
        
        correct_cols <- c("Acc", "Feature")
        colnames(full_file_feature)[1:2] <- correct_cols
        validate(
          need(class(full_file_feature)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
            need(all(unlist(lapply(full_file_feature[,1:2], class))=="character"), paste0("The input file provided is not an appropriate format. Please be sure the first two columns are characters."))%then%
            need(length(colnames(full_file_feature))>3, "Please provide more than 1 sample")
        )
        
        nums <- data.matrix(full_file_feature[,start_col:ncol(full_file_feature)])
        
        checker <- apply(nums, 1, is.integer)
        validate(
          need(all(checker==T), "It appears the provided input 
               is either not a count matrix or another inapppropriate format")
        )
  
        full_file_feature <- subset(full_file_feature, Feature != "NO_NAME")
        full_file_feature$Feature <- gsub("[^[:alnum:]']", "_", full_file_feature$Feature)
        full_file_feature
      }
    }
    full_file_feature
  })
  
  full_file_taxa <- reactive({
    ## now taxa
    if (input$input_type == "EBI"){
      start_col = 2
      inFile_taxa <- input$taxa1

      if (is.null(inFile_taxa)) {
        return(NULL)}
    full_file <- try(
      {fread(inFile_taxa$datapath, header=TRUE, sep='\t')})
    
    correct_cols <- c("Taxa")
    colnames(full_file)[1] <- correct_cols
    validate(
      need(class(full_file)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
        need(unlist(lapply(full_file[,1], class) =="character"), paste0("The input file provided is not an appropriate format. Please be sure the first column is characters."))
    )
    
    full_file <- subset(full_file, Taxa != "NO_NAME")
    full_file$Taxa <- gsub(";k__", "_k__", full_file$Taxa)
    full_file
    }
  })
  
  # Preview Full File
  output$contents <- DT::renderDataTable({
    validate(
      need(full_file_feature(),"")
    )
    full_file <- full_file_feature()
    if (nrow(full_file)<50){
      if(ncol(full_file)<10){
        full_file_show <- full_file[,1:ncol(full_file)]
      }
      else{
        full_file_show <- full_file[,1:10]
      }
    }
    else{
      if(ncol(full_file)<10){
        full_file_show <- full_file[1:50, 1:ncol(full_file)]
      } else {
        full_file_show <- full_file[1:50, 1:10]
      }
    }
    datatable(full_file_show) %>% formatStyle(T, color="white", background="#2c3e4f")
  })
  
  output$contents_taxa <- DT::renderDataTable({
    validate(
      need(full_file_taxa(),"")
    )
    full_file <- full_file_taxa()
    if (nrow(full_file)<50){
      full_file_show <- full_file[,1:ncol(full_file)]
    }
    else{
      full_file_show <- full_file[1:50, 1:ncol(full_file)]
    }
    datatable(full_file_show) %>% formatStyle(T, color="white", background="#2c3e4f")
  })
  
  output$temp_info1 <- renderUI({
    if (input$input_type == "Biobakery"){
      m1 <- c("Biobakery format requires an output file from the Biobakery HUMANn2 Tool")
      m2 <- c("'Features' in the below format can be Gene Families, Pathways or GO Terms")
      m3 <- c("For formatting help please see the sample input on the Home Tab")
      m4 <- c("Please ensure the first three columns are labeled 'Acc', 'Feature', 'Taxa'")
      
      HTML(paste('<br/>', m1, m2, m3, m4, '<br/>', sep = '<br/>'))
    } else if (input$input_type == "EBI"){
      m1 <- c("EBI format requires 2 inputs") 
      m2 <- c("Feature input requires an output file from the EBI Metagenomics Pipeline")
      m3 <- c("'Features' in the below format are typically Interpro Identifiers")
      m4 <- c("Please ensure the first two columns are labeled 'Acc', 'Feature'")
      
      HTML(paste('<br/>', m1, m2, m3, m4, '<br/>', sep = '<br/>'))
    }
  })
  
  output$temp_info2 <- renderUI({
    m5 <- c("Taxa input requires the EBI taxonomic output in the below format")
    m6 <- c("Please ensure the first column is labeled 'Taxa'")
    
    HTML(paste('<br/>', m5, m6, sep = '<br/>'))
  })
  
  output$temp_df <- DT::renderDataTable({
    if (input$input_type == "Biobakery"){
      Acc <- paste0("0000", 1:5)
      Feature <- paste0("Feature", 1:5)
      Taxa <- c(paste0("Bug", 1:3), paste0("Bug", 1:2))
      num_mat <- matrix(NA,5,5)
      for (i in 1:ncol(num_mat)){
        num_mat[,i] <- sample(1000,5)
      }
      colnames(num_mat) = paste0("Sample", 1:5)
      sample_df <- cbind(Acc, Feature, Taxa, num_mat)
    }
    else if (input$input_type == "EBI"){
      Acc <- paste0("0000", 1:5)
      Feature <- paste0("Feature", 1:5)
      num_mat <- matrix(NA,5,5)
      for (i in 1:ncol(num_mat)){
        num_mat[,i] <- sample(1000,5)
      }
      colnames(num_mat) = paste0("Sample", 1:5)
      sample_df <- cbind(Acc, Feature, num_mat)
    }
    datatable(sample_df) %>% formatStyle(T, color="white", background="#2c3e4f")
  })
  
  output$temp_taxa <- renderDataTable({
    Taxa <- c(paste0("Bug", 1:3), paste0("Bug", 1:2))
    num_mat <- matrix(NA,5,5)
    for (i in 1:ncol(num_mat)){
      num_mat[,i] <- sample(1000,5)
    }
    colnames(num_mat) = paste0("Sample", 1:5)
    sample_df <- cbind(Taxa, num_mat)
    datatable(sample_df) %>% formatStyle(T, color="white", background="#2c3e4f")
  })
  
  observe({
    if (input$input_type == "Biobakery"){
      if (input$testme){
        output$preview_shower <- renderUI({
          dataTableOutput('contents')
        })
      }
      else if(!is.null(input$file1)){
        output$preview_shower <- renderUI({
          dataTableOutput('contents')
        })
      }
      else if(is.null(input$file1)){
        output$preview_shower <- renderUI({
          fluidPage(fluidRow(
            uiOutput('temp_info1'),
            tags$head(tags$style(
              "#temp_info1{color: #df691a; font-size: 18px}")),
            dataTableOutput('temp_df')
          ))
        })
      }
    }
    if (input$input_type == "EBI"){
      if (!is.null(input$features1)){
        output$preview_shower <- renderUI({
          fluidRow(
            dataTableOutput('contents'),
            dataTableOutput('contents_taxa'))
        })
      }
      else{
        output$preview_shower <- renderUI({
          fluidPage(fluidRow(
            uiOutput('temp_info1'),
            tags$head(tags$style(
              "#temp_info1{color: #df691a; font-size: 18px}")),
            dataTableOutput('temp_df'),
            uiOutput('temp_info2'),
            tags$head(tags$style(
              "#temp_info2{color: #df691a; font-size: 18px}")),
            dataTableOutput('temp_taxa')
          ))
        })
      }
    }
  })
  
  observe({
    if (input$testme){
      updateRadioButtons(session, "input_type", selected = "Biobakery")
    }
  })
  
  ### Generate prerequisites for plotting

  # Generate dataframe collapsed by Gene Family
  acc_pre_filt <- reactive({
    full_file <- full_file_feature()
    col_num <- ncol(full_file)
    if (input$testme){
      start_col = 4
      
      DT <- data.frame(full_file[,start_col:col_num], check.names = F)
      
      DT2 <- aggregate(DT, list(full_file$Feature), sum)
      colnames(DT2)[1] = "Feature"
      DT2 = subset(DT2, Feature != "filler")
      
      DT2
    } else {
      if (input$input_type == "Biobakery"){
        start_col = 4
        
        DT <- data.frame(full_file[,start_col:col_num], check.names = F)
        DT2 <- aggregate(DT, list(full_file$Feature), sum)
        colnames(DT2)[1] = "Feature"
        DT2 = subset(DT2, Feature != "filler")

        DT2
      }
      if (input$input_type == "EBI"){
        start_col = 3
        
        DT <- data.frame(full_file[,start_col:col_num], check.names = F)
        DT2 <- aggregate(DT, list(full_file$Feature), sum)
        colnames(DT2)[1] = "Feature"
        DT2 = subset(DT2, Feature != "filler")
        
        DT2
      }
    }
    DT2
  })
  
  acc_full <- reactive({
    DT2 <- acc_pre_filt()
    
    if (input$filter_level != 0){
      feat_var <- rowVars(as.matrix(DT2[,-1]))
      keep_quant <- quantile(feat_var, input$filter_level)
      keep_filt <- which(feat_var > keep_quant)
      DT2 <- DT2[keep_filt,]
    } else {
      DT2 <- DT2
    }
    
    DT2
  })
  
  # acc full but numeric values only
  acc_nums <- reactive({
    accs <- acc_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Feature
    accs3 = as.matrix(accs2)
    accs3
  })

  
  ## tell user how much was filtered
  
  output$filter_message1 <- renderUI({
    m1 <- c("Filtering of low variance features is highly 
        recommended in order to speed up differential abundance calculations.")
    
    HTML(paste(m1))
  })
  
  observeEvent(nrow(full_file_feature()>0), {
    output$filter_message2 <- renderUI({
      
      validate(
        need(nrow(full_file_feature()) > 0, "")
      )
      
      old <- nrow(acc_pre_filt())
      new <- nrow(acc_full())
      
      m2 <- paste0("Original number of unique features: ", old)
      m3 <- paste0("Number of unique features after filtering: ", new)
      
      
      HTML(paste('<br/>', m2, m3, '<br/>', sep = '<br/>'))
    })
  })
  
  observeEvent(nrow(full_file_feature()>0), {
    output$filter_message3 <- renderUI({
      
      validate(
        need(nrow(full_file_feature()) > 0, "")
      )
      
      HTML(paste(m_count, sep = '<br/>'))
    })
  })
  
  
  # Generate dataframe collapsed by Species
  spec_full <- reactive({
    if (input$testme){
      start_col = 4
      full_file <- full_file_feature()
      col_num <- ncol(full_file)
      
      DT <- data.frame(full_file[,start_col:col_num], check.names=F)
      DT2 <- aggregate(DT, list(full_file$Taxa), sum)
      colnames(DT2)[1] = "Taxa"
      DT2 = subset(DT2, Taxa != "filler")
      DT2
    }
    if (input$input_type == "Biobakery"){
      start_col = 4
      full_file <- full_file_feature()
      col_num <- ncol(full_file)
      
      DT <- data.frame(full_file[,start_col:col_num], check.names=F)
      DT2 <- aggregate(DT, list(full_file$Taxa), sum)
      colnames(DT2)[1] = "Taxa"
      DT2 = subset(DT2, Taxa != "filler")
      DT2
    }
    if (input$input_type == "EBI"){
      start_col = 2
      validate(
        need(full_file_taxa(), "Taxa File Uploads are required for EBI inputs")
      )
      full_file <- full_file_taxa()
      col_num <- ncol(full_file)
      
      DT <- data.frame(full_file[,start_col:col_num], check.names=F)
      DT2 <- aggregate(DT, list(full_file$Taxa), sum)
      colnames(DT2)[1] = "Taxa"
      DT2 = subset(DT2, Taxa != "filler")
      DT2
    }
    DT2
  })

  spec_nums <- reactive({
    accs <- spec_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Taxa
    accs3 = as.matrix(accs2)
    rownames(accs3) <- rownames(accs2)
    accs3
  })

  # Generate selection list of gene families
  acc_select <- reactive({
    accs <- acc_full()
    accs$Feature
  })
   
  # When file is uploaded update choice of gene families
  observe({
    if (input$testme) {
      updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
      updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
      updateNumericInput(session, 'taxaDims', value = 6, min = 6, max = 7, step = 1)
      updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
    }
    else {
      if (input$input_type == "Biobakery"){
        if(is.null(input$file1)){}
        else {
          updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
          updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
          updateNumericInput(session, 'taxaDims', value = 6, min = 6, max = 7, step = 1)
          updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
        }
      }
      if (input$input_type == "EBI"){
        if(is.null(input$features1)){}
        else {
          updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
          updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
          updateNumericInput(session, 'taxaDims', value = 5, min = 1, max = 7, step = 1)
          updateSelectizeInput(session,'sig_select', choices = sig_tab(), server = TRUE)
        }
      }
    }
  })

  ## Group Selection Elements ##


  ###input grouped samples based on number of groups, output reordered matrix###
  results <- c()
  makeReactiveBinding('results')

  # create number of  columns based on input number of groups
  observe({
    output$inputGroup = renderUI({
      if (input$testme){
        updateNumericInput(session, 'numInputs', value = 4)
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
        fluidRow(selectInput("test_groups", label = "",choices = group1,
                             multiple = TRUE, selected = group1),
                 selectInput("test_groups", label = "", choices = group2,
                             multiple = TRUE, selected = group2),
                 selectInput("test_groups", label = "", choices = group3,
                             multiple = TRUE, selected = group3),
                 selectInput("test_groups", label = "", choices = group4,
                             multiple = TRUE, selected = group4))
      } else {
        if (input$input_type == "Biobakery"){
          validate(
            need(input$file1, 'Please provide a file in the Upload Tab'))
        }
        if (input$input_type == "EBI"){
          validate(
            need(input$features1, 'Please provide a feature file in the Upload Tab'))
        }
        input_list <- lapply(1:input$numInputs, function(i) {
          inputName <- paste0("input", i)
          sampleUploadUI(inputName)
          })
      }
    })
    # Make list of resulting group allocations
    results <<- lapply(1:input$numInputs, function(i) {
      inputName <- paste0("input", i)
      callModule(sampleUpload, inputName, acc_nums)
    })
  })

  # Generate Space to Label Groups (based on input$numInputs)
  output$group_pre = renderUI({
    lapply(1:input$numInputs, function(i) {
      inputName <- paste0("group", i)
      textInput(inputName, label = "Group Name")
    })
  })


  # Ensure groups are named even if not input by user
  group_names <- reactive({
    if (input$testme){
      groups <- c("Arm", "Vagina", "Saliva", "Stool")
    }
    else{
      groups <- sapply(1:input$numInputs, function(i){
        cc<-input[[paste0("group", i)]][1]
        cc
      })
    }
    return(groups)
  })

  new_group_names <- reactive({
    groups <- sapply(1:input$numInputs, function(i){
      cc <- group_names()[i]
      if (cc==''){
        cc = paste0("Group", i)
      }
      cc
    })
    return(groups)
  })


  #Retain group allocations
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
      result_call <- try(lapply(1:input$numInputs, function(i) {
        results[[i]]()}))
      result_call
    }
  })

  ## Store dimensions of each group for plottng
  group_dims <- reactive({
    req(grouped_samps())
    tl <- sapply(1:input$numInputs, function(i){
      sapply(grouped_samps()[[i]], length)
    })
    sample_num <- c(tl)
    return(sample_num)
  })

  # Reorder input matrix based on group allocations
  reorder_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = acc_nums()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    rownames(exp2) <- acc_full()$Feature
    return(exp2)
  })

  # Generate warning if sample is assigned to multiple groups
  output$group_warning <- renderText({
    validate(
      need(length(unlist(grouped_samps()))>0, 'Select groups!')
    )
    group_list = c(unlist(grouped_samps()))
    duplicates = c(duplicated(group_list))
    dups <- group_list[which(duplicates==TRUE)]
    if ('TRUE' %in% duplicates) {
      message = paste0("Warning: ", dups, " assigned to more than one group!")
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

  #### Begin Plotting ! ####

  #
  #
  #
  #
  #
  #
  ## Ultimate ColorBar
  color_select <- reactive({
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    ##### make color bar to add to the plotly object!
    names <- bin
    uniq_names <- unique(names)
    
    cols_keep2 <- c("#f70c1c", "#6d89e8", "#91dd68", "#482e91", "#fc9207",
                    "#fcdb23", "#e87ac3", "#5beabd", "#01871e", "#a0080f")
    
    cols_cols <- cols_keep2[1:length(uniq_names)]
    cols_cols
    
  })
  
  colorer <- reactive({
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    cols_cols <- color_select()
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    ##### make color bar to add to the plotly object!
    names <- bin
    uniq_names <- unique(names)
    
    meta <- data.frame(names)
    
    for (i in 1:length(uniq_names)){
      new_bin <- as.integer(meta$names == uniq_names[i])
      meta <- cbind(meta, new_bin)
    }
    
    series_mat <- meta[,-1]
    
    if (input$numInputs > 1){
      
      colnames(series_mat) <- uniq_names
      #series_mat <- t(series_mat)
      series_mat[series_mat==0] <- NA
      
      g1 = t(as.matrix(series_mat[,1]))
      g1col <- data.frame(x = c(0,1), y = c(cols_cols[1], cols_cols[1]))
      colnames(g1col) <- NULL
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        tickcolor = 'white',
        showgrid = FALSE
      )
      
      colorer <- plot_ly(
        type = "heatmap"
      ) %>% add_trace(
        x=1:length(bin),
        z = g1,
        colorscale = g1col, showscale = F,
        hoverinfo = 'all'
      ) %>%
        layout(yaxis = ax, xaxis = ax)
      
      for (i in 2:length(uniq_names)){
        g = t(as.matrix(series_mat[,i]))
        gcol <- data.frame(x = c(0,1), y = c(cols_cols[i], cols_cols[i]))
        colnames(gcol) <- NULL
        
        colorer <- colorer %>% add_trace(
          x=1:length(bin),
          z = g,
          colorscale = gcol, showscale = F,
          hoverinfo = 'all'
        )
      }
    }
    else {
      names(series_mat) <- uniq_names
      series_mat[series_mat==0] <- NA
      
      g1 = t(as.matrix(series_mat))
      g1col <- data.frame(x = c(0,1), y = c(cols_cols[1], cols_cols[1]))
      colnames(g1col) <- NULL
      
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        tickcolor = 'white',
        showgrid = FALSE
      )
      
      colorer <- plot_ly(
        type = "heatmap"
      ) %>% add_trace(
        x=1:length(bin),
        z = g1,
        colorscale = g1col, showscale = F,
        hoverinfo = 'all'
      ) %>%
        layout(yaxis = ax, xaxis = ax)
    }
    colorer
  })
  
  
  ##### color bar legend element!
  
  colorer_key <- reactive({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 40),
            legend.title = element_text(size = 42))
    
    legend <- g_legend(gg) 
    grid.arrange(legend) 
  })
  
  ### downloadable legend!
  output$legend_download <- downloadHandler(
    filename = function() { paste("legend_download", '.png', sep='') },
    content = function(file) {
      leg_plot <- colorer_key()
      leg_length = input$numInputs
      
      png(file, height = 15, width = 30, units = 'cm', res = 300)
      grid.arrange(leg_plot)
      dev.off()
    }
  )
  
  
  output$key0_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "") %then%
        need(!is.null(unlist(grouped_samps())), "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    validate(
      need(length(Group) == length(cols_cols), "")
    )
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend) 
  }, bg="transparent")
  
  output$key0 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key0_plot", height = pheight)
  })
  
  
  output$key1_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "") %then%
        need(nrow(spec_taxa_data())> 0, "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend) 
  }, bg="transparent")
  
  output$key1 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key1_plot", height = pheight)
  })
  
  output$key2_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend) 
  }, bg="transparent")
  
  output$key2 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key2_plot", height = pheight)
  })
  
  output$key3_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key3 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key3_plot", height = pheight)
  })
  
  output$key4_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "") %then%
        need(input$numInputs > 1, "Differential Abundance not calculated because less than 2 groups are selected") %then%
        need(length(da_feat())> 0, "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key4 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key4_plot", height = pheight)
  })
  
  output$key5_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "") %then%
        need({
          event.data <- event_data("plotly_click", source = "g_exp")
          !is.null(event.data)}, "") %then%
        need(input$input_type == "Biobakery", "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key5 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key5_plot", height = pheight)
  })
  
  
  output$key6_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "") %then%
        need({
          event.data <- event_data("plotly_click", source = "g_select")
          !is.null(event.data)}, "") %then%
        need(input$input_type == "Biobakery", "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 28))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key6 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key6_plot", height = pheight)
  })
  
  
  output$key7_plot <- renderPlot({
    validate(
      need(!is.null(new_group_names()[[1]]), "") %then%
        need(length(input$acc_list) > 0, "")
    )
    groupings = new_group_names()
    cols_cols <- color_select()
    
    n <- length(groupings)
    Group = groupings
    
    plotplot <- data.frame(Group,cols_cols)
    plotplot$Group <- factor(plotplot$Group, levels = groupings)
    
    gg = ggplot(plotplot, aes(1:n, 1:n, color = Group)) +
      geom_point() +
      guides(color = guide_legend(direction = "horizontal",
                                  override.aes = list(shape = 15, size=9),
                                  ncol  = 4))+
      scale_color_manual(values = cols_cols) + 
      theme(legend.key = element_rect(fill = "transparent"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 28))
    
    legend <- g_legend(gg) 
    grid.draw(legend)
  }, bg="transparent")
  
  output$key7 <- renderUI({
    pheight = paste0('50px')
    
    if (input$numInputs > 4){
      pheight = paste0('100px')
    }
    if (input$numInputs > 8){
      pheight = paste0('150px')
    }
    
    plotOutput("key7_plot", height = pheight)
  })
  
  #
  #
  #
  #
  #
  ##
  ## Gene Plots
  ##

  # Gene Explore Heatmap!
  feat_mat <- reactive({
    
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }
    
    feat_plot_df = reorder_mat()
    feat_var2 <- feat_plot_df
    feat_var2
  })
  
  da_feat_stat <- reactive({
    feat_order <- feat_mat()

    groupings = new_group_names()
    grouping_nums = group_dims()
    
    feat_trans <- feat_order
    cc <- data.frame(feat_trans)
    
    RA <- unlist(feat_order)
    RA_clr <- unlist(cc)
    sample_id <- rep(colnames(feat_order), each = nrow(feat_order))
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    conds <- group_titles
    
    combos <- combn(unique(conds), 2)
    stat_results = data.frame(we.ep = NA, we.eBH = NA, wi.ep = NA, wi.eBH = NA,
                              rab.all = NA, rab.win.a = NA, rab.win.b = NA,
                              diff.btw = NA, diff.win = NA, effect = NA, overlap = NA,
                              Int = NA, Feature = NA)
    time_ep = NULL
    for (i in 1:ncol(combos)){
      withProgress({
        start = Sys.time()
        ab <- combos[,i]
        a <- ab[1]
        b <- ab[2]
        a_sel = which(conds == a)
        b_sel = which(conds == b)
        
        feat_subset <- cbind(feat_order[,a_sel], feat_order[,b_sel])
        cond_subset <- c(conds[a_sel], conds[b_sel])
        x <- aldex.clr(feat_subset, cond_subset, mc.samples=16, verbose=TRUE)
        x.tt <- try(aldex.ttest2(x, cond_subset, paired.test = F))
        x.effect <- try(aldex.effect(x, cond_subset, include.sample.summary=FALSE, verbose=TRUE))
        
        validate(
          need(class(x.tt)!= "try-error", 
               paste0("Differential Abundance failed due to the following error: ", x.tt))%then%
            need(class(x.effect)!= "try-error", 
                 paste0("Differential Abundance failed due to the following error: ", x.effect))
        )
        
        x.all <- data.frame(x.tt, x.effect, stringsAsFactors=FALSE)
        x.all$Int <- rep(paste0(a,"-", b), nrow(x.all))
        x.all$Feature <- rownames(x.all)
        names(x.all) = names(stat_results)
        stat_results <- rbind(stat_results, x.all)
        end <- Sys.time()
        time_ep <- paste0(": ", round(end-start,1), "s per step")
      }, message = paste0("Calculating Differential Abundance ", 
                          i, " of ", ncol(combos),
                          time_ep)) 
    }
    print(dim(stat_results))
    
    #stat_results <- da_feat_stat()
    up_effect <- max(abs(stat_results$effect), na.rm=T)
    up_pval <- max(abs(stat_results$we.eBH), na.rm=T)
    
    print(up_pval)
    print(round(up_pval+0.05, 2))
    
    print(up_effect)
    
    updateSliderInput(session, "feat_pval_up", max = round(up_pval+0.05, 2))
    updateNumericInput(session, "feat_pval_num", max = round(up_pval+0.05, 2))
    updateSliderInput(session, "feat_effect", max = ceiling(up_effect))
    updateNumericInput(session, "feat_effect_num", max = ceiling(up_effect))
    
    stat_results
  })
  
  output$feat_pval_distr <- renderPlot({
    stat_results <- da_feat_stat()
    par(mfrow=c(1,2), mar = c(1,1,1,1))
    cd <- density(stat_results$we.eBH, na.rm = T)
    plot(cd,
         xlim = c(0, 
                  round(max(stat_results$we.eBH, na.rm = T)+0.05, 2)),
         frame.plot = F, ann = F, axes = F, col = 'orange', lwd = 5)
    abline(v=input$feat_pval_up, col = 'red', lwd = 5)
    
    ce <- density(abs(stat_results$effect), na.rm = T)
    plot(ce, 
         xlim = c(0,
                  ceiling(max(abs(stat_results$effect), na.rm=T))),
         frame.plot = F, ann = F, axes = F, col = 'orange', lwd = 5)
    abline(v=input$feat_effect, col = 'dodgerblue2', lwd = 5)
  }, bg = 'transparent')
  
  output$feat_pval_scatter <- renderPlot({
    stat_results <- da_feat_stat()
    plot(stat_results$we.eBH, abs(stat_results$effect), pch = 16,
         xlab = "Adjusted p Value", ylab = "Effect Size",
         xlim = c(0, round(max(stat_results$we.eBH, na.rm = T)+0.05, 2)),
         ylim = c(0, ceiling(max(abs(stat_results$effect), na.rm=T))))
    rect(0,input$feat_effect, input$feat_pval_up, ceiling(max(abs(stat_results$effect), na.rm=T)),
         lwd = 1.5, col = NA, border = 'orange')
  })
  
  da_feat <- reactive({
    stat_results <- da_feat_stat()
    
    keepers_df <- subset(stat_results, we.eBH < input$feat_pval_up &
                               wi.eBH < input$feat_pval_up)
    keepers <- unique(subset(keepers_df, effect > input$feat_effect |
                               effect < -input$feat_effect)$Feature)
    
    print(length(keepers))
    validate(
      need(length(keepers) > 0, "No signficant features were detected!")
    )
    keepers
  })
  
  
  output$gene_da_placehold <- renderUI({
    m1 <- c("Differentially Abundant Features across all samples")
    HTML(paste(m1, '<br/>', sep = '<br/>'))
  })
  
  
  gene_da_plotly <- reactive({
    if(input$numInputs > 1){
      da_gene <- da_feat()
      feat_order <- feat_mat()
      go_show <- feat_order[da_gene,]
    } else {
      feat_order <- feat_mat()
      go_show <- feat_order
    }
    

    groupings = new_group_names()
    grouping_nums = group_dims()
    
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    go_show_nums2 <- log10(go_show + 1)
    
    gg_nums <- go_show_nums2
    
    validate(
      need(nrow(gg_nums) > 1, "At least two features required for plotting, try adjusting the sliders above!")
    )
    feat_clust <- hclust(dist(gg_nums))
    
    gg_nums_feat <- gg_nums[feat_clust$order,]
    gg_nums_samp <- gg_nums_feat
    
    gg_names <- factor(colnames(gg_nums_samp), levels = colnames(gg_nums_samp))
    gg_feat <- factor(rownames(gg_nums_samp), levels = rownames(gg_nums_samp))
      
    x <- sweep(gg_nums_samp, 1L, rowMeans(gg_nums_samp, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    
    Z_scored_CPM <- unlist(x)
    names <- rep(gg_names, each = nrow(gg_nums))
    feat <- rep(gg_feat, ncol(gg_nums))
    groups <- rep(bin, each = nrow(gg_nums))
    
    
    plotter <- data.frame(Z_scored_CPM, names, groups, feat)
    
    validate(
      need(all(is.na(plotter$Z_scored_CPM))!= T, "Must provide more than 1 unique Sample")
    )
    
    ###
    library(plotly)
    
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      tickcolor = 'white',
      showgrid = FALSE
    )
    p <- plot_ly(source = "g_exp",
                 x = plotter$names, y = plotter$feat,
                 text = plotter$groups,
                 z = plotter$Z_scored_CPM, type = "heatmap",
                 hoverinfo = 'y+text',
                 colorbar = list(x = -0.3, 
                                 xanchor = 'left',
                                 tickmode='array',
                                 tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                                 ticktext = c("low", "high"), tickfont = list(size=18))) %>%
      layout(yaxis = ax, xaxis = ax)
    
    ####
    p3 <- subplot(colorer(),
                  p, nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = TRUE)
    p3
    
  })
    
  observeEvent({
    input$feat_pval_up
  },{
    
    output$gene_da <- renderPlotly({
      plott <- gene_da_plotly()
      plott
    })
  })
  
  plot3_df <- reactive({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    validate(
      need(is.null(event.data) == F, "") %then%
        need(event.data$y %in% full_file_feature()$Feature, "") %then%
        need(input$input_type == "Biobakery", "Specific Taxa Assignments not availble for EBI inputs")
    )
    
    select_gfam <- event.data$y
    
    full_full <- full_file_feature()
    col_keep <- colnames(full_full[,4:ncol(full_full)])
    full_select <- data.frame(subset(full_full, Feature %in% select_gfam))
    colnames(full_select) <- c("Acc", "Feature", "Taxa", col_keep)


    go_nums <- full_select[,4:ncol(full_select)]
    go_reorder <- go_nums[,unlist(grouped_samps())]

    RA <- unlist(go_reorder)
    
    samp_order <- factor(colnames(go_reorder), 
                         levels = c(colnames(go_reorder)))
    Sample_id <- rep(samp_order, each = nrow(go_reorder))
    
    Taxa <- rep(full_select$Taxa, ncol(go_reorder))
  

    ## groups
    groupings = new_group_names()
    grouping_nums = group_dims()
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    Group = rep(bin, each = nrow(go_reorder))
    ###
    
    plotter <- data.frame(Sample_id, Taxa, Group, RA)
    plotter 
  })
  
  plot3_ggplot <- reactive({
    plotter <- plot3_df()
    
    event.data <- event_data("plotly_click", source = "g_exp")
    select_gfam <- event.data$y
    
    #### filter to speed up plotting!
    plotter2 <- aggregate(plotter$RA, list(plotter$Taxa), sum)
    
    summer <- sum(plotter2[,2])
    plotter2$prop <- plotter2[,2]/summer
    
    spec_up <- subset(plotter2, prop > 0.01)
    spec_down <- subset(plotter2, prop < 0.01)
    
    high_tax <- as.character(unlist(spec_up[,1]))
    low_tax <- as.character(unlist(spec_down[,1]))
    
    spec_high <- subset(plotter, Taxa %in% high_tax)
    spec_low <- subset(plotter, Taxa %in% low_tax)
    
    spec_low$Taxa <- factor(rep("Less than 1%", nrow(spec_low)))
    
    plotter2 <- rbind(spec_high, spec_low)
    
    ####
    plotter <- plotter2

    bug_sorter <- plotter[order(plotter$RA, decreasing = T),]
    bug_sort <- as.character(unique(bug_sorter$Taxa))
    plotter$Taxa <- factor(plotter$Taxa, levels = bug_sort)
    
    tax_sel <- ggplot(plotter, aes(x=Sample_id, y = RA, fill = Taxa, text = Group)) + 
      stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
      scale_fill_manual(values = randomColor(length(unique(plotter$Taxa))),
                        na.value = "grey") +
      ggtitle(paste0("Taxa Contribution: ", select_gfam)) + 
      xlab("Sample") + ylab("Relative Abundance") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 15),
            axis.title = element_text(size = 15))
    tax_sel
  })
  
  output$Plot3 <- renderPlotly({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    validate(
      need(is.null(event.data) == F, "") %then%
        need(event.data$y %in% full_file_feature()$Feature, "")
    )
    tax_sel <- plot3_ggplot()+theme(legend.position='none')
    
    taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
    taxly <- taxly %>% layout(margin=list(l = 100), yaxis=list(tickprefix=" "))
    
    
    p3 <- subplot(colorer(), taxly,
                  nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T, titleY = T)
    p3
  })

  
  
  #####
  
  da_feat_stat_mat <- reactive({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    validate(
      need(is.null(event.data) == F, "")
    )
    
    validate(
      need(event.data$y %in% da_feat_stat()$Feat, "")
    )
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    #### now select feat
    
    select_gfam <- event.data$y
    
    statter <- da_feat_stat()
    statter2 <- subset(statter, Feature %in% select_gfam)
    
    total_ints <- c()
    total_int <- combn(groupings, 2)
    for (i in 1:ncol(total_int)){
      int <- total_int[,i]
      ints <- paste0(int[1], "-", int[2])
      total_ints <- c(total_ints, ints)
    }
    
    have <- statter2$Int
    add <- setdiff(total_ints, have)
    have_ints <- c(have, add)
    
    num_fill <- round(statter2$we.eBH, 4)
    if (nrow(statter2) != ncol(combn(groupings, 2))){
      how_short <- ncol(combn(groupings, 2)) - nrow(statter2)
      add <- rep(NA, how_short)
      num_fill <- c(num_fill, add)
    }
    
    name_df <- data.frame(name=have_ints, num = num_fill)
    rownames(name_df) <- name_df$name
    name_df = name_df[total_ints,]
    
    group_num <- length(groupings)
    group1 <- gsub('-.*', "", name_df$name[1])
    group_rest <- gsub(".*-", "", name_df$name[1:(group_num-1)])
    group_names <- c(group1, group_rest)
    
    num_fill2 <- name_df$num
    
    resm <- matrix(NA, group_num, group_num)
    resm[lower.tri(resm) ] <- num_fill2
    resm <- t(resm)
    resm[lower.tri(resm) ] <- num_fill2
    rownames(resm) <- group_names
    colnames(resm) <- group_names
    print(resm)
    resm
  })
  
  
  observe({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
    } else {
      output$da_feat_stat_heat <- renderPlot({
        resm <- da_feat_stat_mat()
        
        resm_lab <- resm
        resm_lab[resm_lab < 0.00001] <- "<0.00001"
        
        breaker <- seq(0, 1, by = 0.005)
        coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
                   c(colorRampPalette(c("white"))(n=189)))
        
        select_feat <- event.data$y
        
        if (input$numInputs == 2){
          resm[1,2] = resm[1,2] + 0.0000001
        }
        if (sum(resm, na.rm = T) == 0){
          resm[1,2] = 1e-7
        }
        
        par(cex.main=0.8)
        heatmap.2(data.matrix(resm),
                  cellnote = resm_lab,
                  density.info="none",
                  trace=c("none"),
                  notecol="black",
                  key.title = NULL,
                  sepcolor="black",
                  breaks = breaker,
                  col = coler,
                  dendrogram = 'none',
                  Rowv=F,
                  Colv=F,
                  margins=c(10,10),
                  cexRow=1.2,
                  cexCol=1.2#,
        )
        title(paste0("Significance between groups for \n", select_feat), line= -2.5)
      })
    }
  })
  
  
  observe({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
    } else {
      output$da_feat_stat_tab <- renderTable({
        resm <- da_feat_stat_mat()
        resm_lab <- resm
        resm_lab[resm_lab < 0.00001] <- "<0.00001"
        resm_lab <- cbind(Names = rownames(resm_lab), resm_lab)
        resm_lab
      }, caption = {
        select_taxa <- event.data$y
        return(select_taxa)
      })
    }
  })
  
  output$feat_map_description <- renderUI({
    event.data <- event_data("plotly_click", source = "g_exp")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
      m1 <- c("Click the Heatmap for Differential Abundance Results")
      m2 <- c("Additional feature specific Taxa Information will appear if available")
      HTML(paste('<br/>', '<br/>', m1, m2, sep = '<br/>'))
    } else {
      return(NULL)
    }
  })
  
  observe({
    validate(
      need(input$numInputs > 1, "Need at least 2 groups for ALDEX Test")
    )
      output$da_feat_stat_ui <- renderUI({
        fluidPage(
          uiOutput("feat_map_description"),
          tags$head(tags$style(
            "#feat_map_description {color:#df691a; font-size:18px}")),
          plotOutput("da_feat_stat_heat"),
          fluidRow(fluidPage(downloadButton("feature_stat_download", "Download Heatmap"),
                    downloadButton("feature_stat_results", "Download Stats")))
        )
      })
  })
  
  
  ### download uis
  
  output$gene_explore_download <- downloadHandler(
    filename = function() { paste("feature_explore", '.png', sep='') },
    content = function(file) {
      p3 <- gene_da_plotly()
      
      p3$width = 1200
      p3$height = 800
      
      export(p3, file = file)
      
    })
  
  
  
  
  output$gene_explore_taxa_download <- downloadHandler(
    filename = function() { paste("taxa_explore", '.png', sep='') },
    content = function(file) {
      
      tax_sel <- plot3_ggplot()+theme(legend.position='none',
                                      axis.text.y = element_text(size = 18),
                                      axis.title = element_text(size = 18))
      
      taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
      taxly <- taxly %>% layout(margin=list(l = 100), yaxis=list(tickprefix=" "))
      
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T, titleY = T)
      p3
      
      p3$width = 1200
      p3$height = 800
    
      export(p3, file = file)
      
    }
  )
  
  gene_explore_taxa_legend <- reactive({
    tax_sel <- plot3_ggplot()
    legend1 <- g_legend(tax_sel+guides(fill=guide_legend(ncol=3)))
    legend1
  })
  
  output$gene_explore_taxa_legend_download <- downloadHandler(
    filename = function() { paste("gene_explore_taxa_legend", '.png', sep='') },
    content = function(file) {
      species = length(unique(plot3_df()$Taxa))
      png(file, width = 35, height = 0.3*(species/1), units ='cm', res = 300)
      grid.draw(gene_explore_taxa_legend())
      dev.off()
    }
  )
  
  
  
  output$feature_stat_download <- downloadHandler(
    filename = function() { paste("differential_features_heat", '.png', sep='') },
    content = function(file) {
      event.data <- event_data("plotly_click", source = "g_exp")
      resm <- da_feat_stat_mat()
      
      resm_lab <- resm
      resm_lab[resm_lab < 0.00001] <- "<0.00001"
      
      breaker <- seq(0, 1, by = 0.005)
      coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
                 c(colorRampPalette(c("white"))(n=189)))
      
      select_gfam <- event.data$y

      if (input$numInputs == 2){
        resm[1,2] = resm[1,2] + 0.0000001
      }
      if (sum(resm, na.rm = T) == 0){
        resm[1,2] = 1e-7
      }
      
      dim_fact <- dim(resm)[1]
      if (dim_fact < 6) {
        fill_size <- 2-(0.17*dim_fact)
      } else {
        fill_size <- 2-(0.19*dim_fact)
      }
      
      png(file, width = 15, height = 12, units = 'cm', res = 300)
      par(cex.main=0.8)
      v = heatmap.2(data.matrix(resm),
                    cellnote = resm_lab,
                    notecex = fill_size,
                    density.info="none",
                    trace=c("none"),
                    notecol="black",
                    key.title = NULL,
                    breaks = breaker,
                    col = coler,
                    dendrogram = 'none',
                    Rowv=F,
                    Colv=F,
                    margins=c(7,7),
                    cexRow=1.2,
                    cexCol=1.2
      )
      v
      title(paste0("Significance between groups for \n", select_gfam), line= -1.5)
      dev.off()
    }
  )
  
  output$feature_stat_results <- downloadHandler(
    filename = function() { paste("differential_features_aldex", '.txt', sep='') },
    content = function(file) {
      stat_results <- da_feat_stat()
      #keepers <- subset(stat_results, we.eBH < input$feat_pval_up &
      #                    wi.eBH < input$feat_pval_up)
      write.table(stat_results, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  
  ####### Query Your Data ######
  
  
  #### new plot1 for selectable features
  
  ####
  
  plot1_df <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }
    
    validate(
      need(length(input$acc_list) > 0, "")
    )
    
    da_gene <- input$acc_list
    feat_order <- reorder_mat()
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    go_show <- subset(feat_order, rownames(feat_order) %in% da_gene)
    
    go_show_nums2 <- log10(go_show + 1)
    
    x <- sweep(go_show_nums2, 1L, rowMeans(go_show_nums2, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    
    Z_scored_CPM <- unlist(x)
    

    names <- rep(colnames(go_show_nums2), each = nrow(go_show_nums2))
    feat <- rep(rownames(go_show_nums2), ncol(go_show_nums2))
    groups <- rep(bin, each = nrow(go_show_nums2))
    
    plotter <- data.frame(Z_scored_CPM, names, groups, feat)
    
    plotter$names <- factor(plotter$names, levels = colnames(go_show_nums2))
    
    plotter
    
  })
  
  plot1_plotly <- reactive({
      plotter <- plot1_df()

    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      tickcolor = 'white',
      showgrid = FALSE
    )
    p <- plot_ly(source = "g_select",
                 x = plotter$names, y = plotter$feat,
                 text = plotter$groups,
                 z = plotter$Z_scored_CPM, type = "heatmap",
                 hoverinfo = 'y+text',
                 colorbar = list(x = -0.3, 
                                 xanchor = 'left',
                                 tickmode='array',
                                 tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                                 ticktext = c("low", "high"), tickfont = list(size=18))) %>%
      layout(yaxis = ax, xaxis = ax)
    
    ####
    p3 <- subplot(colorer(),
                  p, nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T)
    
    p3
    
  })
  
  output$plot1 <- renderPlotly({
    plot1_plotly <- plot1_plotly()
    plot1_plotly
  })
  
  output$plot1_placehold <- renderUI({
    if (length(input$acc_list) < 1){
      m1 <- c("Select Features to view relative abundance across samples")
      HTML(paste(m1, '<br/>', sep = '<br/>'))
    } else{
      return(NULL)
    }
  })
  
  
  ##### new spec select plot
  
  spec_select_df <- reactive({
    event.data <- event_data("plotly_click", source = "g_select")
    
    # If NULL dont do anything
    validate(
      need(is.null(event.data) == F, "") %then%
        need(event.data$y %in% input$acc_list, "")
    )
   
    select_gfam <- event.data$y
    
    full_full <- full_file_feature()
    col_keep <- colnames(full_full[,4:ncol(full_full)])
    full_select <- data.frame(subset(full_full, Feature %in% select_gfam))
    colnames(full_select) <- c("Acc", "Feature", "Taxa", col_keep)

    validate(
      need(input$input_type == "Biobakery", "Specific Taxa Assignments not availble for EBI inputs")
    )
    
    go_nums <- full_select[,4:ncol(full_select)]
    go_reorder <- go_nums[,unlist(grouped_samps())]
    
    RA <- unlist(go_reorder)
    
    samp_order <- factor(colnames(go_reorder), 
                         levels = c(colnames(go_reorder)))
    Sample_id <- rep(samp_order, each = nrow(go_reorder))
    
    Taxa <- rep(full_select$Taxa, ncol(go_reorder))
    
    
    ## groups
    groupings = new_group_names()
    grouping_nums = group_dims()
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    Group = rep(bin, each = nrow(go_reorder))
    
    plotter <- data.frame(Sample_id, Taxa, Group, RA)
    
  })
  
  spec_select_ggplot <- reactive({
    
    plotter <- spec_select_df()
    
    event.data <- event_data("plotly_click", source = "g_select")
    select_gfam <- event.data$y
    
    #### filter to speed up plotting!!!!
    plotter2 <- aggregate(plotter$RA, list(plotter$Taxa), sum)
    
    summer <- sum(plotter2[,2])
    plotter2$prop <- plotter2[,2]/summer
    
    spec_up <- subset(plotter2, prop > 0.01)
    spec_down <- subset(plotter2, prop < 0.01)
    
    high_tax <- as.character(unlist(spec_up[,1]))
    low_tax <- as.character(unlist(spec_down[,1]))
    
    spec_high <- subset(plotter, Taxa %in% high_tax)
    spec_low <- subset(plotter, Taxa %in% low_tax)
    
    spec_low$Taxa <- factor(rep("Less than 1%", nrow(spec_low)))
    
    plotter2 <- rbind(spec_high, spec_low)
    
    ####
    plotter <- plotter2
    
    bug_sorter <- plotter[order(plotter$RA, decreasing = T),]
    bug_sort <- as.character(unique(bug_sorter$Taxa))
    plotter$Taxa <- factor(plotter$Taxa, levels = bug_sort)
    
    tax_sel <- ggplot(plotter, aes(x=Sample_id, y = RA, fill = Taxa, text = Group)) + 
      stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
      scale_fill_manual(values = randomColor(length(unique(plotter$Taxa))),
                        na.value = "grey") +
      ggtitle(paste0("Taxa Contribution: ", select_gfam)) +
      xlab("Sample") + ylab("Relative Abundance") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank())
    tax_sel
  })
  
  
  output$spec_select <- renderPlotly({
    tax_sel <- spec_select_ggplot() + theme(legend.position='none',
                                            axis.text.y = element_text(size = 15),
                                            axis.title = element_text(size = 15))
    
    taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
    taxly <- taxly %>% layout(margin=list(l = 100), yaxis=list(tickprefix=" "))
    
    
    p3 <- subplot(colorer(), taxly,
                  nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = T, titleY = T)
    p3
  })
  
  output$curr_select_search <- renderText({
    event.data <- event_data("plotly_click", source = "g_select")
    select_gfam <- event.data$y
    message <- paste0("Current Selection from Plot: ", select_gfam)
    message
  })
  
  output$curr_select_exp <- renderText({
    event.data <- event_data("plotly_click", source = "g_exp")
    select_gfam <- event.data$y
    message <- paste0("Current Selection from Plot: ", select_gfam)
    message
  })
  
  output$total_select_exp <- renderText({
    thing <- da_feat()
    paste0("Number of features plotted: ", length(thing))
  })
  
  output$curr_taxa_exp <- renderText({
    event.data <- event_data("plotly_click", source = "t_exp")
    select_taxa <- event.data$y
    message <- paste0("Current Selection from Plot: ", select_taxa)
    message
  })
  
  output$total_taxa_exp <- renderText({
    thing <- da_taxa()
    paste0("Number of features plotted: ", length(thing))
  })
  
  select_stat_plot <- reactive({
    event.data <- event_data("plotly_click", source = "g_select")
    validate(
      need(!is.null(event.data), "") #%then%
        #need(##calc diff exp??)##
    )

    statter2 <- da_feat_stat()
    statter2 <- statter2[-1,]
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    #### now select feat
    
    select_gfam <- event.data$y
    
    statter2 <- subset(statter2, Feature %in% select_gfam)
    
    total_ints <- c()
    total_int <- combn(groupings, 2)
    for (i in 1:ncol(total_int)){
      int <- total_int[,i]
      ints <- paste0(int[1], "-", int[2])
      total_ints <- c(total_ints, ints)
    }
    
    have <- statter2$Int
    add <- setdiff(total_ints, have)
    have_ints <- c(have, add)
    
    num_fill <- round(statter2$we.eBH, 4)
    if (nrow(statter2) != ncol(combn(groupings, 2))){
      how_short <- ncol(combn(groupings, 2)) - nrow(statter2)
      add <- rep(NA, how_short)
      num_fill <- c(num_fill, add)
    }
    
    name_df <- data.frame(name=have_ints, num = num_fill)
    rownames(name_df) <- name_df$name
    name_df = name_df[total_ints,]
    
    group_num <- length(groupings)
    group1 <- gsub('-.*', "", name_df$name[1])
    group_rest <- gsub(".*-", "", name_df$name[1:(group_num-1)])
    group_names <- c(group1, group_rest)
    
    num_fill2 <- name_df$num
    
    resm <- matrix(NA, group_num, group_num)
    resm[lower.tri(resm) ] <- num_fill2
    resm <- t(resm)
    resm[lower.tri(resm) ] <- num_fill2
    rownames(resm) <- group_names
    colnames(resm) <- group_names
    print(resm)
    resm
    
    resm_lab <- resm
    resm_lab[resm_lab < 0.00001] <- "<0.00001"
    
    breaker <- seq(0, 1, by = 0.005)
    coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
               c(colorRampPalette(c("white"))(n=189)))
    
    if (input$numInputs == 2){
      resm[1,2] = resm[1,2] + 0.0000001
    }
    if (sum(resm, na.rm = T) == 0){
      resm[1,2] = 1e-7
    }
    
    par(cex.main=0.8)
    
    v = heatmap.2(data.matrix(resm),
              cellnote = resm_lab,
              density.info="none",
              trace=c("none"),
              notecol="black",
              key.title = NULL,
              breaks = breaker,
              col = coler,
              dendrogram = 'none',
              Rowv=F,
              Colv=F,
              margins=c(10,10),
              cexRow=1.2,
              cexCol=1.2
    )
    v
    title(paste0("Significance between groups for \n", select_gfam), line= -2.5)
  })
  
  output$select_stat_heat <- renderPlot({
    select_mat <- select_stat_plot()
    select_mat
  })
  
  
  output$sel_map_description <- renderUI({
    validate(
      need(length(input$acc_list)>0, "")
    )
    event.data <- event_data("plotly_click", source = "g_select")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
      m1 <- c("Click the Heatmap for Differential Abundance Results")
      m2 <- c("Additional feature specific Taxa Information will appear if available")
      HTML(paste('<br/>', '<br/>', m1, m2, '<br/>', sep = '<br/>'))
    } else {
      return(NULL)
    }
  })
  
  
  observe({
    validate(
      need(input$numInputs > 1, "Need at least 2 groups for ALDEX Test")
    )
    output$select_feat_stat_ui <- renderUI({
      fluidPage(
        uiOutput("sel_map_description"),
        tags$head(tags$style(
          "#sel_map_description {color:#df691a; font-size:18px}")),
        plotOutput("select_stat_heat"),
        fluidRow(fluidPage(downloadButton("sel_stat_download", "Download Plot"),
                 downloadButton("sel_stat_results", "Download Stats")))
      )
    })
  })
  
  
  #####
  ### download uis
  
  output$sel_explore_download <- downloadHandler(
    filename = function() { paste("select_explore", '.png', sep='') },
    content = function(file) {
      p3 <- plot1_plotly()
      
      p3$width = 1200
      p3$height = 800
      
      export(p3, file = file)
      
    })
  
  
  
  output$sel_explore_taxa_download <- downloadHandler(
    filename = function() { paste("select_taxa", '.png', sep='') },
    content = function(file) {
      
      tax_sel <- spec_select_ggplot()+theme(legend.position='none',
                                            axis.text.y = element_text(size = 18),
                                            axis.title = element_text(size = 18))
      
      taxly <- ggplotly(tax_sel, tooltip = c('fill', 'text'))
      taxly <- taxly %>% layout(margin=list(l = 100), yaxis=list(tickprefix=" "))
      
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T, titleY = T)
      p3
      
      p3$width = 1200
      p3$height = 800
      
      export(p3, file = file)
      
    }
  )
  
  sel_explore_taxa_legend <- reactive({
    tax_sel <- spec_select_ggplot()
    legend1 <- g_legend(tax_sel+guides(fill=guide_legend(ncol=3)))
    legend1
  })
  
  output$sel_explore_taxa_legend_download <- downloadHandler(
    filename = function() { paste("select_taxa_legend", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_select_df()$Taxa))
      png(file, width = 35, height = 0.3*(species), units ='cm', res = 300)
      grid.draw(sel_explore_taxa_legend())
      dev.off()
    }
  )
  
  
  output$sel_stat_download <- downloadHandler(
    filename = function() { paste("select_differential_features_heat", '.png', sep='') },
    content = function(file) {
      event.data <- event_data("plotly_click", source = "g_select")
      validate(
        need(!is.null(event.data), "") #%then%
        #need(##calc diff exp??)##
      )
      
      statter2 <- da_feat_stat()
      statter2 <- statter2[-1,]
      
      groupings = new_group_names()
      grouping_nums = group_dims()
      
      #### now select feat
      
      select_gfam <- event.data$y
      
      statter2 <- subset(statter2, Feature %in% select_gfam)
      
      total_ints <- c()
      total_int <- combn(groupings, 2)
      for (i in 1:ncol(total_int)){
        int <- total_int[,i]
        ints <- paste0(int[1], "-", int[2])
        total_ints <- c(total_ints, ints)
      }
      
      have <- statter2$Int
      add <- setdiff(total_ints, have)
      have_ints <- c(have, add)
      
      num_fill <- round(statter2$we.eBH, 4)
      if (nrow(statter2) != ncol(combn(groupings, 2))){
        how_short <- ncol(combn(groupings, 2)) - nrow(statter2)
        add <- rep(NA, how_short)
        num_fill <- c(num_fill, add)
      }
      
      name_df <- data.frame(name=have_ints, num = num_fill)
      rownames(name_df) <- name_df$name
      name_df = name_df[total_ints,]
      
      group_num <- length(groupings)
      group1 <- gsub('-.*', "", name_df$name[1])
      group_rest <- gsub(".*-", "", name_df$name[1:(group_num-1)])
      group_names <- c(group1, group_rest)
      
      num_fill2 <- name_df$num
      
      resm <- matrix(NA, group_num, group_num)
      resm[lower.tri(resm) ] <- num_fill2
      resm <- t(resm)
      resm[lower.tri(resm) ] <- num_fill2
      rownames(resm) <- group_names
      colnames(resm) <- group_names
      print(resm)
      resm
      
      resm_lab <- resm
      resm_lab[resm_lab < 0.00001] <- "<0.00001"
      
      breaker <- seq(0, 1, by = 0.005)
      coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
                 c(colorRampPalette(c("white"))(n=189)))
      
      if (input$numInputs == 2){
        resm[1,2] = resm[1,2] + 0.0000001
      }
      if (sum(resm, na.rm = T) == 0){
        resm[1,2] = 1e-7
      }
      
      dim_fact <- dim(resm)[1]
      if (dim_fact < 6) {
        fill_size <- 2-(0.17*dim_fact)
      } else {
        fill_size <- 2-(0.19*dim_fact)
      }
      
      png(file, width = 15, height = 12, units = 'cm', res = 300)
      par(cex.main=0.8)
      v = heatmap.2(data.matrix(resm),
                    cellnote = resm_lab,
                    notecex = fill_size,
                    density.info="none",
                    trace=c("none"),
                    notecol="black",
                    key.title = NULL,
                    breaks = breaker,
                    col = coler,
                    dendrogram = 'none',
                    Rowv=F,
                    Colv=F,
                    margins=c(7,7),
                    cexRow=1.2,
                    cexCol=1.2
      )
      v
      title(paste0("Significance between groups for \n", select_gfam), line= -1.5)
      
      dev.off()
    }
  )
  
  output$sel_stat_results <- downloadHandler(
    filename = function() { paste("select_differential_features_aldex", '.txt', sep='') },
    content = function(file) {
      event.data <- event_data("plotly_click", source = "g_select")
      validate(
        need(!is.null(event.data), "") #%then%
        #need(##calc diff exp??)##
      )
      
      statter2 <- da_feat_stat()
      statter2 <- statter2[-1,]
      
      groupings = new_group_names()
      grouping_nums = group_dims()
      
      #### now select feat
      
      select_gfam <- event.data$y
      
      statter2 <- subset(statter2, Feature %in% select_gfam)
      
      total_ints <- c()
      total_int <- combn(groupings, 2)
      for (i in 1:ncol(total_int)){
        int <- total_int[,i]
        ints <- paste0(int[1], "-", int[2])
        total_ints <- c(total_ints, ints)
      }
      
      have <- statter2$Int
      add <- setdiff(total_ints, have)
      have_ints <- c(have, add)
      
      num_fill <- round(statter2$we.eBH, 4)
      if (nrow(statter2) != ncol(combn(groupings, 2))){
        how_short <- ncol(combn(groupings, 2)) - nrow(statter2)
        add <- rep(NA, how_short)
        num_fill <- c(num_fill, add)
      }
      
      name_df <- data.frame(name=have_ints, num = num_fill)
      rownames(name_df) <- name_df$name
      name_df = name_df[total_ints,]
      colnames(name_df) <- c("Interaction", "Adjusted_p_value")
      
      write.table(name_df, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  output$expression_download <- downloadHandler(
    filename = function() { paste("gene_family_abundance", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = expression_plot(), device = 'png',
             width = 35, height = 20, units = "cm")
    }
  )
   
  ## Plots for Taxa! ##
  
  ## Taxa instructions
  
  output$ex_delimiter <- renderUI({
    str1 = c("Select Taxa Level of Interest")
    str3 = c("1=Kingdom, 2=Phylum, 3=Class...7=Species")
    HTML(paste(str1, str3, '<br/>', sep = '<br/>'))
  })
  
  
  ## UI based on group nums
  
  output$da_taxa_ui <- renderUI({
    if (input$numInputs > 1){
      fluidPage(
        fluidRow(uiOutput("taxa_selectors")),
        fluidRow(textOutput("curr_taxa_exp")),
        tags$head(tags$style("#curr_taxa_exp{font-size: 20px}")),
        fluidRow(textOutput("total_taxa_exp")),
        tags$head(tags$style("#total_taxa_exp{font-size: 20px}")),
        fluidRow(column(12, uiOutput("key2"))),
        fluidRow(column(12, uiOutput("da_taxa_heat_UI"))),
        fluidRow(downloadButton("species_heat_download", "Download Heatmap")),
        fluidRow(column(8, uiOutput("da_taxa_stat_ui"))),
        fluidRow(downloadButton("species_stat_download", "Download Heatmap"),
                 downloadButton("species_stat_results", "Download Stats"))
      )
    }
  })
  
  
  # Taxa explore #

  # matrix reordered by group and collapsed by species
  reorder_spec_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = spec_nums()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    rownames(exp2) <- spec_full()$Taxa
    return(exp2)
  })

  # Explore Taxa dataframe
  spec_taxa_data <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }
    
    withProgress({

      spec_plot_df = reorder_spec_mat()
      specs = spec_plot_df
      
      specs1 <- cbind(rownames(specs), specs)
      colnames(specs1)[1] <- "Taxa"
      
      if (input$input_type == "EBI"){
        specs1$Taxa <- gsub(";k__", "_k__", specs1$Taxa)
        
        tax_sep <- strsplit(as.character(specs1$Taxa), ";")
        ebi_new <- data.frame(Kingdom = NA, Phylum = NA, 
                              Class = NA, Order = NA, Family = NA,
                              Genus = NA, Species = NA)
        
        for (i in 1:length(tax_sep)){
          curr <- tax_sep[[i]]
          k = curr[1]
          p = gsub("p__", "", curr[2])
          c = gsub("c__", "", curr[3])
          o = gsub("o__", "", curr[4])
          f = gsub("f__", "", curr[5])
          g = gsub("g__", "", curr[6])
          s = gsub("s__", "", curr[7])
          
          taxs <- data.frame(Kingdom = k, Phylum = p, Class = c,
                             Order = o, Family = f, Genus = g, Species = s)
          ebi_new <- rbind(ebi_new, taxs)
        }
        
        ebi_new <- ebi_new[-1,]
        ebi_new[is.na(ebi_new)] <- ""
        
        new_spec <- cbind(ebi_new[,input$taxaDims], specs1[,-1])
        colnames(new_spec)[1] = "Taxa"
        
        
        ### blanks?
        new_new_spec <- c()
        for (i in 1:length(new_spec$Taxa)){
          curr <- as.character(new_spec$Taxa[i])
          if (curr == ""){
            curr = "undetermined"
          }
          new_new_spec <- c(new_new_spec, curr)
        }
        ##
        specs1$Taxa <- new_new_spec
        specs <- specs1
        
      }
      if (input$input_type == "Biobakery"){
        if (input$taxaDims == 6){
          specs1$Taxa <- gsub("\\..*", "", specs1$Taxa)
        }
        specs <- specs1
      }
      
      
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
      
      groupings = new_group_names()
      grouping_nums = group_dims()
      
      group_titles = c()
      for (i in 1:length(groupings)) {
        title = groupings[i]
        reps = grouping_nums[i]
        title_rep = rep(title, reps)
        group_titles = c(group_titles, title_rep)
      }
      
      group_titles2 = rep(group_titles, spec_row)
      
      samp_order <- factor(colnames(spec_RelExp), 
                           levels = c(colnames(spec_RelExp)))
      sample_id <- rep(samp_order, nrow(spec_RelExp))
      
      spec_spec <- rep(specs$Taxa, each = spec_col)
      
      spec_g_data = data.frame(spec_datas,
                               sample_id, group_titles2, spec_spec)
      colnames(spec_g_data) = c("Relative_Abundance", "Sample_num", "Group", "Taxa")
      
      spec_g_data
    }, message = "Collecting Taxanomic Data")
  })

  # Explore Taxa Species Plot Object
  species_explore_plot <- reactive({
    withProgress({
      spec_g_data <- spec_taxa_data()
      uniq_spec_num = length(unique(spec_g_data$Taxa))
      
      spec_g_data$Taxa <- reorder(spec_g_data$Taxa, -spec_g_data$Relative_Abundance)
      
      spec_g_data2 <- aggregate(spec_g_data$Relative_Abundance, list(spec_g_data$Taxa), sum)
      
      summer <- sum(spec_g_data2[,2])
      spec_g_data2$prop <- spec_g_data2[,2]/summer
  
      spec_up <- subset(spec_g_data2, prop > 0.01)
      spec_down <- subset(spec_g_data2, prop < 0.01)
      
      high_tax <- as.character(unlist(spec_up[,1]))
      low_tax <- as.character(unlist(spec_down[,1]))
      
      spec_high <- subset(spec_g_data, Taxa %in% high_tax)
      spec_low <- subset(spec_g_data, Taxa %in% low_tax)
      
      spec_low$Taxa <- factor(rep("Less than 1%", nrow(spec_low)))
      
      spec_g_data_filt <- rbind(spec_high, spec_low)
      
      exp_plot = ggplot(spec_g_data_filt, aes(x = Sample_num,
                                              y = Relative_Abundance,
                                              fill = Taxa, text = Group)) +
        stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
        ggtitle("Total Taxonomic Contribution") +
        scale_fill_manual(values = randomColor(uniq_spec_num)) +
        ylab("Relative Abundance") + xlab("Sample")
    }, message = "Organizing Taxa")
    exp_plot
  })

  output$species_explore_placehold <- renderUI({
    m1 <- c("Taxonomic Distribution across all samples")
    HTML(paste(m1, '<br/>', sep = '<br/>'))
  })
  
  output$species_explore <- renderPlotly({
    withProgress({
      taxa_all <- species_explore_plot() + 
        theme(legend.position = 'none',
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.title = element_text(size = 15))#legend.text = element_text(size = 8))
      
      taxly <- ggplotly(taxa_all, tooltip = c("fill", "text"))
      
      taxly <- taxly %>% layout(margin=list(l = 100), yaxis=list(tickprefix=" "))
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T, titleY = T)
      p3
    }, message = "Rendering")
  })

  species_explore_legend <- reactive({
    legend1 <- g_legend(species_explore_plot()+guides(fill=guide_legend(ncol=3)))
    legend1
  })

  output$species_download <- downloadHandler(
    filename = function() { paste("taxa_explore", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_taxa_data()$Taxa))
      taxa_all <- species_explore_plot() + 
        theme(legend.position = 'none',
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 18),
              axis.title = element_text(size = 18))#legend.text = element_text(size = 8))
      
      taxly <- ggplotly(taxa_all, tooltip = c("fill", "text"))
      taxly <- taxly %>% layout(margin=list(l = 100), yaxis=list(tickprefix=" "))
      
      p3 <- subplot(colorer(), taxly,
                    nrows = 2, margin = c(0,0,-0.01,0),
                    heights = c(0.1, 0.9), shareX = T, titleY = T)
      
      p3$width = 1200
      p3$height = 800
      
      export(p3, file = file)
      
    }
  )

  output$species_legend_download <- downloadHandler(
    filename = function() { paste("taxa_explore_legend", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_taxa_data()$Taxa))
      png(file, width = 35, height = 0.3*(species/3), units ='cm', res = 300)
      grid.draw(species_explore_legend())
      dev.off()
    }
  )

  output$species_raw_data <- downloadHandler(
    filename = function() { paste("taxa_explore_raw_data", '.txt', sep='') },
    content = function(file) {
      spec_g_data <- spec_taxa_data()
      uniq_spec_num = length(unique(spec_g_data$Taxa))

      spec_g_data$Taxa <- reorder(spec_g_data$Taxa, -spec_g_data$Relative_Abundance)
      write.table(spec_g_data, file, quote=F, sep='\t', row.names = F)
    }
  )

  observe({
    output$taxa_selectors <- renderUI({
      fluidRow(
        column(8, fluidPage(
          fluidRow(
            column(4, sliderInput("taxa_pval_up", "Adjusted p Value", min = 0, max = 0.5,
                                  step = 0.001, value = 0.05)),
            column(2, numericInput("taxa_pval_num", " ", min = 0, max = 0.1, 
                                   step = 0.001, value = 0.05)),
            column(4, sliderInput("taxa_effect", "Effect Size", min = 0, max = 10,
                                  step = 0.25, value = 5)),
            column(2, numericInput("taxa_effect_num", " ", min = 0, max = 10,
                                  step = 0.25, value = 5))
          ),
          fluidRow(column(12, plotOutput("da_pval_distr", height = 150)))
        )),
        column(4, fluidPage(fluidRow(plotOutput("da_pval_scatter", height = 260))))
      )
    })
    
    observe({
      validate(need(input$taxa_pval_num, "wait"))
      
      updateSliderInput(session, "taxa_pval_up", value = input$taxa_pval_num)
      updateSliderInput(session, "taxa_effect", value = input$taxa_effect_num)
      
      
    })
    
    output$feat_selectors <- renderUI({
      fluidRow(
        column(8, fluidPage(
          fluidRow(
            column(4, sliderInput("feat_pval_up", "Adjusted p Value", min = 0, max = 0.1,
                                  step = 0.001, value = 0.05)),
            column(2, numericInput("feat_pval_num", " ", min = 0, max = 0.1, 
                                   step = 0.001, value = 0.05)),
            column(4, sliderInput("feat_effect", "Effect Size", min = 0, max = 10,
                                  step = 0.25, value = 5)),
            column(2, numericInput("feat_effect_num", " ", min = 0, max = 10,
                                   step = 0.25, value = 5))
          ),
          fluidRow(column(12, plotOutput("feat_pval_distr", height = 150)))
        )),
        column(4, fluidPage(fluidRow(plotOutput("feat_pval_scatter", height = 260))))
      )
    })
    
    observe({
       validate(need(input$feat_pval_num, "wait"))
      
       updateSliderInput(session, "feat_pval_up", value = input$feat_pval_num)
       updateSliderInput(session, "feat_effect", value = input$feat_effect_num)
     })
  
  })
    

 ##### Test Differential Abundance for Taxa
  
  da_taxa_stat <- reactive({
    if (input$input_type == "Biobakery"){
      if (input$testme){
        validate(
          need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
        )
      }
      else {validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
          need(unlist(grouped_samps()), "Please select samples in the Group Tab")
      )
      }
    } else {
      validate(
        need(input$taxa1, "Please provide a Taxa File in the Upload Tab")%then%
          need(unlist(grouped_samps()), "Please select samples in the Group Tab")
      )
    }
    spec_order <- reorder_spec_mat()
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    group_titles = c()
    for (i in 1:length(groupings)) {
      title = groupings[i]
      reps = grouping_nums[i]
      title_rep = rep(title, reps)
      group_titles = c(group_titles, title_rep)
    }
    
    conds <- group_titles
    
    combos <- combn(unique(conds), 2)
    stat_results = data.frame(we.ep = NA, we.eBH = NA, wi.ep = NA, wi.eBH = NA,
                              rab.all = NA, rab.win.a = NA, rab.win.b = NA,
                              diff.btw = NA, diff.win = NA, effect = NA, overlap = NA,
                              Int = NA, Feature = NA)
    for (i in 1:ncol(combos)){
      withProgress({
        ab <- combos[,i]
        a <- ab[1]
        b <- ab[2]
        a_sel = which(conds == a)
        b_sel = which(conds == b)
        
        feat_subset <- cbind(spec_order[,a_sel], spec_order[,b_sel])
        cond_subset <- c(conds[a_sel], conds[b_sel])
        x <- aldex.clr(feat_subset, cond_subset, mc.samples=16, verbose=TRUE)
        x.tt <- try(aldex.ttest2(x, cond_subset, paired.test = F))
        x.effect <- try(aldex.effect(x, cond_subset, include.sample.summary=FALSE, verbose=TRUE))
        
        validate(
          need(class(x.tt)!= "try-error", 
               paste0("Differential Abundance failed due to the following error: ", x.tt))%then%
            need(class(x.effect)!= "try-error", 
                 paste0("Differential Abundance failed due to the following error: ", x.effect))
        )
        
        x.all <- data.frame(x.tt, x.effect, stringsAsFactors=FALSE)
        x.all$Int <- rep(paste0(a,"-", b), nrow(x.all))
        x.all$Feature <- rownames(x.all)
        names(x.all) = names(stat_results)
        stat_results <- rbind(stat_results, x.all)
      }, message = paste0("Calculating Differential Abundance ", i, " of ", ncol(combos))) 
    }
    up_effect <- max(abs(stat_results$effect), na.rm=T)
    up_pval <- max(abs(stat_results$we.eBH), na.rm=T)
    
    updateSliderInput(session, "taxa_pval_up", max = round(up_pval+0.05, 2))
    updateNumericInput(session, "taxa_pval_num", max = round(up_pval+0.05, 2))
    updateSliderInput(session, "taxa_effect", max = ceiling(up_effect))
    updateNumericInput(session, "taxa_effect_num", max = ceiling(up_effect))
    
    stat_results
  })
  
  output$da_pval_distr <- renderPlot({
    stat_results <- da_taxa_stat()
    par(mfrow=c(1,2), mar = c(1,1,1,1))
    cd <- density(stat_results$we.eBH, na.rm = T)
    plot(cd,
         xlim = c(0, 
                  round(max(stat_results$we.eBH, na.rm = T)+0.05, 2)),
         frame.plot = F, ann = F, axes = F, col = 'orange', lwd = 5)
    abline(v=input$taxa_pval_up, col = 'red', lwd = 5)
    
    ce <- density(abs(stat_results$effect), na.rm = T)
    plot(ce, 
         xlim = c(0,
                  ceiling(max(abs(stat_results$effect), na.rm=T))),
         frame.plot = F, ann = F, axes = F, col = 'orange', lwd = 5)
    abline(v=input$taxa_effect, col = 'dodgerblue2', lwd = 5)
  }, bg = 'transparent')
  
  output$da_pval_scatter <- renderPlot({
    stat_results <- da_taxa_stat()
    print(max(stat_results$effect))
    print(min(stat_results$we.ebh))
    plot(stat_results$we.eBH, abs(stat_results$effect), pch = 16,
         xlab = "Adjusted p Value", ylab = "Effect Size",
         xlim = c(0, round(max(stat_results$we.eBH, na.rm = T)+0.05, 2)),
         ylim = c(0, ceiling(max(abs(stat_results$effect), na.rm=T))))
    rect(0,input$taxa_effect, input$taxa_pval_up, ceiling(max(abs(stat_results$effect), na.rm=T)),
         lwd = 1.5, col = NA, border = 'orange')
  })
  
  da_taxa <- reactive({
    stat_results <- da_taxa_stat()
    keepers_df <- subset(stat_results, wi.eBH < input$taxa_pval_up &
                               we.eBH < input$taxa_pval_up)
    
    keepers <- unique(subset(keepers_df, effect > input$taxa_effect |
                               effect < -input$taxa_effect)$Feature)
    
    print(keepers_df)
    print(length(keepers))
    validate(
      need(length(keepers) > 0, "No signficant features were detected!")
    )
    keepers
  })
  
  da_taxa_heat_plotly <- reactive({
    da_taxa <- da_taxa()
    spec_order <- reorder_spec_mat()
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    ### colors
    bin <- c()
    for (i in 1:length(groupings)){
      use <- rep(groupings[i], grouping_nums[i])
      bin <- c(bin, use)
    }
    
    go_show <- spec_order[da_taxa,]

    go_show_nums2 <- log10(go_show + 1)
    
    gg_nums <- go_show_nums2
    
    validate(
      need(nrow(gg_nums) > 1, "At least two features required for plotting, try adjusting the sliders above!")
    )
    feat_clust <- hclust(dist(gg_nums))
    
    gg_nums_feat <- gg_nums[feat_clust$order,]
    gg_nums_samp <- gg_nums_feat
    
    gg_names <- factor(colnames(gg_nums_samp), levels = colnames(gg_nums_samp))
    gg_feat <- factor(rownames(gg_nums_samp), levels = rownames(gg_nums_samp))
    
    x <- sweep(gg_nums_samp, 1L, rowMeans(gg_nums_samp, na.rm = T), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = T)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    
    Z_scored_CPM <- unlist(x)

    
    names <- rep(gg_names, each = nrow(gg_nums))
    feat <- rep(gg_feat, ncol(gg_nums))
    groups <- rep(bin, each = nrow(gg_nums))
    
    
    plotter <- data.frame(Z_scored_CPM, names, groups, feat)
    ###
    ##
    ### make gmain in plotly!!!!
    library(plotly)
    
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      tickcolor = 'white',
      showgrid = FALSE
    )
    p <- plot_ly(source = "t_exp",
                 x = plotter$names, y = plotter$feat,
                 text = plotter$groups,
                 z = plotter$Z_scored_CPM, type = "heatmap",
                 hoverinfo = 'y+text',
                 colorbar = list(x = -0.3, 
                                 xanchor = 'left',
                                 tickmode='array',
                                 tickvals = c(min(plotter$Z_scored_CPM),max(plotter$Z_scored_CPM)),
                                 ticktext = c("low", "high"), tickfont = list(size=18))) %>%
      layout(yaxis = ax, xaxis = ax, margin=list(l = 100))
    
    
    ####
    p3=subplot(colorer(),
                  p, nrows = 2, margin = c(0,0,-0.01,0),
                  heights = c(0.1, 0.9), shareX = TRUE)
    
    p3
    
  })
  
  output$da_taxa_heat <- renderPlotly({
    p3 <- da_taxa_heat_plotly()
    p3
  })
  
  output$da_taxa_heat_UI <- renderUI({
    plotlyOutput("da_taxa_heat", height = 450)
  })
  
  da_taxa_stat_mat <- reactive({
    event.data <- event_data("plotly_click", source = "t_exp")
    
    validate(
      need(is.null(event.data) == F, "")
    )
    
    validate(
      need(event.data$y %in% da_taxa_stat()$Feature, "")
    )
    
    groupings = new_group_names()
    grouping_nums = group_dims()
    
    #### now select bug
    
    select_taxa <- event.data$y
    
    statter <- da_taxa_stat()
    statter2 <- subset(statter, Feature %in% select_taxa)
    
    total_ints <- c()
    total_int <- combn(groupings, 2)
    for (i in 1:ncol(total_int)){
      int <- total_int[,i]
      ints <- paste0(int[1], "-", int[2])
      total_ints <- c(total_ints, ints)
    }
    
    have <- statter2$Int
    add <- setdiff(total_ints, have)
    have_ints <- c(have, add)

    num_fill <- round(statter2$we.eBH, 4)
    if (nrow(statter2) != ncol(combn(groupings, 2))){
      how_short <- ncol(combn(groupings, 2)) - nrow(statter2)
      add <- rep(NA, how_short)
      num_fill <- c(num_fill, add)
    }
    
    name_df <- data.frame(name=have_ints, num = num_fill)
    rownames(name_df) <- name_df$name
    name_df = name_df[total_ints,]
    
    group_num <- length(groupings)
    group1 <- gsub('-.*', "", name_df$name[1])
    group_rest <- gsub(".*-", "", name_df$name[1:(group_num-1)])
    group_names <- c(group1, group_rest)
    
    num_fill2 <- name_df$num
    
    resm <- matrix(NA, group_num, group_num)
    resm[lower.tri(resm) ] <- num_fill2
    resm <- t(resm)
    resm[lower.tri(resm) ] <- num_fill2
    rownames(resm) <- group_names
    colnames(resm) <- group_names
    print(resm)
    resm
  })
  
  
  
  da_taxa_resm <- reactive({
    event.data <- event_data("plotly_click", source = "t_exp")
    resm <- da_taxa_stat_mat()
    
    resm_lab <- resm
    resm_lab[resm_lab < 0.00001] <- "<0.00001"
    
    breaker <- seq(0, 1, by = 0.005)
    coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
               c(colorRampPalette(c("white"))(n=189)))
    
    select_taxa <- event.data$y
    new_name <- unlist(strsplit(select_taxa, ";"))
    new_name <- new_name[length(new_name)]
    
    if (input$numInputs == 2){
      resm[1,2] = resm[1,2] + 0.0000001
    }
    if (sum(resm, na.rm = T) == 0){
      resm[1,2] = 1e-7
    }
    
    par(cex.main=0.8)
    v = heatmap.2(data.matrix(resm),
              cellnote = resm_lab,
              density.info="none",
              trace=c("none"),
              notecol="black",
              key.title = NULL,
              breaks = breaker,
              col = coler,
              dendrogram = 'none',
              Rowv=F,
              Colv=F,
              margins=c(10,10),
              cexRow=1.2,
              cexCol=1.2
    )
    v
    title(paste0("Significance between groups for \n", new_name), line= -2.5)
  })

  
    
  output$da_taxa_stat_heat <- renderPlot({
    event.data <- event_data("plotly_click", source = "t_exp")
    validate(
      need(is.null(event.data) == F, "")
    )
    plott <- da_taxa_resm()
    plott
  })

  output$taxa_map_description <- renderUI({
    event.data <- event_data("plotly_click", source = "t_exp")
    
    # If NULL dont do anything
    if (is.null(event.data) == T){
      m1 <- c("Click the Heatmap for Differential Abundance Results")
      HTML(paste('<br/>', m1, '<br/>', sep = '<br/>'))
    } else {
      return(NULL)
    }
  })
  
  output$da_taxa_stat_ui <- renderUI({
    fluidRow(fluidPage(
      uiOutput("taxa_map_description"),
      tags$head(tags$style(
        "#taxa_map_description {color:#df691a; font-size:18px}")),
      plotOutput("da_taxa_stat_heat")
    ))
  })
  
  
  output$species_heat_download <- downloadHandler(
    filename = function() { paste("differential_taxa", '.png', sep='') },
    content = function(file) {
      species = length(unique(spec_taxa_data()$Taxa))
      
      p3 <- da_taxa_heat_plotly()

      p3$width = 1200
      p3$height = 800
      
      export(p3, file = file)
    }
  )
  
  output$species_stat_download <- downloadHandler(
    filename = function() { paste("differential_taxa_heat", '.png', sep='') },
    content = function(file) {
      event.data <- event_data("plotly_click", source = "t_exp")
      resm <- da_taxa_stat_mat()
      
      resm_lab <- resm
      resm_lab[resm_lab < 0.00001] <- "<0.00001"
      
      breaker <- seq(0, 1, by = 0.005)
      coler <- c(c(colorRampPalette(c("red", "white"))(n=11)),
                 c(colorRampPalette(c("white"))(n=189)))
      
      select_taxa <- event.data$y
      new_name <- unlist(strsplit(select_taxa, ";"))
      new_name <- new_name[length(new_name)]
      
      if (input$numInputs == 2){
        resm[1,2] = resm[1,2] + 0.0000001
      }
      if (sum(resm, na.rm = T) == 0){
        resm[1,2] = 1e-7
      }
      
      dim_fact <- dim(resm)[1]
      if (dim_fact < 6) {
        fill_size <- 2-(0.17*dim_fact)
      } else {
        fill_size <- 2-(0.19*dim_fact)
      }
      
      png(file, width = 15, height = 12, units = 'cm', res = 300)
      par(cex.main=0.8)
      v = heatmap.2(data.matrix(resm),
                    cellnote = resm_lab,
                    notecex = fill_size,
                    density.info="none",
                    trace=c("none"),
                    notecol="black",
                    key.title = NULL,
                    breaks = breaker,
                    col = coler,
                    dendrogram = 'none',
                    Rowv=F,
                    Colv=F,
                    margins=c(7,7),
                    cexRow=1.2,
                    cexCol=1.2
      )
      v
      title(paste0("Significance between groups for \n", new_name), line= -1.5)
      dev.off()
    }
  )
  
  output$species_stat_results <- downloadHandler(
    filename = function() { paste("differential_taxa_aldex", '.txt', sep='') },
    content = function(file) {
      stat_results <- da_taxa_stat()
      keepers <- subset(stat_results, we.eBH < input$taxa_pval_up &
                          wi.eBH < input$taxa_pval_up)
      write.table(keepers, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )

  # 
  # 
  #### Correlation Plot #####

  ## Correlation selections
  sig_tab <- reactive({
    gfam_DF = acc_full()
    samp_paths = gfam_DF
    samp_paths$Feature
  })

  
  # All sample Correlation Plot #

  
  actual_corr_plot <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    }
    else{
      if (input$input_type == "Biobakery"){
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, '')
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

    heyo <- subset(samp_paths, header_DF %in% corr_list)
    colnames(heyo)[1] = "Feature"
    
    col = (ncol(samp_paths))-1
    row = nrow(heyo)

    heyo_small <- heyo
    rownames(heyo_small) = heyo$Feature
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
                   trace=c("none"),
                   breaks = seq(-1, 1, by = 0.02),
                   col=my_palette,
                   dendrogram="both",
                   margins=c(5,5),
                   cexRow=1,
                   cexCol=1.2,
                   na.color="gray60"
    )
    v
  })

  
  output$corr_placehold <- renderUI({
    if (length(input$sig_select) < 2){
      m1 <- c("Select at least two Features for cross-correlation analysis")
      HTML(paste(m1, '<br/>', sep = '<br/>'))
    } else{
      return(NULL)
    }
  })
  
  output$corr_info <- renderUI({
    if (length(input$sig_select) > 1){
      m1 <- c("Feature number key can be found below generated plots")
      HTML(paste(m1, '<br/>', sep = '<br/>'))
    } else{
      return(NULL)
    }
  })
  
  output$corr_plot <- renderPlot({
    actual_corr_plot()
  })

  actual_corr_names <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    }
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, '')
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
    
    heyo <- subset(samp_paths, header_DF %in% corr_list)
    colnames(heyo)[1] = "Feature"

    actual_corr_names <- heyo$Feature
    actual_corr_names
  })


  # Generate list of correlation matrices for each Group
  group_corr_plist <- reactive({
    if (input$testme) {}
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, '')
    )
    corr_list2 = input$sig_select

    groupings = new_group_names()
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

      heyo <- subset(samp_paths, header_DF %in% corr_list2)
      colnames(heyo)[1] = "Feature"

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small<-heyo
      rownames(heyo_small) = heyo$Feature
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

  # Generate list of correlation significance symbols for each group
  group_sym_plist <- reactive({
    if (input$testme) {}
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, '')
    )
    corr_list2 = input$sig_select

    groupings = new_group_names()
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
      
      heyo <- subset(samp_paths, header_DF %in% corr_list2)
      colnames(heyo)[1] = "Feature"

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small <- heyo
      rownames(heyo_small) = heyo$Feature
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


  output$sig_message <- renderText({
    paste("* = p < 0.05", "** = p < 0.01", "*** = p < 0.001", sep='\n')
  })

  # Generate Table of Selected Gene Families
  corr_label_table <- reactive({
    if (input$testme) {}
    else{
      if (input$input_type == "Biobakery"){
        validate(
          need(input$file1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
      if (input$input_type == "EBI"){
        validate(
          need(input$features1, "Please provide a file in the Upload Tab") %then%
            need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
      }
    }

    validate(
      need(length(input$sig_select) > 1, '')
    )

    Gene_Families = actual_corr_names()
    nums = length(Gene_Families)
    Label = paste(1:nums)
    Label_Matrix = data.frame(Gene_Families, Label)
    Label_Matrix
  })

  output$corr_labels <- renderTable({
    corr_label_table()
  })

  output$corr_table_download <- downloadHandler(
    filename = function() { paste("correlation_sample_key", '.txt', sep='') },
    content = function(file) {
      corr_tab = corr_label_table()
      corr_tab2 = data.frame(corr_tab)
      write.table(corr_tab2, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )

  observeEvent({
    input$sig_select
    }, {
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "")
      )
    }
      else{
        if (input$input_type == "Biobakery"){
          validate(
            need(input$file1, "Please provide a file in the Upload Tab") %then%
              need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
        }
        if (input$input_type == "EBI"){
          validate(
            need(input$features1, "Please provide a file in the Upload Tab") %then%
              need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
        }
      }

    if (input$numInputs > 1) {
     output$group_corrs <- renderUI({ get_plot_output_list(max_plots,
                                                            input$numInputs,
                                                            group_corr_plist(),
                                                            group_sym_plist(),
                                                            new_group_names())
      })
    }
  })


  # Download full correlation heatmap
  output$corr_download <- downloadHandler(
    filename = function() { paste("all_gene_family_correlation", '.png', sep='') },
    content = function(file) {
      png(file, width = 20, height = 15, units ='cm', res = 300)
      corr_list = input$sig_select

      gfam_DF1 = acc_full()
      col_num = ncol(gfam_DF1)
      row_num = nrow(gfam_DF1)
      reorder_samps = gfam_DF1[,2:col_num]
      reorder_samps = reorder_samps[,c(unlist(grouped_samps()))]
      header_DF = gfam_DF1[,1]

      gfam_DF = cbind(header_DF, reorder_samps)
      samp_paths = gfam_DF
      heyo <- subset(samp_paths, header_DF %in% corr_list)
      colnames(heyo)[1] = "Feature"

      col = (ncol(samp_paths))-1
      row = nrow(heyo)

      heyo_small<-heyo
      rownames(heyo_small) = heyo$Feature
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
      
      dim_fact <- dim(corr_mat)[1]
      if (dim_fact > 10) {
        fill_size <- 2-(0.8*dim_fact)
      } else if (dim_fact > 20) {
        fill_size <- 2-(0.12*dim_fact)
      } else {
        fill_size <- 2-(0.1*dim_fact)
      }
      
      v <- heatmap.2(corr_mat,
                     cellnote = sym_mat, notecol="black", notecex = fill_size,
                     main="Correlation Values across \nAll Samples", #heatmap title
                     density.info="none",
                     breaks = seq(-1, 1, by = 0.02),
                     col=my_palette,
                     dendrogram="both",
                     margins=c(7,7),
                     cexRow=1,
                     cexCol=1.2,
                     trace=c("none"),
                     na.color="gray60"
      )
      v
      dev.off()
    }
  )

  ## UI Elements for downloading group correlation ##

  output$group_download <- renderUI({
    lapply(1:input$numInputs, function(i) {
      display_name = new_group_names()
      downloadButton(paste0("downloadData", i), paste("Download", display_name[i], sep=" "))
    })
  })

  observe({
    lapply(1:input$numInputs, function(i) {
      output[[paste0("downloadData", i)]] <- downloadHandler(
        filename = function() { paste(new_group_names()[i], "_correlation", '.png', sep='') },
        content = function(file) {
          png(file, width = 20, height = 15, units ='cm', res = 300)
          download_plot_output_list(max_plots,
                                    1,
                                    group_corr_plist()[i],
                                    group_sym_plist()[i],
                                    new_group_names()[i])
          dev.off()
        }
      )
    })
  })
  
}
