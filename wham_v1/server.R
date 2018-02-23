
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
          "*Citation Info*")
  })
  
  
  ## Load in the data or try a sample dataset
  full_file <- reactive({
    if (input$testme) {
      full_file <- fread("sample_input.tsv", header=TRUE, sep=input$sep)
    }
    else {
      inFile <- input$file1
      if (is.null(inFile)) {
        return(NULL)}
      full_file <- try(
        {fread(inFile$datapath, header=TRUE, sep=input$sep)})
      
      correct_cols <- c("Acc", "Gene_Family", "Species")
      validate(
        need(class(full_file)!="try-error", paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab.")) %then%
          need(all(colnames(full_file)[1:3]==correct_cols), paste0("The input file provided is not an appropriate format. Please view the sample input file provided in the Home Tab."))
      )
    }
    nums <- data.matrix(full_file[,4:ncol(full_file)])
    rownames(nums) <- rownames(full_file)
    if (input$filter_level == 0){
      full_file <- full_file
    }
    else{
      keep_rows = rownames((nums[apply(nums==0,1,sum)<=input$filter_level*ncol(nums),]))
      full_file <- full_file[as.numeric(keep_rows),]
    }
    full_file <- subset(full_file, Gene_Family != "NO_NAME")
    
    full_file$Gene_Family <- gsub("[^[:alnum:]']", "_", full_file$Gene_Family)
    
    full_file
  })
  
  # Preview Full File
  output$contents <- renderDataTable({
    validate(
      need(full_file(),"")
    )
    full_file <- full_file()
    if (nrow(full_file)<50){
      full_file_show <- full_file[,1:10]
    }
    else{
      full_file_show <- full_file[1:50, 1:10]
    }
    full_file_show
  })
  
  ### Generate prerequisites for plotting
  
  # Generate dataframe collapsed by Gene Family
  acc_full <- reactive({
    full_file <- full_file()
    col_num <- ncol(full_file)
    DT <- data.table(full_file[,4:col_num])
    DT$Gene_Family = full_file$Gene_Family
    DT2 <- DT[, lapply(.SD,sum), by = "Gene_Family"]
    DT2
  })
  # acc full but numeric values only
  acc_nums <- reactive({
    accs <- acc_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Gene_Family
    accs3 = as.matrix(accs2)
    accs3
  })
  
  # Generate dataframe collapsed by Species
  spec_full <- reactive({
    full_file <- full_file()
    col_num <- ncol(full_file)
    DT <- data.table(full_file[,4:col_num])
    DT$Species = full_file$Species
    DT2 <- DT[, lapply(.SD,sum), by = "Species"]
    DT2
  })
  
  spec_nums <- reactive({
    accs <- spec_full()
    accs2 = accs[,-1]
    rownames(accs2) = accs$Species
    accs3 = as.matrix(accs2)
    rownames(accs3) <- rownames(accs2)
    accs3
  })
  
  # Generate selection list of gene families
  acc_select <- reactive({
    accs <- acc_full()
    accs$Gene_Family = paste(" ", accs$Gene_Family, " ", sep = "") #ensure unique gene families
    accs$Gene_Family
  })
  
  # When file is uploaded update choice of gene families
  observe({
    if (input$testme) {
      updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
      updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
      }
    else {
      if(is.null(input$file1)){}
      else {
        updateSelectizeInput(session,'acc_list', choices = acc_select(), server = TRUE)
        updateSelectizeInput(session,'excluder', choices = acc_select(), server = TRUE)
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
        fluidRow(selectInput("test_groups", label = "",choices = group1, 
                             multiple = TRUE, selected = group1),
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
      textInput(inputName, label = "Group Prefix")
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
    #validate(need(group_names(), ""))
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
      #print(acc_nums()[1:2,1:2])
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
    return(exp2)
  })
  
  # Generate warning if sample is assigned to multiple groups
  #observeEvent(input$input1, {
  output$group_warning <- renderText({
    validate(
      need(length(unlist(grouped_samps()))>0, 'Select groups!')
    )
    group_list = c(unlist(grouped_samps()))
    duplicates = c(duplicated(group_list))
    if ('TRUE' %in% duplicates) {
      message = c("Warning: A sample has been assigned to more than one group!")
    } else {
      message = c(" ")
    }
    message
  })
  #})
  
  output$tutorialGroup <- renderText({
    if (input$testme){
      message = c("Groups have been named and assigned!\n Please continue to the next Tab")
      message
    }
  })
  
  #### Begin Plotting ! ####
  
  ##
  ## Gene Plots
  ##
  
  # Gene Explore Plot
  
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

    groupings = new_group_names()
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
    colnames(g_data) = c("Relative_Abundance", "Gene_Family", "Sample_num", "Group")
    #g_data without zeros
    
    g_data2 = subset(g_data, Relative_Abundance > 0, select=c(Relative_Abundance,Gene_Family,Sample_num, Group))
    
    if (length(input$excluder)>0){
      exclude_list = input$excluder
      g_data2 = g_data2[grep(paste(exclude_list, collapse="|"), g_data2$Gene_Family, invert=TRUE), ]
    }
    
    uniq_gfam_num = length(unique(g_data2$Gene_Family))
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector[10] = "blue2" 
    
    if (uniq_gfam_num > 74) {
      col_vector_edit = sample(col_vector, uniq_gfam_num, TRUE)
    } else {
      col_vector_edit = col_vector
    }
    
    exp_plot = ggplot(g_data2, aes(x = Group, y = Relative_Abundance, fill = Gene_Family)) +
      #geom_bar(position = 'fill', stat = 'identity') +
      stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
      scale_fill_manual(values = col_vector_edit) +
      ylab("Relative Abundance") +
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
    instruction = c("Hover over bars to view Gene Family Name")
    instruction
  })
  
  # Gene Search Plot
  # Gene search dataframe for plot
  exp_plot <- reactive({
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    }
    else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
          need(unlist(reorder_mat()), "Please select samples in the Group Tab")
      )
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

    groupings = new_group_names()
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
    sample_num = paste(1:(col))
    samp=strtoi(sample_num)
    isamp = rep(samp,row)
    g_data = data.frame(datas, Gfam, isamp, group_titles2)
    colnames(g_data) = c("Exp", "Gfam", "Sample_num", "groups")
    #g_data without zeros
    
  
    g_data2 = subset(g_data, select=c(Exp,Gfam,Sample_num, groups))
    g_data2
  })
  
  # Gene search plot object
  expression_plot <- reactive({
    g_data2 = data.frame(exp_plot())
    class(g_data2$Exp) = "numeric"
    
    g_data2$Gfam2 <- gsub(" ", "", g_data2$Gfam)
    g_data2$Gfam2 <- str_trunc(g_data2$Gfam2, 30, "center")
    
    g_data2$group_gfam <- paste(g_data2$group, g_data2$Gfam2, sep=":")
    
    # necessary to keep boxplot boxes at correct width, psuedo data at 1,000,000
    combos <- unique(g_data2$group_gfam)
    g_data_new <- data.frame()
    for (i in 1:length(combos)){
      checker <- subset(g_data2, group_gfam == combos[i])
      if (sum(checker$Exp) > 0){
        new_df <- checker
      } else {
        new_df <- checker[1,]
        new_df$Exp = 1000000
      }
      g_data_new <- rbind(g_data_new, new_df)
    }
    colnames(g_data_new) = colnames(g_data2)
    
    g_data_new <- subset(g_data_new, Exp > 0)
    glogs <- g_data_new$Exp[g_data_new$Exp < 1000000]
    gmax <- max(glogs)
    gmin <- min(glogs)
    
    if (input$xy_switch){
      plot_exp = ggplot(g_data_new, aes(x = Gfam2, y = Exp, fill = groups)) +
        geom_boxplot(outlier.shape=3) +
        ggtitle("Logarthimic Gene Abundance") +
        coord_cartesian(ylim = c(gmin, gmax)) +
        scale_y_log10() +
        xlab("Gene Family") +
        ylab("Log Relative Abundance") +
        guides(color=FALSE, fill = guide_legend(title = "Group")) +
        theme(#axis.text.x = element_blank(), 
          #plot.margin = unit(c(.1, .1, .1, .1), "cm"),
          plot.title = element_text(hjust = 0, size = 22),
          axis.title.y = element_text(size = 22),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 22),
          axis.text.x = element_text(size = 18, 
                                     angle = 20+5*(length(unique(g_data_new$Gfam2))), 
                                     hjust=1),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18))
      plot_exp
    } else {
      plot_exp = ggplot(g_data_new, aes(x = groups, y = Exp, fill = Gfam2)) +
        geom_boxplot(outlier.shape=3) +
        ggtitle("Logarthimic Gene Abundance") +
        coord_cartesian(ylim = c(gmin, gmax)) +
        scale_y_log10() +
        xlab("Groups") +
        ylab("Log Relative Abundance") +
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
    filename = function() { paste("gene_family_abundance", '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = expression_plot(), device = 'png', 
             width = 35, height = 20, units = "cm")
    }
  )
  
  # Significance Table for Gene search t-tests
  #
  #change text annoations to "<0.0001"
  expression_table <- reactive({
    g_data2 = data.frame(exp_plot())
    g_data2 = g_data2[order(g_data2$groups),]
    
    Gfam_uniq = unique(g_data2$Gfam)
    group_uniq = as.character(unique(g_data2$groups))
    if (length(group_uniq) > 1) {
      hio=c()
      for (i in Gfam_uniq){
        acc_split = g_data2[grep(i, g_data2$Gfam),]
        acc_split[c("Exp")][is.na(acc_split[c("Exp")])] <- 0
        hi = pairwise.t.test(acc_split$Exp, acc_split$groups, p.adjust = 'BH')
        hiP <- hi$p.value
        nas <- rep(NA, length(group_uniq))
        hiP <- rbind(nas[-1], hiP)
        hiP <- cbind(hiP, nas)
        hiP[upper.tri(hiP)] = t(hiP)[upper.tri(hiP)]
        hiP <- round(hiP, 4)
        hiP[hiP < 0.0001] <- "<0.0001"
        hiP <- data.frame(hiP)
        hio = rbind.fill(hio, hiP)
      }
      hi3 = data.frame(rep(Gfam_uniq, each = length(group_uniq)))
      colnames(hi3) = "Gfam"
      hi3$group_comparisons = c(rep(group_uniq, length(Gfam_uniq)))
      hi3 = cbind.fill(hi3, hio)
      colnames(hi3) <- c("Gfam", "Group Comparison", group_uniq)
    }
    else if (length(group_uniq)==1){
      hi3 <- data.frame(Seletion=1:length(Gfam_uniq), Gene_Family=Gfam_uniq)
      hi3
    }
    hi3
  })
  
  #actual numbers for correct heat map colors
  expression_table_orig <- reactive({
    g_data2 = data.frame(exp_plot())
    g_data2 = g_data2[order(g_data2$groups),]
    
    Gfam_uniq = unique(g_data2$Gfam)
    group_uniq = as.character(unique(g_data2$groups))
    if (length(group_uniq) > 1) {
      hio=c()
      for (i in Gfam_uniq){
        acc_split = g_data2[grep(i, g_data2$Gfam),]
        acc_split[c("Exp")][is.na(acc_split[c("Exp")])] <- 0
        hi = pairwise.t.test(acc_split$Exp, acc_split$groups, p.adjust = 'BH')
        hiP <- hi$p.value
        nas <- rep(NA, length(group_uniq))
        hiP <- rbind(nas[-1], hiP)
        hiP <- cbind(hiP, nas)
        hiP[upper.tri(hiP)] = t(hiP)[upper.tri(hiP)]
        hiP <- round(hiP, 5)
        hiP <- data.frame(hiP)
        hio = rbind.fill(hio, hiP)
      }
      hi3 = data.frame(rep(Gfam_uniq, each = length(group_uniq)))
      colnames(hi3) = "Gfam"
      hi3$group_comparisons = c(rep(group_uniq, length(Gfam_uniq)))
      hi3 = cbind.fill(hi3, hio)
      colnames(hi3) <- c("Gfam", "Group Comparison", group_uniq)
    }
    else if (length(group_uniq)==1){
      hi3 <- data.frame(Seletion=1:length(Gfam_uniq), Gene_Family=Gfam_uniq)
      hi3
    }
    hi3
  })
  
  output$exp_table = renderTable({
    expression_table()
  })
  
  output$expression_table_download <- downloadHandler(
    filename = function() { paste("gene_family_abundance_table", '.txt', sep='') },
    content = function(file) {
      exp = expression_table()
      exp2 = data.frame(exp)
      write.table(exp2, file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
  
  output$exp_heat <- renderUI({
    g_data2 = data.frame(exp_plot())
    g_data2 = g_data2[order(g_data2$groups),]
    
    Gfam_uniq = unique(g_data2$Gfam)
    group_uniq = as.character(unique(g_data2$groups))
    exp_table <- expression_table_orig()
    exp_labs <- expression_table()
    if (length(group_uniq) > 1) {
      get_exp_heat_list(max_plots, length(Gfam_uniq), exp_table, exp_labs)
    }
  })
  
  output$exp_heat_download <- renderUI({
    g_data2 = data.frame(exp_plot())
    g_data2 = g_data2[order(g_data2$groups),]
    
    Gfam_uniq = unique(g_data2$Gfam)
    group_uniq = as.character(unique(g_data2$groups))
    if (length(group_uniq) > 1) {
      lapply(1:length(Gfam_uniq), function(i) {
        display_name = Gfam_uniq
        downloadButton(paste0("downloadExp", i), paste("Download", display_name[i], sep=" "))
      })
    }
  })
  
  #download gene search t-test heat maps with variable width according to Gene name
  observe({
    g_data2 = data.frame(exp_plot())
    g_data2 = g_data2[order(g_data2$groups),]
    
    Gfam_uniq = unique(g_data2$Gfam)
    group_uniq = as.character(unique(g_data2$groups))
  
    
    lapply(1:length(Gfam_uniq), function(i) {
      if (nchar(as.character(Gfam_uniq[i])) > 30){
        fixed_width <- 0.4*nchar(as.character(Gfam_uniq[i]))
      } else {fixed_width <- 15}
      
      if (nchar(as.character(Gfam_uniq[i])) > 30){
        fixed_height <- 0.35*nchar(as.character(Gfam_uniq[i]))
      } else {fixed_height <- 15}
      
      output[[paste0("downloadExp", i)]] <- downloadHandler(
        filename = function() { paste(Gfam_uniq[i], "_p_value_heat", '.png', sep='') },
        content = function(file) {
          png(file, width = fixed_width, 
              height = fixed_height, units ='cm', res = 300)
          exp_table <- expression_table_orig()
          exp_labs <- expression_table()
          download_exp_heat_list(max_plots, i, exp_table, exp_labs)
          dev.off()
        }
      )
    })
  })
  
  ## Plots for Taxa! ##

  # Taxa explore #  
  
  # matrix reordered by group and collapsed by species
  reorder_spec_mat <- reactive({
    req(grouped_samps())
    exprs_reorder = spec_nums()[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    return(exp2)
  })
  
  # Explore Taxa dataframe
  spec_taxa_data <- reactive({
    # validate(
    #   need(input$taxaSep, "Please provide a Character Delimiter")
    # )
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "Please visit the group tab to verify group assignment")
      )
    } else {
      validate(
        need(input$file1, "Please provide a file in the Upload Tab") %then%
        need(unlist(grouped_samps()), "Please select samples in the Group Tab"))
    }
    
    spec_plot_df = reorder_spec_mat()
    specs = spec_plot_df
    
    specs <- cbind(rownames(specs), specs)
    colnames(specs)[1] <- "Species"

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
    
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)
    
    spec_spec <- rep(specs$Species, each = spec_col)
    
    spec_g_data = data.frame(spec_datas, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Relative_Abundance", "Sample_num", "Group", "Species")
    
    taxa_sub_pattern = paste0("\\", input$taxaSep, ".*")
    spec_g_data$Genus = unlist(gsub(taxa_sub_pattern,"",as.character(spec_g_data$Species))) 
    spec_g_data$Genus <- str_trunc(spec_g_data$Genus, 70, "center")
    
    spec_g_data$Species = unlist(gsub(".*-","",as.character(spec_g_data$Species)))
    spec_g_data$Species <- str_trunc(spec_g_data$Species, 70, "center")
    
    spec_g_data
  })
  
  # Explore Taxa Genus Plot object
  taxa_explore_plot <- reactive({
    spec_g_data <- spec_taxa_data()
    uniq_genus_num = length(unique(spec_g_data$Genus))
    
    spec_g_data$Genus <- reorder(spec_g_data$Genus, -spec_g_data$Relative_Abundance)
    
    spec_g_data2 <- aggregate(spec_g_data$Relative_Abundance, list(spec_g_data$Genus), sum)
    
    summer <- sum(spec_g_data2[,2])
    spec_g_data2$prop <- spec_g_data2[,2]/summer
    spec_g_data3 <- subset(spec_g_data2, prop < input$upper_limit)
    spec_g_data4 <- subset(spec_g_data3, prop > input$lower_limit)
    
    keep_taxa <- as.character(unlist(spec_g_data4[,1]))
    spec_g_data_filt <- subset(spec_g_data, Genus %in% keep_taxa)
    
    exp_plot = ggplot(spec_g_data_filt, aes(x = Group, y = Relative_Abundance, fill = Genus)) +
      stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
      scale_fill_manual(values = randomColor(uniq_genus_num)) +
      ylab("Relative Abundance")
    exp_plot
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
      genera <- length(unique(spec_taxa_data()$Genus))
      png(file, width = 40, height = 0.25*genera, units ='cm', res = 300)
      par(mfrow = c(1,1))
      grid.draw(taxa_explore_legend())
      dev.off()
    }
  )
  
  output$genus_raw_data <- downloadHandler(
    filename = function() { paste("taxa_explore_raw_data", '.txt', sep='') },
    content = function(file) {
      spec_g_data <- spec_taxa_data()
      uniq_spec_num = length(unique(spec_g_data$Genus))
      
      spec_g_data$Genus <- reorder(spec_g_data$Genus, -spec_g_data$Relative_Abundance)
      write.table(spec_g_data, file, quote=F, sep='\t', row.names = F)
    }
  )
  
  # Explore Taxa Species Plot Object
  species_explore_plot <- reactive({
    spec_g_data <- spec_taxa_data()
    uniq_spec_num = length(unique(spec_g_data$Species))
    
    spec_g_data$Species <- reorder(spec_g_data$Species, -spec_g_data$Relative_Abundance)
    
    spec_g_data2 <- aggregate(spec_g_data$Relative_Abundance, list(spec_g_data$Species), sum)
    
    summer <- sum(spec_g_data2[,2])
    spec_g_data2$prop <- spec_g_data2[,2]/summer
    spec_g_data3 <- subset(spec_g_data2, prop < input$upper_limit)
    spec_g_data4 <- subset(spec_g_data3, prop > input$lower_limit)
    
    keep_taxa <- as.character(unlist(spec_g_data4[,1]))
    spec_g_data_filt <- subset(spec_g_data, Species %in% keep_taxa)

    exp_plot = ggplot(spec_g_data_filt, aes(x = Group, 
                                       y = Relative_Abundance, 
                                       fill = Species)) +
      stat_summary(fun.y = "mean", geom = "bar", position = "fill")+
      scale_fill_manual(values = randomColor(uniq_spec_num)) +
      ylab("Relative Abundance")
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
      species = length(unique(spec_taxa_data()$Species))
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
      species = length(unique(spec_taxa_data()$Species))
      png(file, width = 45, height = 0.25*species, units ='cm', res = 300)
      grid.draw(species_explore_legend())
      dev.off()
    }
  )
  
  output$species_raw_data <- downloadHandler(
    filename = function() { paste("species_explore_raw_data", '.txt', sep='') },
    content = function(file) {
      spec_g_data <- spec_taxa_data()
      uniq_spec_num = length(unique(spec_g_data$Species))
      
      spec_g_data$Species <- reorder(spec_g_data$Species, -spec_g_data$Relative_Abundance)
      write.table(spec_g_data, file, quote=F, sep='\t', row.names = F)
    }
  )
  
  # UI elements if additional taxa layer present
  output$TaxaDimExp <- renderUI({
    if (input$taxaDims > 1){
      validate(
        need(input$taxaSep, "Please provide a Character Delimiter")
      )
      mainPanel(fluidRow(plotlyOutput("taxa_explore")),
                fluidRow(downloadButton("taxa_download", "Download Plot"),
                         downloadButton("taxa_legend_download", "Download Legend"),
                         downloadButton("genus_raw_data", "Download Table")),
                width = 12)
    }
  })
  
  output$ex_delimiter <- renderText({
    paste("For example in the following string...",
          "g__Acinetobacter.s__Acinetobacter_sp_NIPH_284",
          "The genus Acinetobacter is separated from the species sp_NIPH_284 by a period.",
          "Therefore the delimiter is a period", sep='\n')
  })
  
  #### Taxa Search UI elements! ####
  
  # Taxa Search Plot
  
  #reorder full file with group allocations
  reorder_full_mat <- reactive({
    req(grouped_samps())
    full_nums <- full_file()[,4:ncol(full_file())]
    exprs_reorder = full_nums[,c(unlist(grouped_samps()))]
    exp2 = data.frame(exprs_reorder)
    return(exp2)
  })
  
  # Taxa Search Dataframe
  spec <- reactive({
    validate(
      need(input$acc_list, 'Please Select at least one Gene Family in the Gene Search Tab')
    )

    spec_list = input$acc_list
    acc_column = as.vector(acc_select())
    gfam_df = data.frame(full_file())
    
    print(dim(full_file()))
    print(dim(reorder_full_mat()))
    
    spec_plot_df = reorder_full_mat()
    spec_plot_df$Acc = gfam_df$Acc
    spec_plot_df$Gene_Family = paste(" ", gfam_df$Gene_Family, " ", sep='')
    spec_plot_df$Species = gfam_df$Species

    specs = spec_plot_df[which(spec_plot_df$Gene_Family %in% spec_list), ]
    
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

    spec_Accs = as.vector(t(specs[1+spec_col]))
    spec_Acc = rep(spec_Accs, each = spec_col)
    spec_sample_num = paste(1:(spec_col))
    spec_samp=strtoi(spec_sample_num)
    spec_isamp = rep(spec_samp,spec_row)

    spec_specs = as.vector(t(specs[3+spec_col]))
    spec_spec = rep(spec_specs, each = spec_col)

    spec_gfams0 = as.vector(t(specs[2+spec_col]))
    spec_gfam = rep(spec_gfams0, each = spec_col)
    

    spec_g_data = data.frame(spec_datas, spec_Acc, spec_gfam, 
                             spec_isamp, group_titles2, spec_spec)
    colnames(spec_g_data) = c("Exp","Acc", "Gfam", "Sample_num", "Groups", "Species")
    
    spec_g_data
  })
  
  # Generate necessary color for Taxa Search Plot
  spec_colors <- reactive({
    spec_df <- spec()
    most_specs <- max(table(spec_df$Gfam))
    cols <- randomColor(most_specs)
    cols
  })
  
  # Generate necessary UI elements for correct number of Taxa Search Plots
  observeEvent({
    input$acc_list
    new_group_names()
    }, {
    if (input$testme) {}
    else{
      validate(
        need(input$file1, "Please provide a file in the Upload Tab")
        )
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
          png(file, width = 25+0.3*specs_nums, height = 30, units ='cm', res = 300)
          plotter<-download_spec_plot_list(max_plots, 
                                  1,
                                  #subset(spec, spec$Gfam==input$acc_list[i]),
                                  spec_curr,
                                  spec_colors())
          print(plotter)
          dev.off()
        }
      )
    })
  })
  
  
  #### Correlation Plot #####
  
  ## Correlation selections
  sig_tab <- reactive({
    gfam_DF = acc_full()
    samp_paths = gfam_DF
    #samp_paths = gfam_DF[apply(gfam_DF==0,1,sum)<=(0.90*length(group_names())),]
    samp_paths$Gene_Family <- paste(" ", samp_paths$Gene_Family, " ", sep="")
    samp_paths$Gene_Family
  })
  
  ## Update selections ##
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
  
  # All sample Correlation Plot #
  
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
    samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
    
    
    heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene_Family), ]
    
    col = (ncol(samp_paths))-1
    row = nrow(heyo)
    
    heyo_small <- heyo
    rownames(heyo_small) = heyo$Gene_Family
    heyo_small = heyo_small[,-1]
    
    heyo_side = data.frame(t(heyo_small))
    colnames(heyo_side) = rownames(heyo_small)
    
    print(heyo$Gene_Family)
    
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
  
  actual_corr_names <- reactive({
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
    samp_paths$Gene_Family = paste(" ", samp_paths$Gene_Family, " ", sep="")
    
    
    heyo = samp_paths[grep(paste(corr_list, collapse = '|'), samp_paths$Gene_Family), ]
    actual_corr_names <- heyo$Gene_Family
    actual_corr_names
  })
  
  
  # Generate list of correlation matrices for each Group
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

  # Generate list of correlation significance symbols for each group
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
      validate(
        need(input$file1, "") %then%
        need(unlist(grouped_samps()), ""))
    }
    
    validate(
      need(length(input$sig_select) > 1, 'Select at least two Gene Families!')
    )
    
    Gene_Families = actual_corr_names()
    nums = length(Gene_Families)
    Label = paste(1:nums)
    Label_Matrix = data.frame(Gene_Families, Label)
    #print(Label_Matrix)
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
    new_group_names()
    input$sig_select
    }, {
    if (input$testme) {
      validate(
        need(length(group_dims())==4, "")
      )
    }
    else{
      validate(
        need(input$file1, "") 
      )
    }
    
    if (input$numInputs > 1) {
      #checker = input$sig_select
      #display_name = group_names()
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
          png(file, width = 25, height = 20, units ='cm', res = 300)
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