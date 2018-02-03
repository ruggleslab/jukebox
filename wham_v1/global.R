## global functions

library(ggplot2)
library(plyr)
library(rowr)
library(plotly)
library(grid)
library(base)
library(rPython)
library(data.table)
library(stringr)
library(randomcoloR)
library(gplots)
library(tidyr)
library(dplyr)

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
      spec$Species <- reorder(spec$Species, -spec$Exp)
      spec$Abundance <- spec$Exp

      p = ggplot(spec, aes(x=Groups, y=Abundance, fill = Species)) +
        #geom_bar(position = "fill", stat='identity') +
        stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
        ggtitle(paste("Species Contribution for", gfam_uniq[i], sep =' ')) +
        ylab("Relative Abundance") +
        xlab("Group") +
        theme(plot.title = element_text(size = 11),
              legend.position="right",
              #legend.key.width = unit(0.1, "cm"),
              #legend.key.height = unit(0.1, "cm"),
              legend.text = element_text(size=10),
              legend.title = element_text(size=10),
              axis.title.y = element_text(size = 11),
              axis.title.x = element_text(size = 11),
              axis.text.y = element_text(size = 11),
              axis.text.x = element_text(size = 11)) +
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
  spec$Species <- reorder(spec$Species, -spec$Exp)
  spec$Abundance <- spec$Exp
  
  p = ggplot(spec, aes(x=Groups, y=Abundance, fill = Species)) +
    #geom_bar(position = "fill", stat='identity') +
    stat_summary(fun.y = "mean", geom = "bar", position = "fill") +
    ggtitle(paste("Species Contribution for \n", spec$Gfam[1], sep =' ')) +
    ylab("Relative Abundance") +
    xlab("Group") +
    guides(fill = guide_legend(nrow = 40))+
    theme(plot.title = element_text(size = 22),
          legend.position="right",
          #legend.key.width = unit(0.1, "cm"),
          #legend.key.height = unit(0.1, "cm"),
          legend.text = element_text(size=10),
          legend.title = element_text(size=22),
          axis.title.y = element_text(size = 22),
          axis.title.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          axis.text.x = element_text(size = 22)) +
    scale_fill_manual(values = supplied_col_vector)
  # })
  # 
  # do.call(tagList, plot_output_list) # needed to display properly.
  # 
  # return(plot_output_list)
  p
}

#### exp table heatmaps!

get_exp_heat_list <- function(max_plots, input_n, input_df, lab_df) {
  # Insert plot output objects the list
  plot_output_list <- lapply(1:input_n, function(i) {
    plotname <- paste("plot", i, sep="")
    savename <- paste("save", i, sep="")
    plot_output_object <- plotOutput(plotname, height = 280, width = 300)
    plot_output_object <- renderPlot({
      spec_df <- input_df
      Gfam_uniq = as.character(unique(spec_df$Gfam))
      exp_table <- spec_df
      new_df <- subset(exp_table, Gfam == Gfam_uniq[i])
      new_mat <- data.matrix(new_df[,3:ncol(new_df)])
      
      exp_lab <- lab_df
      lab_df2 <- subset(exp_lab, Gfam == Gfam_uniq[i])
      lab_mat <- as.matrix(lab_df2[,3:ncol(lab_df2)])

      rownames(new_mat) <- new_df[,2]
      mat_title <- Gfam_uniq[i]
      par(cex.main=0.7)
      heatmap.2(new_mat,
                cellnote = lab_mat,
                notecol="black",
                density.info="none",
                main=mat_title,
                #key = TRUE,
                #keysize = 1.0,
                breaks = seq(0, 0.05, by = 0.0005),
                col=c(colorRampPalette(c("red", "white"))(n=100)),
                dendrogram = 'none',
                Rowv=F,
                Colv=F,
                margins=c(10,10),
                cexRow=1.2,
                cexCol=1.2,
                trace=c("none"),
                na.color="gray60"
      )
    })
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
}

download_exp_heat_list <- function(max_plots, input_n, input_df, lab_df) {
  # Insert plot output objects the list
  i <- input_n
  plotname <- paste("plot", i, sep="")
  savename <- paste("save", i, sep="")
  spec_df <- input_df
  Gfam_uniq = as.character(unique(spec_df$Gfam))
  exp_table <- spec_df
  new_df <- subset(exp_table, Gfam == Gfam_uniq[i])
  new_mat <- data.matrix(new_df[,3:ncol(new_df)])
  
  exp_lab <- lab_df
  lab_df2 <- subset(exp_lab, Gfam == Gfam_uniq[i])
  lab_mat <- as.matrix(lab_df2[,3:ncol(lab_df2)])

  rownames(new_mat) <- new_df[,2]
  mat_title <- Gfam_uniq[i]
  par(cex.main=0.7)
  heatmap.2(new_mat,
            cellnote = lab_mat,
            notecol="black",
            density.info="none",
            main=mat_title,
            #key = TRUE,
            #keysize = 1.0,
            breaks = seq(0, 0.05, by = 0.0005),
            col=c(colorRampPalette(c("red", "white"))(n=100)),
            dendrogram="none",
            Rowv=F,
            Colv=F,
            margins=c(10,10),
            cexRow=1.2,
            cexCol=1.2,
            trace=c("none"),
            na.color="gray60"
  )
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
  
  #makes a list of grouped samples that can be used to reorder
  #original matrix
  listcond <- reactive({
    cond2_t = input$exptallSamples
    transcription_conditions = list(cond2_t)
    return(transcription_conditions)
  })
  return(listcond)
}
