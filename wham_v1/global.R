## global functions

library(ggplot2)
library(rowr)
library(plotly)
library(grid)
library(gridExtra)
library(data.table)
library(randomcoloR)
library(gplots)
library(viridis)
library(matrixStats)
library(webshot)
library(ALDEx2)
library(DT)



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
    
    dim_fact <- dim(plot_mat)[1]
    if (dim_fact > 10) {
      fill_size <- 2-(0.8*dim_fact)
    } else if (dim_fact > 20) {
      fill_size <- 2-(0.12*dim_fact)
    } else {
      fill_size <- 2-(0.1*dim_fact)
    }
    
    heatmap.2(plot_mat,
              cellnote = sym_mat, notecol="black", notecex = fill_size,
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

### heatmap.3 source code courtesy of obigriffith
## https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 <-
  function(x,
           Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
           distfun = dist,
           hclustfun = hclust,
           dendrogram = c("both","row", "column", "none"),
           symm = FALSE,
           scale = c("none","row", "column"),
           na.rm = TRUE,
           revC = identical(Colv,"Rowv"),
           add.expr,
           breaks,
           symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
           col = "heat.colors",
           colsep,
           rowsep,
           sepcolor = "white",
           sepwidth = c(0.05, 0.05),
           cellnote,
           notecex = 1,
           notecol = "cyan",
           na.color = par("bg"),
           trace = c("none", "column","row", "both"),
           tracecol = "cyan",
           hline = median(breaks),
           vline = median(breaks),
           linecol = tracecol,
           margins = c(5,5),
           ColSideColors,
           RowSideColors,
           side.height.fraction=0.3,
           cexRow = 0.2 + 1/log10(nr),
           cexCol = 0.2 + 1/log10(nc),
           labRow = NULL,
           labCol = NULL,
           key = TRUE,
           keysize = 1.5,
           density.info = c("none", "histogram", "density"),
           denscol = tracecol,
           symkey = max(x < 0, na.rm = TRUE) || symbreaks,
           densadj = 0.25,
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           lmat = NULL,
           lhei = NULL,
           lwid = NULL,
           ColSideColorsSize = 1,
           RowSideColorsSize = 1,
           KeyValueName="Value",...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        #max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        max.raw = 1
        min.raw = 0
        #min.raw <- -max.raw
        #tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        #tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        tmpbreaks[1] <- 0
        tmpbreaks[length(tmpbreaks)] <- 1
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }



scale_fill_viridis <- function (..., alpha=1, discrete=TRUE) {
  if (discrete) {
    discrete_scale("fill", "viridis", viridis_pal(alpha), ...)
  } else {
    scale_fill_gradientn(colours = viridis(256, alpha), ...)
  }
}



##plotly
custom_export <- function (p = last_plot(), file = "plotly.png", selenium = NULL, 
          ...) 
{
  fileType <- tolower(tools::file_ext(file))
  if (!fileType %in% c("jpeg", "png", "webp", "svg", "pdf")) {
    stop("File type ", fileType, " not supported", call. = FALSE)
  }
  if (is.webgl(p) && fileType %in% "pdf") {
    stop("A personal (or professional) plan is required to export WebGL to pdf:\\n", 
         "https://plot.ly/products/cloud/", call. = FALSE)
  }
  use_webshot <- !is.webgl(p) && fileType %in% c("jpeg", "png", 
                                                 "pdf")
  if (!use_webshot) {
    cmd <- sprintf("function(el, x) {\\n        var gd = document.getElementById(el.id); \\n        Plotly.downloadImage(gd, {format: '%s', width: %s, height: %s, filename: '%s'});\\n      }", 
                   fileType, p$width %||% p$layout$width %||% 800, p$height %||% 
                     p$layout$height %||% 600, tools::file_path_sans_ext(file))
    p <- htmlwidgets::onRender(p, cmd)
  }
  f <- basename(tempfile("plotly", ".", ".html"))
  on.exit(unlink(f), add = TRUE)
  html <- htmlwidgets::saveWidget(p, f)
  if (use_webshot) {
    try_library("webshot", "export")
    return(webshot::webshot(f, file, ...))
  }
  if (inherits(selenium, "rsClientServer")) {
    selenium$client$navigate(paste0("file://", normalizePath(f)))
  }
  else {
    stop("Must provide an object of class 'rsClientServer' to the `selenium` ", 
         "argument to export this plot (see examples section on `help(export)`)", 
         call. = FALSE)
  }
  message(sprintf("Success! Check your downloads folder for a file named: '%s'", 
                  file))
  invisible(file)
}



#
#
#
#

##### aldex functions

aldex.ttest2 <- function (clr, conditions, paired.test = FALSE, hist.plot = FALSE) 
{
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  conditions <- as.factor(conditions)
  levels <- levels(conditions)
  if (length(conditions) != numConditions(clr)) {
    stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):", 
               length(conditions), "len(names(clr)):", numConditions(clr)))
  }
  if (length(levels) != 2) 
    stop("only two condition levels are currently supported")
  levels <- vector("list", length(levels))
  names(levels) <- levels(conditions)
  sets <- names(levels)
  setAsBinary <- as.numeric(conditions == sets[1])
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  wi.p.matrix <- as.data.frame(matrix(1, nrow = feature.number, 
                                      ncol = mc.instances))
  wi.BH.matrix <- wi.p.matrix
  we.p.matrix <- wi.p.matrix
  we.BH.matrix <- wi.p.matrix
  print("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for (mc.i in 1:mc.instances) {
    numTicks <- progress(mc.i, mc.instances, numTicks)
    t.input <- sapply(mc.all, function(y) {
      y[, mc.i]
    })
    #for (j in 1:nrow(t.input)){
    #  wi.p.matrix[j, mc.i] <- wilcox.test(t.input[j,1:length(setA)], t.input[j,c(length(setA)+1):c(length(setA)+length(setB))], 
    #                                      paired = paired.test)$p.value
    #  
    #  wi.BH.matrix[j, mc.i] <- p.adjust(wi.p.matrix[j, mc.i], 
    #                                    method = "BH")
    #  
    #  #we.p.matrix[, mc.i] <- t.fast(t.input, setAsBinary, paired.test)
    #  
    #  we.p.matrix[j, mc.i] <- t.test(t.input[j,1:length(setA)], t.input[j,c(length(setA)+1):c(length(setA)+length(setB))], 
    #                                 paired = paired.test)$p.value
    #  
    #  we.BH.matrix[j, mc.i] <- p.adjust(we.p.matrix[j, mc.i], 
    #                                    method = "BH")
    #}

    for (j in 1:nrow(t.input)){
    wi.p.matrix[j, mc.i] <- wilcox.test(t.input[j,1:length(setA)], t.input[j,c(length(setA)+1):c(length(setA)+length(setB))], 
                                        paired = paired.test)$p.value
    
    #we.p.matrix[, mc.i] <- t.fast(t.input, setAsBinary, paired.test)
    
    we.p.matrix[j, mc.i] <- t.test(t.input[j,1:length(setA)], t.input[j,c(length(setA)+1):c(length(setA)+length(setB))], 
                                   paired = paired.test)$p.value
    }
  
  
    wi.BH.matrix[, mc.i] <- p.adjust(wi.p.matrix[, mc.i], method = "BH")
    we.BH.matrix[, mc.i] <- p.adjust(we.p.matrix[, mc.i], method = "BH")
  }
  if (hist.plot == TRUE) {
    par(mfrow = c(2, 2))
    hist(we.p.matrix[, 1], breaks = 99, main = "Welch's P values Instance 1")
    hist(wi.p.matrix[, 1], breaks = 99, main = "Wilcoxon P values Instance 1")
    hist(we.BH.matrix[, 1], breaks = 99, main = "Welch's BH values Instance 1")
    hist(wi.BH.matrix[, 1], breaks = 99, main = "Wilcoxon BH values Instance 1")
    par(mfrow = c(1, 1))
  }
  we.ep <- rowMeans(we.p.matrix)
  we.eBH <- rowMeans(we.BH.matrix)
  wi.ep <- rowMeans(wi.p.matrix)
  wi.eBH <- rowMeans(wi.BH.matrix)
  z <- data.frame(we.ep, we.eBH, wi.ep, wi.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}



progress <- function(i, k, numTicks){
  
  if(i == 1) numTicks <- 0
  
  if(numTicks == 0) cat("|-")
  
  while(i > numTicks*(k/40)){
    
    cat("-")
    if(numTicks == 10) cat("(25%)")
    if(numTicks == 20) cat("(50%)")
    if(numTicks == 30) cat("(75%)")
    numTicks <- numTicks + 1
  }
  
  if(i == k) cat("-|\n")
  
  return(numTicks)
}
