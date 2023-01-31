#' @title Convert \code{\link{phyloseq-class}} object to long data format
#' @description An alternative to psmelt function from 
#'              \code{\link{phyloseq-class}} object.
#' @param x \code{\link{phyloseq-class}} object
#' @param x \code{\link{phyloseq-class}} object
#' @param sample.column A single character string specifying name
#'                      of the column to hold sample names.
#' @param feature.column A single character string specifying name
#'                      of the column to hold OTU or ASV names.
#'                      
#' @examples
#' data("dietswap")
#' ps.melt <- psmelt2(dietswap, sample.column="SampleID", 
#'                    feature.column="Feature") 
#' head(ps.melt)                                         
#' @return A \code{tibble} in long format
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
#' @export                     
png.merge_pseq_list <- function(pseq_list, label="Site"){
  
  
  
  if(FALSE){
    load("./data/pseq_list.rda")
    # pseq_list: list(data=pseq_total_phylum, joint=jpseq.jive, individual=apseq.jive)
    label="Site"
  }
  
  N <- length(pseq_list)
  
  if( !is.null(names(pseq_list)) ){
    Name <- names(pseq_list)
  } else {
    Name <- LETTERS[seq_len(N)]
  }
  
  
  otu_list <- lapply(pseq_list, otu_table)
  tax_list <- lapply(pseq_list, tax_table)
  sample_list <- lapply(pseq_list, sample_data)
  
  # taxa_all <- unique(do.call("c", lapply(otu_list, taxa_names)))
  
  
  # start <- proc.time()
  otu_new <- lapply(1:N, function(idx){
    x_new <- otu_list[[idx]] %>% { cbind.data.frame(Taxa=rownames(.), .) }
    colnames(x_new)[-1] <- paste0(Name[idx],"_",colnames(x_new)[-1])
    x_new
  })
  tax_new <- lapply(1:N, function(idx){
    tax_list[[idx]] %>% { cbind.data.frame(Taxa=rownames(.), .) }
  })
  
  otu_all <- otu_new %>% reduce(full_join, by="Taxa")
  otu_all[is.na(otu_all)] <- 0
  tax_all <- tax_new %>% reduce(full_join)
  tax_all[is.na(tax_all)] <- 0
  
  # end <- proc.time()
  # end - start
  
  
  
  for( idx in 1:N ){
    sample_new <- cbind.data.frame(Site=names(sample_list)[idx], sample_list[[idx]])
    colnames(sample_new) <- c(label, colnames(sample_list[[idx]]))
    
    if( idx == 1 ){
      sample_all <- sample_new
    } else {
      sample_all <- rbind.data.frame( sample_all, sample_new )
    }
  }
  
  rownames(otu_all) <- otu_all[,1]
  otu_all <- as.matrix(otu_all[,-1])
  # colnames(otu_all) <- gsub("[a-z]+\\_", "", colnames(otu_all))
  
  my_sample_data <- sample_data( sample_all )
  my_taxonomyTable <- tax_table( as.matrix(tax_all[,-1]) )
  my_otu_table <- otu_table( otu_all, taxa_are_rows = TRUE )
  
  rownames(my_taxonomyTable) <- rownames(my_otu_table)
  my_sample_data@row.names <- colnames(my_otu_table)
  
  my_phyloseq <- phyloseq(my_otu_table, my_taxonomyTable, my_sample_data)
  
  my_phyloseq
}















make_barplot2 <- function (dfm, group_by, direction="vertical") {
  
  Tax <- Sample <- Abundance <- NULL
  
  # Provide barplot
  dfm <- dfm %>% arrange(Tax)  # Show Taxs always in the same order
  dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
  
  p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) +
    geom_bar(position="stack", stat="identity") +
    scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample)
  
  # Name appropriately
  p <- p + labs(y = "Abundance")
  
  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5,
                                          hjust=0))
  p <- p + guides(fill=guide_legend(reverse=FALSE))
  
  if (!is.null(group_by)) {
    if( direction == "vertical" ){
      p <- p + facet_grid(Group1~Group2)
    } else if( direction == "horizontal" ) {
      p <- p + facet_grid(Group2~Group1, drop = TRUE,
                          space = "free", scales = "free") 
    }
    
    
  }
  list(df=dfm, plot=p)
}








#' @export 
png.top_taxa_list <- function(LIST, cutoff=0.0001, n=NULL){
  
  if( is.null(n) ){
    n <- 0.5*length(taxa(LIST[[1]]))
  }
  
  lapply( LIST, function(x) taxa_sums(x) / ncol(x@otu_table) ) %>% do.call("c", .) %>% sort(decreasing=TRUE) %>% {.[.>cutoff]} %>% {names(.)[!duplicated(names(.))]} %>% head(n)
  
}




#' @export png.PlotComposition
png.PlotComposition <- function(pseq, taxa=NULL, group_by="Site", group.level=c("urine", "serum", "stool", "stoolp", "stools"), sample.pattern="[a-z]+\\_", filename="./plot.pdf", height=5, width=10, legend.ncol=10){
  
  
  # tmp.list[[idx]] %>% { prune_taxa( taxa_sums(.)/N > 0.001, . ) } %>% taxa
  # tmp.list[[idx]] %>% top_taxa(14)
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  TOP_taxa <- pseq@tax_table[top_taxa(pseq, n = 1), "unique"]
  
  fit.comp <- pseq %>% 
    plot_composition(sample.sort = TOP_taxa, otu.sort = "abundance", group_by=group_by, group.level=group.level, direction="vertical", sample.pattern=sample.pattern)
  
  fit.comp$df %>% head %>% print
  
  if( is.null(taxa) ){
    taxa <- levels(fit.comp$df$Tax)
  }
  
  if( length(taxa) > 74 ) warnings("The maximum number of colors is 74.")
  
  
  p <- fit.comp$plot +
    scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) + 
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=16) +
    labs(x = "Samples", y = "Relative Abundance",
         title = "Relative Abundance data") + 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          legend.position = "bottom") + 
    guides(fill = guide_legend(ncol = legend.ncol))
  
  
  
  pdf(file=filename, height=height, width=width)
  print(p)
  dev.off()
  
  
  invisible(fit.comp)
  
}









#' @export png.PlotComposition2
png.PlotComposition2 <- function(pseq, taxa=NULL, group_by=c("Structure", "Site"), group.level=list( c("data", "joint", "individual"), c("urine", "serum", "stools") ), filename="./plot.pdf", sample.sort="top", height=5, width=10, legend.ncol=10){
  
  if(FALSE){
    pseq <- decomp_ajive
    level="Phylum"
    filename="./plot.pdf"
    height=5
    width=10
    legend.ncol=10
    sample.sort="top"
  }
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  
  if(sample.sort == "top"){
    TOP_taxa <- pseq@tax_table[top_taxa(pseq, n = 1), "unique"]
  } else {
    TOP_taxa <- NULL
  }
  
  
  
  if( is.null(taxa) ){
    fit.comp0 <- pseq %>% phyloseq::subset_samples(Structure == "data") %>% 
      plot_composition(sample.sort = TOP_taxa, otu.sort = "abundance", group_by=c("Site"), group.level=c("urine", "serum", "stool", "stoolp", "stools"), group.type="vertical", sample.pattern="[a-z]+\\_[a-z]+\\_")
    
    taxa <- levels(fit.comp0$df$Tax)
  }
  
  
  fit.comp <- pseq %>% 
    plot_composition2(sample.sort = TOP_taxa, otu.sort = "abundance", group_by=group_by, group.level=group.level, group.type="vertical", sample.pattern="[a-z]+\\_[a-z]+\\_")
  
  
  fit.comp$df %>% head %>% print
  
  
  if( length(taxa) > 74 ) warnings("The maximum number of colors is 74.")
  
  
  p <- fit.comp$plot +
    scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) + 
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=16) +
    labs(x = "Samples", y = "Relative Abundance",
         title = "Relative Abundance data") + 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          legend.position = "bottom") + 
    guides(fill = guide_legend(ncol = legend.ncol))
  
  
  
  pdf(file=filename, height=height, width=width)
  print(p)
  dev.off()
  
  
  invisible(fit.comp)
  
}











plot_composition2 <- function(x,
                              sample.sort=NULL,
                              otu.sort=NULL,
                              x.label="sample",
                              plot.type="barplot",
                              verbose=FALSE, 
                              average_by=NULL,
                              group_by = NULL, 
                              group.level = NULL,
                              direction = "vertical", 
                              sample.pattern = NULL, ...) {
  
  if(FALSE){
    x <- decomp_jive
    sample.sort = TOP_taxa
    otu.sort = "abundance"
    group_by = c("Structure", "Site")
    group.level = list( c("data", "joint", "individual"), c("urine", "serum", "stool", "stoolp", "stools") )
    direction = "vertical"
    sample.pattern="[a-z]+\\_"
    
    x.label="sample"
    plot.type="barplot"
    verbose=FALSE
    average_by=NULL
    
  }
  
  
  # Avoid warnings
  Sample <- Abundance <- Taxon <- Group <- Tax <-
    horiz <- value <- scales <- ID <- 
    meta <- OTU <- taxic <- otu.df <- taxmat <-  new.tax<- NULL
  if (!is.null(x@phy_tree)){
    x@phy_tree <- NULL
  }
  
  xorig <- x
  if (verbose) {message("Pick the abundance matrix taxa x samples")}
  abu <- abundances(x)
  if (verbose) {message("Average the samples by group")}
  group <- NULL
  if (!is.null(average_by)) {
    dff <- as.data.frame(t(abu))
    dff$group <- sample_data(x)[[average_by]]
    if (is.numeric(dff$group) || is.character(dff$group)) {
      dff$group <- factor(dff$group, levels=sort(unique(dff$group)))
    }
    # Remove samples with no group info
    dff <- dff %>% filter(!is.na(group))
    dff$group <- droplevels(dff$group)
    # av <- ddply(dff, "group", colwise(mean))
    av <- aggregate(. ~ group, data = dff, mean)
    rownames(av) <- as.character(av$group)
    av$group <- NULL
    abu <- t(av)  # taxa x groups
  }
  
  if (verbose) {message("Sort samples")}
  if (is.null(sample.sort) || sample.sort == "none" ||
      !is.null(average_by)) {
    # No sorting sample.sort <- sample_names(x)
    sample.sort <- colnames(abu)
  } else if (length(sample.sort) == 1 && sample.sort %in% taxa(xorig)) {
    tax <- sample.sort
    sample.sort <- rev(sample_names(x)[order(abundances(x)[tax,])])
  } else if (length(sample.sort) == 1 &&
             sample.sort %in% names(sample_data(x)) && 
             is.null(average_by)) {
    # Sort by metadata field
    sample.sort <- rownames(sample_data(x))[order(
      sample_data(x)[[sample.sort]])]
  } else if (all(sample.sort %in% sample_names(x)) & is.null(average_by)) {
    # Use predefined order
    sample.sort <- sample.sort
  } else if (length(sample.sort) == 1 && sample.sort == "neatmap") {
    sample.sort <- neatsort(x, method="NMDS", distance="bray",
                            target="sites", first=NULL)
  } else if (is.vector(sample.sort) && length(sample.sort) > 1) {
    sample.sort <- sample_names(x)[sample.sort]            
  } else if (!sample.sort %in% names(sample_data(x))) {
    warning(paste("The sample.sort argument", sample.sort,
                  "is not included in sample_data(x). 
            Using original sample ordering."))
    sample.sort <- sample_names(x)
  }
  
  # Sort taxa
  if (is.null(otu.sort) || otu.sort == "none") {
    # No sorting
    otu.sort <- taxa(x)
  } else if (length(otu.sort) == 1 && otu.sort == "abundance2") {
    otu.sort <- rev(c(rev(names(sort(rowSums(abu)))[seq(1, nrow(abu), 2)]),
                      names(sort(rowSums(abu)))[seq(2, nrow(abu), 2)]))
  } else if (length(otu.sort) == 1 && otu.sort == "abundance") {
    otu.sort <- rev(names(sort(rowSums(abu))))
  } else if (length(otu.sort) == 1 && otu.sort %in%
             colnames(tax_table(x))) {
    otu.sort <- rownames(sample_data(x))[order(tax_table(x)[[otu.sort]])]
  } else if (all(otu.sort %in% taxa(x))) {
    # Use predefined order
    otu.sort <- otu.sort
  } else if (length(otu.sort) == 1 && otu.sort == "neatmap") {
    otu.sort <- neatsort(x, method="NMDS", distance="bray",
                         target="species", first=NULL)
  }
  
  # Abundances as data.frame dfm <- psmelt(x)
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("Tax", "Sample", "Abundance")
  dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels=otu.sort)
  
  if (!is.null(group_by)) {
    if (!is.null(average_by)) {
      dfm$Group1 <- meta(x)[[group_by[1]]][match(as.character(dfm$Sample),
                                                 meta(x)[[average_by]])]
      dfm$Group2 <- meta(x)[[group_by[2]]][match(as.character(dfm$Sample),
                                                 meta(x)[[average_by]])]
    }else{
      dfm$Group1 <- meta(x)[[group_by[1]]][match(as.character(dfm$Sample),
                                                 sample_names(x))]
      dfm$Group2 <- meta(x)[[group_by[2]]][match(as.character(dfm$Sample),
                                                 sample_names(x))]
    }
  }
  
  if(!is.null(group.level)){
    dfm$Group1 <- factor(dfm$Group1, levels=group.level[[1]])
    dfm$Group2 <- factor(dfm$Group2, levels=group.level[[2]])
  }
  
  # SampleIDs for plotting
  if (x.label %in% colnames(sample_data(x)) & is.null(average_by)) {
    
    meta <- sample_data(x)
    dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample), x.label]))
    
    # Sort the levels as in the original metadata
    if (is.factor(meta[, x.label])) {            
      lev <- levels(meta[, x.label])            
    } else {            
      lev <- unique(as.character(unname(unlist(meta[, x.label]))))
    }       
    dfm$xlabel <- factor(dfm$xlabel, levels=lev)        
  } else {       
    dfm$xlabel <- dfm$Sample        
  }
  
  
  if( !is.null(sample.pattern) ){
    dfm$Sample <- gsub(sample.pattern, "", dfm$Sample)
    dfm$xlabel <- gsub(sample.pattern, "", dfm$xlabel)
  }
  
  
  if (verbose) {message("Construct the plots")}   
  if (plot.type == "barplot") {
    p <- make_barplot2(dfm, group_by, direction=direction)
  } else if (plot.type == "heatmap") {
    p <- make_heatmap1(x, otu.sort, sample.sort, verbose)        
  } else if (plot.type == "lineplot") {
    p <- make_lineplot1(dfm) 
  } else {
    stop("plot.type argument not recognized")
  }
  
  p
  
}











#' @export png.PlotComposition3
png.PlotComposition3 <- function(pseq, taxa=NULL, sample.sort="unique", filename="./plot.pdf", height=5, width=10, legend.ncol=10){
  
  # png.PlotComposition2(pseq_comp, taxa=TopN.Taxa, sample.sort=level,
  #                      filename=paste0("./tmp_", idx, ".",level,".pdf"),
  #                      height=10, width=21, legend.ncol=10)
  
  # tmp.list[[idx]] %>% { prune_taxa( taxa_sums(.)/N > 0.001, . ) } %>% taxa
  # tmp.list[[idx]] %>% top_taxa(14)
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  TOP_taxa <- pseq@tax_table[top_taxa(pseq, n = 1), sample.sort]
  
  fit.comp <- pseq %>%
    plot_composition(sample.sort = TOP_taxa, otu.sort = "abundance")
  
  # fit.comp$df %>% head %>% print
  
  if( is.null(taxa) ){
    taxa <- levels(fit.comp$df$Tax)
  }
  
  if( length(taxa) > 74 ) warnings("The maximum number of colors is 74.")
  
  
  p <- fit.comp$plot +
    scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) +
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=16) +
    labs(x = "Samples", y = "Relative Abundance",
         title = "Relative Abundance data") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = legend.ncol))
  
  
  
  pdf(file=filename, height=height, width=width)
  print(p)
  dev.off()
  
  
  invisible(fit.comp)
  
}
