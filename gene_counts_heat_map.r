# Author:  Rachel Ehrlich
# This program makes heat maps from Roary's output gene counts.
# Command line inputs are the roary output folder, nickname for the run,
# and the output directory
# Example usage:
# Rscript gene_counts_heat_map.r roary_dir nickname outdir

library(stringr)
library(gplots)
library(RColorBrewer)

# Input is a string from Roary's gene_presence_absence.csv
# Output is the copy number for the corresponding gene/strain
# Copies are separated by tabs
GetCopyNumber <- function(cell) {
  cell
  if (cell == ""){
    return(0)
  } else {
    return(1 + str_count(cell, "\t"))
  }
}

# Makes a heatmap showing gene presence (1+ copies) and
# absence (0 copies) for the gene.copy.num matrix.  prefix
# is for the output file name
MakeGenePossHeatMap <- function(gene.copy.num, prefix){
  two.color.heatmap <- gene.copy.num
  two.color.heatmap[two.color.heatmap > 0] <- 1
  png(paste(prefix, 'gene_poss.png', sep = ""))

  heatmap.2(as.matrix(two.color.heatmap), col=brewer.pal(3,"BuGn"),
            trace="none", margins=c(5,15))
  graphics.off()
}

# Makes a heatmap showing gene counts for the gene.copy.num
# matrix.  prefix is for the output file name
MakeGeneCountsHeatMap <- function(gene.copy.num, prefix){
  png(paste(prefix, 'gene_counts.png', sep = ""))
  heatmap.2(as.matrix(gene.copy.num), trace="none", margins=c(5,15))
  graphics.off()
}

# Input: cols are strains, rows are genes
# Returns a matrix containing only distributed genes
FilterCoreGenes <- function(df) {

  is.core <- apply(df, 1, FUN=function(x) {min(x) == 1})
  filtered <- df[!is.core, ]
  return(filtered)
}

# Input is the roary output directory
# Returns a data frame whose cols are strains, rows are genes and entries are
# the copy number
GetGeneDataCopies <- function(roary.output){
  gene.data.full <- read.csv(paste(roary.output, "gene_presence_absence.csv",
                                   sep="/"), row.names=1, stringsAsFactors=F)

  gene.data.poss <- gene.data.full[, 11:ncol(gene.data.full)]
  gene.data.copies <- as.data.frame(lapply(gene.data.poss,
                                           FUN=function(x)
                                           {sapply(x, FUN=GetCopyNumber)}))
  row.names(gene.data.copies) <- row.names(gene.data.poss)
  return(gene.data.copies)
}

# a data frame whose cols are strains, rows are genes and entries are
# the copy number and a prefix for the output file
# Calls functions to make a heat map of gene possession and of gene counts
MakeHeatMaps <- function(gene.counts, prefix){
  gene.counts <- t(gene.counts)
  MakeGenePossHeatMap(gene.counts, prefix)
  MakeGeneCountsHeatMap(gene.counts, prefix)
}

main <- function(){
  args <- commandArgs(TRUE)
  roary.output <- args[1]
  nickname <- args[2]
  outdir <- args[3]

  gene.data.copies <- GetGeneDataCopies(roary.output)

  prefix <- paste(outdir, "/", nickname, "_heatmap_all_", sep="")
  MakeHeatMaps(gene.data.copies, prefix)

  gene.data.copies.distrib <- FilterCoreGenes(gene.data.copies)

  prefix <- paste(outdir, "/", nickname, "_heatmap_no_core_", sep="")
  MakeHeatMaps(gene.data.copies.distrib, prefix)
}

main()
