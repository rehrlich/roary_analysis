library(stringr)
library(gplots)
source('plotTree.R')


main <- function(){
  args <- commandArgs(TRUE)
  outdir <- args[1]
  nickname <- args[2]
  my.tree <- args[3]

   for (table_type in c('similarity', 'difference', 'pair_unique')){
    prefix <- paste(outdir, nickname, "_pairwise_", table_type, sep="")
    in.file <- paste(prefix, "_table.txt", sep="")
    out.file <- paste(prefix, "_heatmap.pdf", sep="")
    pairwise <- read.table(in.file, row.names=1, header = T, stringsAsFactors=F,
                           sep="\t", fill = T)

    missing.col <- row.names(pairwise)[1]
    missing.row <- colnames(pairwise)[1]
    pairwise[missing.row, ] <- NA
    pairwise[, missing.col] <- NA

    for(r in row.names(pairwise)){
      for(c in colnames(pairwise)){
        if(is.na(pairwise[r, c]) & r != c){
          pairwise[r, c] <- pairwise[c, r]
        }
      }
    }

    pairwise[is.na(pairwise)] = 0

    #colnames(pairwise)[names(pairwise) == ref] <- paste(ref, ".ref", sep="")
    #rownames(pairwise)[rownames(pairwise) == ref] <- paste(ref, ".ref", sep="")


    my.heatmap <- plotTree(tree=my.tree,heatmapData=pairwise,ancestral.reconstruction=F,
                           tip.colour.cex=1,cluster=T,tip.labels=T,
                           tipColours=c("black","purple2","skyblue2","grey"),lwd=1,
                           colourNodesBy="location",treeWidth=5,
                           heatmap.colours=c("white","grey","seagreen3","darkgreen"),
                           dataWidth=20,infoCols=NA,outputPDF=out.file)

    st.order <- my.heatmap$strain_order

    pairwise <- pairwise[st.order, ]
    pairwise <- pairwise[, st.order]
    pairwise <- apply(pairwise, 2, rev)
    pairwise <- apply(pairwise, 1, rev)

    num.strains <-dim(pairwise)[1]
    for(r in 1:num.strains){
      for(c in 1:r){
        pairwise[r, c] <- NA
      }
    }

    my.heatmap2 <- plotTree(tree=my.tree,heatmapData=pairwise,ancestral.reconstruction=F,
                            tip.colour.cex=1,cluster=F,tip.labels=T,
                            lwd=1,
                            colourNodesBy="location",treeWidth=5,
                            heatmap.colours=c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"),
                            dataWidth=20,infoCols=NA,outputPDF=out.file)

  }

}
main()