library(ape)

args <- commandArgs(TRUE)
outdir <- args[1]
tr <- read.tree(args[2])
jpeg(paste(outdir, 'tree.jpg', sep='/'))
plot.phylo(tr)
dev.off()

