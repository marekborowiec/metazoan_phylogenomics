# We used this code to plot unrooted
# tree topologies from taxon removal analyses

setwd("/Your/Path/Here")

library("ape")

# load and concatenate trees
# m15 = 'NoOutgr' matrix
# m16 = 'NoAmphi' matrix
# m17 = 'NoOutgrAmphi' matrix
m15_tr <- read.tree("./NoOutgr/RAxML_bipartitions.part_matrix15.out")
m16_tr <- read.tree("./NoAmphi/RAxML_bipartitions.BIPS")
m17_tr <- read.tree("./NoOutgrAmphi/RAxML_bipartitions.part_matrix17.out")

trees <- c(m15_tr, m17_tr, m16_tr)

# open pdf device
cairo_pdf(filename="unrooted2.pdf",    # create SVG file       
          width = 3.5,        
          height = 6,
          pointsize = 12)  

# create plot space
par(mfrow=c(3,1))

# plot each tree with nodes colored 
# according to support values

for (tree in trees) {
  # get support values for the tree
  BS <- as.numeric(tree$node.label)
  
  # create a vector to store colors
  p <- character(length(BS))
  
  # define colors for node labels
  co <- c("black", "red", "blue")
  
  # define what color goes with what support
  p[BS == 100] <- co[1]
  p[BS <= 99 & BS >= 95] <- co[2]
  p[BS < 95] <- co[3]
  
  # this is for the 'null' node of the tree,
  # later removed
  p[1] <- "green"
  
  # plot the tree and node labels
  plot(tree, type="u", lab4ut="axial", direction="l", label.offset=0.01, no.margin=T, cex=0.4)
  nodelabels(pch=21, bg=p, cex=0.5, lwd=0.1)
}

# add legend
legend("bottomleft", legend=expression(100 == BS, 95 <= BS * " <= 99", BS < 95), pch=21, pt.bg=co, bty="n", cex=0.5)

# close pdf device
dev.off()