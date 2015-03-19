library("ape")

# disable scientific notation
options(scipen=999)

# set working directory
setwd("~/Dropbox/Metazoan_Partitions/Metazoan_Project/RAxML_single_genes/")

# creates a list of all tree files in the working directory
infiles <- dir(pattern='*bipartitions.OG*')

extract.lengths <- function(file){
  tr <- read.tree(file)
  # identifies terminal branches with TRUE
  terms <- tr$edge[, 2] <= Ntip(tr)
  terminal.edges <- tr$edge.length[terms]
  # reordering tip labels and assigning them to terminal branches 
  names(terminal.edges) <- tr$tip.label[tr$edge[terms, 2]]
  # putting the data above in a table
  table <- data.frame(rbind(terminal.edges))
  #sorting columns by name
  sort_table <- table[,sort(names(table))]
  # get OG number from filename
  OG_no <- sub("RAxML_bipartitions.([A-Z]+[0-9]+).phy.no-undet.spurious-out-gappyout.fasta.aln.OUT", "\\1", perl=TRUE, x=file)
  # gets average length of tip branches
  avg_tip_length <- mean(terminal.edges)
  # selects tip branches that are 5 times or more longer than avgerage length
  sel <- terminal.edges >= 5 * avg_tip_length
  # add OG number to each selection
  long_br <- cbind(OG_no,terminal.edges[sel])
  # add OG number to each row of tip length table
  og_sort_table <- cbind(OG_no,sort_table)
  # put the long branched selections in a table
  long_br_table <- data.frame(long_br)
  # write table with tip lengths
  write.table(og_sort_table, file="OG_tip_lengths.csv", row.names=F, col.names=T, append=T, quote=F)
  # write table with selected OGs
  write.table(long_br_table, file="OG_long_branch.csv", row.names=T, col.names=F, append=T, quote=F)
}

# loop over all files
lapply(infiles, extract.lengths)