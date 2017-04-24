# disable scientific notation
options(scipen=999)

# install required libraries
# uncomment if you don't have these installed
#install.packages("ape")
#install.packages("seqinr")
#install.packages("data.table")

# load needed libraries
library("ape")
library("seqinr")
library("data.table")

# Set working directories
# and load all tree files with support.

# The relative paths in this code will work if you downloaded
# the Dryad folder and retained its structure.

# Change your working directory to where this script resides
# on your system or use:

#source("your/path_to/Code/gene_stats.R", chdir=T)

# for the relative paths to work correctly.

trees_dir <- file.path("..", "RAxML_single_genes/")
trees_bip <- dir(path=trees_dir, pattern="*bipartitions.OG*")

# load all alignment files in FASTA format
fasta_dir <- file.path("..", "Single_gene_alignments/")
alignments <- dir(path=fasta_dir, pattern="*fasta.aln")


### AVERAGE BOOTSTRAP SUPPORT ###

# This code should work with any Newick tree files
# that have some measure of support at nodes, 
# including PP from Bayesian analysis.

# define a function to calculate average support
Avg_support <- function(file) {
  
  # get OG number from filename
  OG_no <- sub("RAxML_bipartitions.(OG[0-9]+).OUT", "\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # store support values in a vector
  support <- c(as.numeric(tree$node.label))
  # calculate average support
  avg_supp <- mean(support, na.rm=T)
  return(c(OG_no,avg_supp))
  
}

# loop over all files
average_bootstrap <- lapply(trees_bip, Avg_support)
average_bootstrap <- data.frame(matrix(unlist(average_bootstrap), nrow=(length(average_bootstrap)), byrow=T))
colnames(average_bootstrap) <- c("OG", "Average_bootstrap")


### TOTAL AND AVERAGE BRANCH LENGTHS ###

# This takes a Newick tree with branch lengths
# and returns the number of tips for each tree,
# calculates total tree length, average branch length,
# and variance of branch lengths across the tree.

Br_length.trees <- function(file) {
  
  # reads the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # gets number of tips
  no_tips <- length(tree$tip.label)
  # get tree length
  tr_length <- sum(tree$edge.length)
  # calculate avg branch length
  avg_br_length <- mean(tree$edge.length)
  # calculate variance in branch lengths
  var_br_length <- var(tree$edge.length)
  # get OG number from filename
  OG_no <- sub("RAxML_bipartitions.(OG[0-9]+).OUT", "\\1", perl=TRUE, x=file)
  return(c(OG_no,no_tips,tr_length,avg_br_length,var_br_length))
  
}

# loop over all files
br_lengths <- lapply(trees_bip, Br_length.trees)
br_lengths <- data.frame(matrix(unlist(br_lengths), nrow=(length(br_lengths)), byrow=T))
colnames(br_lengths) <- c("OG", "No_of_taxa", "Tree_length", "Average_branch_length", "Variance_in_branch_lengths")

### LONG-BRANCH SCORES ###

# This code calculates long-branch scores 
# as defined by Struck 2014 Evol. Bioinform. 10:51-67
# for each taxon in each locus, as well as standard deviation
# of these scores in each locus

LB_score <- function(file) {
  
  # read the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # make matrix of pairwise distances in branch lengths from the tree
  cophentr <- cophenetic(tree)
  # get OG number from filename
  OG_no <- sub("RAxML_bipartitions.(OG[0-9]+).OUT", "\\1", perl=TRUE, x=file)
  # create empty data frame
  LB_table <- list()
  
  # loop over all taxa in the matrix
  Scoring <- function(taxon) {
    
    # calculate mean of patristic distances for the taxon
    tax_mean <- mean(cophentr[,taxon])
    # calculate Struck's LB score
    LB_score <- ( tax_mean / mean(cophentr) - 1 ) * 100
    # define a row for the table
    return <- c(taxon, tax_mean, LB_score, OG_no)
    
  }
  
  LB_table <- lapply(row.names(cophentr), Scoring)
  LB_table <- data.frame(matrix(unlist(LB_table), nrow=(length(LB_table)), byrow=T))
  colnames(LB_table) <- c("Taxon", "Tax_mean", "LB_score", "OG")

  # calculate standard deviation of LB
  LB_SD <- sd(LB_table$LB_score, na.rm=T)
  return(data.frame(LB_SD, LB_table))
  
}

# loop over all files
LB_scores <- lapply(trees_bip, LB_score)
# compress the list to a data frame
LB_scores <- rbindlist(LB_scores)


### PLOTTING GENE TREES ### 

# This will save png files of 600 x 600 pixel plots
# of all input phylogenetic trees with tip labels
# and support values
# the plots are created in the same folder as tree files

Plot_trees <- function(file) {
  
  # reads the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # extracts plot name (OG) from file name 
  plot.name <- sub("RAxML_bipartitions.(OG[0-9]+).OUT", "\\1", perl=TRUE, x=file)
  # open png file
  png(file=paste(trees_dir, plot.name, "-tree.png", sep=""), width=600, height=600)
  # get sum of branch lengths
  br_length <- sum(tree$edge.length)
  # plot tree with BS support values
  plot.phylo(tree, show.node.label=T)
  # give title as OG number and subtitle as tree length
  title(main=plot.name,sub=paste("Tree length: ", br_length, sep=""))
  # close png file
  dev.off()
  
}

# loop over all files
lapply(trees_bip, Plot_trees)


### SATURATION ###

# The code below takes as input FASTA alignments and tree files with branch lengths
# and calculates regression slope and r-squared for each locus
# also saving regression plots for each locus 
# in a newly created 'Saturation_plots' folder

# get current working directory
work_dir <- getwd()
# one level up
up_work_dir <- sub("/Code", "", work_dir)

# create a folder for saturation plots
dir.create("../Saturation_plots")
sat_dir <- file.path(up_work_dir, "Saturation_plots/")

# get new absolute names of directories
al_dir <- file.path(up_work_dir, "Single_gene_alignments/")
tr_dir <- file.path(up_work_dir, "RAxML_single_genes/")

# define function to perform regression, 
# calculate R-squared, and save saturation plots 
Saturation <- function(seq, tree) {
  
  # read in alignment
  alignment <- read.alignment(file=seq, format="fasta")
  # read in tree
  tree <- read.tree(tree)
  # get OG number from filename
  OG_no <- sub(".*/(OG[0-9]+).fasta.aln", "\\1", perl=TRUE, x=seq)
  # matrix with pairwise identity
  mat <- dist.alignment(alignment, matrix="identity")
  # matrix with uncorrected p-distances
  p_mat <- mat*mat
  # make matrix of pairwise distances in branch lengths from the tree
  cophentr <- cophenetic(tree)  
  # store as matrix objects
  mat_mat <- as.matrix(mat)
  mat_p_mat <- as.matrix(p_mat)
  # order p-distance matrix by names
  mat_p_mat <- mat_p_mat[order(row.names(mat_p_mat)),order(row.names(mat_p_mat))]
  mat_co <- as.matrix(cophentr)
  # order pairwise distances matrix by names
  mat_co <- mat_co[order(row.names(mat_co)),order(row.names(mat_co))]
  # get lower triangulars of both matrices
  branch_dist <- mat_co[lower.tri(mat_co)]
  p_dist <- mat_p_mat[lower.tri(mat_p_mat)]
  # perform simple linear regression
  regress <- lm(p_dist ~ branch_dist)
  # get slope
  slope <- coef(regress)[2]
  # get r-squared
  Rsquared <- summary(regress)$r.squared
  
  # plot branch length pairwise distances on x
  # and uncorrected p-distances on y
  
  # open png file
  png(file=paste(sat_dir, OG_no, "-saturation.png", sep=""), width=600, height=600)
  
  plot(branch_dist, p_dist)
  # add simple linear regression line
  abline(lm(p_dist ~ branch_dist), col="red")
  # give title as OG number and subtitle as tree length
  title(main=OG_no,sub=paste("Slope: ", round(slope, digits=3), " R-squared: ", round(Rsquared, digits=3), sep=""), cex.sub=1.25)
  # close png file
  dev.off()
  
  return(list(OG_no, slope, Rsquared))
  
}

# create a table with file names
files_table <- as.data.frame(cbind(paste(al_dir, alignments, sep=""), paste(tr_dir, trees_bip, sep="")))

saturation_table <- t(mapply(Saturation, as.matrix(files_table$V1), as.matrix(files_table$V2)))
saturation_table <- as.data.frame(saturation_table)
colnames(saturation_table) <- c("OG","Slope","R-squared")
row.names(saturation_table) <- NULL


### ALIGNMENT LENGTH AND MISSING DATA ###

# This code assumes that the format is FASTA and there are 
# 36 taxa in concatenated matrix, but with appropriate adjustments 
# it should handle any of the formats supported by seqinr.
# This code also assumes all gaps are coded as "-".

Missing_data <- function(file) {
  
  # read in alignment in FASTA format
  alignment <- read.alignment(file=paste(fasta_dir, file, sep=""), format="fasta")
  # get OG number from filename
  OG_no <- sub("(OG[0-9]+).fasta.aln", "\\1", perl=TRUE, x=file)
  # alignment length, assuming it is the same across taxa
  al_length <- nchar(alignment$seq[1])
  # count the number of cells in the alignment
  al_total_cells <- sum(sapply(alignment$seq, nchar))
  # count gaps, or any "-" in the alignment
  gaps <- unname((sapply(regmatches(alignment, gregexpr("-", alignment)), length))[3])
  # this part scales missing data to 36 taxa
  al_cells_36tax <- 36 * al_length
  gaps_to_36tax <- ( 36 - alignment$nb ) * al_length
  # count the proportion of missing
  missing <- gaps / al_total_cells 
  # count the proportion of missing scaled to 36 taxa
  missing_scaled <- ( gaps + gaps_to_36tax ) / al_cells_36tax
  
  return(list(OG_no,al_length,missing,missing_scaled))
  
}

# apply to all files
missing <- lapply(alignments, Missing_data)
missing <- rbindlist(missing)
setnames(missing, names(missing), c("OG","Length","Missing","Missing_scaled"))

# putting together all the data except LB scores, 
# which are taxon- and gene-specific
dfs <- list(average_bootstrap, br_lengths, missing, saturation_table)

Multmerge <- function(dfs){
  datalist <- lapply(dfs, function(x){data.frame(x)})
  Reduce(function(x,y) {merge(x,y)}, datalist)
}

All_loci_stats <- Multmerge(dfs)
