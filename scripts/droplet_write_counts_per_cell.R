## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('lattice')
usePackage('heatmap3')
usePackage('Matrix')

## Read in the list of data sets
my_filenames <- list.files(
  'droplet/raw', 
  pattern = ".*\\.mtx"
)

## Initialize list of counts per cell
counts_per_cell <- as.list(my_filenames)
names(counts_per_cell) <- my_filenames

## Extract required data
for (i in 1:length(my_filenames)) {
  # Select data set for this iteration
  my_f <- my_filenames[i]
  
  # Read in data
	dat <- readMM(paste0('droplet/raw/',my_f))
	dat <- as.data.frame(summary(dat))
	
	# Keep track of # counts per cell
	counts_per_cell[[my_f]] <- aggregate(
		dat$x,
		by=list(dat$j),
		FUN=sum
	)
}

# Knit counts_per_cell into one table
for (i in 1:length(counts_per_cell)) {
	 counts_per_cell[[i]]$Group.1 <- paste0(
		substr(
			names(counts_per_cell)[i],
			1,
			nchar(names(counts_per_cell)[i]) - 3
		),
		counts_per_cell[[i]]$Group.1
	)
}
counts_per_cell <- do.call("rbind.data.frame",counts_per_cell)
colnames(counts_per_cell) <- c('cell','n_counts')
counts_per_cell$tissue <- sapply(
	strsplit(
		rownames(counts_per_cell), 
		'-',
		fixed = TRUE
	), 
	function(x) x[1]
)
rownames(counts_per_cell) <- 1:nrow(counts_per_cell)

## Write counts_per_cell to file
write.table(
  x = counts_per_cell, 
  file = "droplet/counts_per_cell.txt", 
  row.names = FALSE,
  col.names = TRUE
)












