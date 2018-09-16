## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('lattice')
usePackage('heatmap3')

## Read in master gene list
genes <- read.table(file = "droplet/genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

## Get list of data sets
my_filenames <- list.files(
  'FACS/FACS', 
  pattern = ".*\\.csv")

## Convert each dense data matrix into a melted sparse matrix
## and save in raw folder
for (i in 1:length(my_filenames)) {
	
	# Get data from file
	my_f = my_filenames[i]
	print(paste("Converting", my_f))
	dat <- read.csv(
		paste0('FACS/FACS/',my_f), 
    header = T, 
    stringsAsFactors = F
	)
	
	# Melt into long form and reformat
	dat = melt(dat)
	colnames(dat) <- c('gene', 'cell','expr')
	dat$gene <- match(tolower(dat$gene), genes)
	dat$cell <- as.character(dat$cell)	
	
	# Remove zero-count rows
	dat <- dat[dat$expr > 0,]
	
	# Write to file
	write.table(
		x = dat,
		file = paste0("FACS/raw/", 
			substr(my_f,1,nchar(my_f) - 4), 
			'.txt'),
		row.names = FALSE
	)
}









