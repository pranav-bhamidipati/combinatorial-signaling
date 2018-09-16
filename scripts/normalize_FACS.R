## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('lattice')
usePackage('heatmap3')

## Get list of raw data to be normalized
in_files <- list.files(
  'FACS/raw', 
  pattern = ".*\\.txt",
	full.names = T)

## Make list of output filenames
out_files <- sub(
	'-counts',
	'_FACS_norm',
	in_files,
	fixed = T
)
out_files <- sub(
	'/raw/',
	'/normalized/',
	out_files,
	fixed = T
)
	
## Read in table of # counts per cell
counts_per_cell <- read.table(
	'FACS/counts_per_cell.txt', 
	stringsAsFactors=F, 
	header=T)
	
## Read in raw data, normalize, write to file
for (i in 1:length(in_files)) {
  f = in_files[i]
	writeLines(paste0("Reading in \"",f,"\""))
	# Read
  dat <- read.table(
    file = f, 
    header = T, 
    stringsAsFactors = F
  )
	# Normalize
	dat$n_counts <- counts_per_cell$n_counts[match(dat$cell,counts_per_cell$cell)]
	dat$expr <- dat$expr / dat$n_counts * 10^6
  dat <- dat[-4]
	writeLines(paste0("Writing to \"",out_files[i],"\""))
	# Write
	write.table(
		x = dat,
		file = out_files[i],
		row.names = F,
		col.names = T
	)
}
	
	
	
	
	
	
	
	
	
	