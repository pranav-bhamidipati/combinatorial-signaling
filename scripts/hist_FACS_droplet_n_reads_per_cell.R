## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('lattice')
usePackage('heatmap3')
usePackage('Matrix')

## Read in data to be plotted
FACS_counts_per_cell <- read.table(
	file = "FACS/counts_per_cell.txt",
	stringsAsFactors = FALSE,
	header = TRUE
)
droplet_counts_per_cell <- read.table(
	file = "droplet/counts_per_cell.txt",
	stringsAsFactors = FALSE,
	header = TRUE
)



### ***		Plots		***


## Plot # counts per cell for FACS data
pdf("FACS/hist_n_counts_per_cell.pdf")
x = FACS_counts_per_cell$n_counts
x = log10(x)
h = hist(
	x = x, 
	breaks=seq(0, 7.4, 0.2), 
	plot = FALSE
)
h$counts <- log10(h$counts + 1)
plot(
	h,
	col="darkslategray", 
	main=paste0(
		"FACS data: # counts per cell", 
		" (n = ",
		nrow(FACS_counts_per_cell),
		")"
	),
	xlim=c(0, 9), 
	xlab = "Log10(counts per cell)",
	ylab = "Frequency (Log10(# cells))",
)
dev.off()

## Plot # counts per cell for droplet data
pdf("droplet/hist_n_counts_per_cell.pdf")
x = droplet_counts_per_cell$n_counts
x = log10(x)
h = hist(
	x = x, 
	breaks=seq(0, 7.4, 0.2), 
	plot = FALSE
)
h$counts <- log10(h$counts + 1)
plot(
	h,
	col="darkslategray", 
	main=paste0(
		"Droplet data: # counts per cell", 
		" (n = ",
		nrow(FACS_counts_per_cell),
		")"
	),
	xlim=c(0, 9), 
	xlab = "Log10(counts per cell)",
	ylab = "Frequency (Log10(# cells))",
)
dev.off()













