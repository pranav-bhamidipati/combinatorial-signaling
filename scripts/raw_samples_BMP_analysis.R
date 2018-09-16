## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('lattice')
usePackage('reshape2')
usePackage('heatmap3')

## Read in master gene list and signaling gene list
genes <- read.table(file = "droplet/genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_genes <- read.table('./BMP_genes.txt', stringsAsFactors = F)
BMP_genes <- tolower(BMP_genes$V1)
BMP_genes_idx <- which(genes %in% BMP_genes)
BMP_genes_annotated <- read.table('./BMP_genes_annotated.txt', 
                                  stringsAsFactors = F,
                                  header = T)
all_BMP_receptors <- BMP_genes_annotated$idx[
  grepl('receptor',BMP_genes_annotated$role)
  ]
all_BMP_ligands <- BMP_genes_annotated$idx[
  grepl('ligand',BMP_genes_annotated$role)
  ]

BMP_receptors <- read.table('./BMP_Receptors.txt', stringsAsFactors = F)
BMP_receptors <- which(genes %in% tolower(unlist(BMP_receptors)))
BMP_ligands <- read.table('./BMP_Ligands.txt', stringsAsFactors = F)
BMP_ligands <- which(genes %in% tolower(unlist(BMP_ligands)))


## Read in the list of data sets
dat_filenames <- list.files(
  'droplet/', 
  pattern = ".*\\.mtx")

## Define the total set of genes to be studied
my_genes <- BMP_genes_annotated$idx[
  grepl(pattern = "ligand", x = BMP_genes_annotated$role) | 
    grepl(pattern = "receptor", x = BMP_genes_annotated$role)
  ]

## Initialize the master data set and cell count vector
master_dat <- data.frame(
  gene=integer(0), 
  cell=integer(0), 
  expr=numeric(0),
  tissue=character(0),
  tissue.cell=character(0)
)

## Extract required data
for (i in 1:length(dat_filenames)) {
  # Select data set for this iteration
  my_dat <- dat_filenames[i]
  
  # Read in data
	dat <- readAllMtx('droplet/', my_dat, fixed = T)[[1]]
  dat <- as.data.frame(summary(dat))
	colnames(dat) <- c('gene', 'cell','expr')
  dat$tissue <- rep(my_dat, nrow(dat))
  dat$tissue.cell <- paste0(dat$tissue,'.',dat$cell)

  # Select data for the genes we want
  dat <- dat[dat$gene %in% my_genes,]
  
  # Add the data to the master data set
  master_dat <- rbind.data.frame(master_dat, dat)
}

## Calculate # counts per gene
n_counts_per_gene <- aggregate(
	master_dat$expr, 
	by=list(master_dat$gene), 
	FUN=sum)
colnames(n_counts_per_gene) = c('gene', 'n_counts')
n_counts_per_gene$gene = genes[n_counts_per_gene$gene]

## Plot frequency distribution of # counts per cell for each gene 
# BMP receptors
pdf('droplet/all_tissues_BMP_ligs_recs/n_reads_per_cell_per_BMP_rec.pdf')
par(mfrow=c(3, 3))
for (g in BMP_receptors) {
	freqs = table(master_dat$expr[master_dat$gene == g])
  barplot(
    freqs, 
    main = paste0(
      genes[g],
      " (n = ", 
      sum(master_dat$gene == g),
      ")"
    ),
    xlab = "# counts",
    ylab = "Frequency(# cells)",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

# BMP ligands
pdf('droplet/all_tissues_BMP_ligs_recs/n_reads_per_cell_per_BMP_lig.pdf')
par(mfrow=c(3, 4))
for (g in BMP_ligands) {
	if (!(g %in% master_dat$gene)) next
	freqs = table(master_dat$expr[master_dat$gene == g])
  barplot(
    freqs, 
    main = paste0(
      genes[g],
      " (n = ", 
      sum(master_dat$gene == g),
      ")"
    ),
    xlab = "# counts",
    ylab = "Frequency(# cells)",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

## Calculate # reads per cell
a <- aggregate(
	master_dat$x,
	by=list(master_dat$j),
	FUN = sum
)
n_reads_per_cell <- unlist(sapply(n_reads_per_cell, function(x) x[2]))

## Plot histogram of # reads per cell
y = hist(n_reads_per_cell, 
	breaks=seq(0, 180000, 2000), 
	plot=F)
pdf('droplet/n_reads_per_cell_pooled.pdf')
hist(
	x = log10(n_reads_per_cell), 
	breaks=100,
	main=paste0(
		"Log10(# counts) per cell\nAll samples pooled (n = ",
		length(n_reads_per_cell),
		")"	
	),
	xlab = "Log10(# counts)",
	ylab = "Frequency (# cells)",
	col="darkslategray"
)
dev.off()

## Get quartiles of # reads per cell for each BMP gene
n_reads_per_BMP_gene <- aggregate(
	dat$expr,
	by = list(dat$gene),
	FUN = quantile
)
n_reads_per_BMP_gene <- cbind(n_reads_per_BMP_gene[1], n_reads_per_BMP_gene[[2]])
colnames(n_reads_per_BMP_gene)[1] <- 'gene'
n_reads_per_BMP_gene$gene = genes[n_reads_per_BMP_gene$gene]











