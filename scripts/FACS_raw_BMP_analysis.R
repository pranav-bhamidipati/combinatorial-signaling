## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('lattice')
usePackage('heatmap3')

## Read in master gene list and signaling gene lists
genes <- read.table(file = "genes.txt", stringsAsFactors = F)
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
	
## Read in FACS annotations and metadata
FACS_annots <- read.csv(
	'FACS/annotations_FACS.csv', 
	stringsAsFactors=F, 
	header=T)

FACS_meta <- read.csv(
	'FACS/metadata_FACS.csv', 
	stringsAsFactors=F, 
	header=T)

## Define the list of data sets
my_filenames <- list.files(
  'FACS/raw', 
  pattern = ".*\\.txt")

## Define the total set of genes to be studied
all_genes <- BMP_genes_annotated$idx[
  grepl(pattern = "ligand", x = BMP_genes_annotated$role) | 
    grepl(pattern = "receptor", x = BMP_genes_annotated$role)
  ]

## Initialize the master data list and counts_per_cell list
master_dat <- as.list(my_filenames)
names(master_dat) <- my_filenames

counts_per_cell <- as.list(my_filenames)
names(counts_per_cell) <- my_filenames

## Extract required data
for (i in 1:length(my_filenames)) {
  # Select data set for this iteration
  my_f <- my_filenames[i]
  
  # Read in data
  dat <- read.table(
    paste0('FACS/raw/',my_f), 
    header = T, 
    stringsAsFactors = F
  )
	
	# Keep track of # counts per cell
	counts_per_cell[[my_f]] <- aggregate(
		dat$expr,
		by=list(dat$cell),
		FUN=sum
	)
	
  # Select data for the genes we want
  dat <- dat[dat$gene %in% all_genes,]
  dat$tissue <- sub("-counts.txt","",my_f)
  
  # Add the data to the master data set list
  master_dat[[my_f]] <- dat
}

# Knit master data set list into one table
master_dat <- do.call("rbind.data.frame",master_dat)

# Knit counts_per_cell into one table
counts_per_cell <- do.call("rbind.data.frame",counts_per_cell)
colnames(counts_per_cell) <- c('cell','n_counts')
counts_per_cell$tissue <- sapply(
	strsplit(
		rownames(counts_per_cell), 
		'-counts.txt'
	), 
	function(x) x[1]
)
rownames(counts_per_cell) <- 1:nrow(counts_per_cell)

## Write master_dat and counts_per_cell to file
write.table(
  x = master_dat, 
  file = "FACS/all_tissues_BMP_ligs_recs/BMP_lig_rec_dat.txt", 
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  x = counts_per_cell, 
  file = "FACS/counts_per_cell.txt", 
  row.names = FALSE,
  col.names = TRUE
)

## Read master_dat and counts_per_cell from file
master_dat <- read.table(
	file = "FACS/all_tissues_BMP_ligs_recs/BMP_lig_rec_dat.txt", 
	header = T,
	stringsAsFactors= F
)

counts_per_cell <- read.table(
	file = "FACS/counts_per_cell.txt", 
	header = T,
	stringsAsFactors= F
)


## After adding all the needed data...

## Define list of tissues
my_tissues = unlist(strsplit(my_filenames, '-counts.txt'))

## Define list of cells
all_cells <- sort(unique(master_dat$cell))

## Calculate total # cells per sample
cells_per_sample <- data.frame(table(counts_per_cell$tissue))
colnames(cells_per_sample) <- c('sample','n_cells')

## Restrict data to specific genes of interest	
my_genes = c(BMP_ligands, BMP_receptors)
my_filter = 0
my_dat <- master_dat[master_dat$expr > my_filter,]
my_dat <- my_dat[my_dat$gene %in% my_genes,]



### ***Plots***

## Plot histograms of expression for each BMP receptor	
pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_rec_expr_raw_1.pdf")
par(mfrow=c(3, 3))
for (rec in BMP_receptors[1:9]) {
  x = my_dat$expr[my_dat$gene == rec]
  x = c(x, numeric(sum(cells_per_sample$n_cells) - length(x)))
  x = log10(x + 1)
	h = hist(
		x = x, 
		breaks=seq(-0.2, 6, 0.2), 
		plot = FALSE
	)
	h$counts <- log10(h$counts + 1)
	plot(
		h,
    col="darkslategray", 
    main=paste0(
      genes[rec], 
      " (n = ",
      sum(cells_per_sample$n_cells),
      ")"
    ),
		xlim=c(0, 6), 
    xlab = "Log10(counts + 1)",
    ylab = "Log10(# cells + 1)"
  )
}
dev.off()

pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_rec_expr_raw_2.pdf")
par(mfrow=c(3, 3))
for (rec in BMP_receptors[10:11]) {
  x = my_dat$expr[my_dat$gene == rec]
  x = c(x, numeric(sum(cells_per_sample$n_cells) - length(x)))
  x = log10(x + 1)
	h = hist(
		x = x, 
		breaks=seq(-0.2, 6, 0.2), 
		plot = FALSE
	)
	h$counts <- log10(h$counts + 1)
	plot(
		h,
    col="darkslategray", 
    main=paste0(
      genes[rec], 
      " (n = ",
      sum(cells_per_sample$n_cells),
      ")"
    ),
		xlim=c(0, 6), 
    xlab = "Log10(counts + 1)",
    ylab = "Log10(# cells + 1)"
  )
}
dev.off()



## Histograms for each ligand
pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_raw_1.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[1:9]) {
  x = my_dat$expr[my_dat$gene == lig]
  x = c(x, numeric(sum(cells_per_sample$n_cells) - length(x)))
  x = log10(x + 1)
  y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
  plot(
    x = seq(-0.25, 6.25, 0.5),
    y = log10(y$counts + 1),  
    type='h', 
    lwd=10, 
    lend=2,
    main=paste0(
      genes[lig], 
      " (n = ",
      sum(cells_per_sample$n_cells),
      ")"
    ),
    xlab = "Log10(expr + 1)",
    ylab = "Log10(# cells + 1)", 
    col="darkslategray", 
  )
}
dev.off()

pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_raw_2.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[10:13]) {
  x = my_dat$expr[my_dat$gene == lig]
  x = c(x, numeric(sum(cells_per_sample$n_cells) - length(x)))
  x = log10(x + 1)
  y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
  plot(
    x = seq(-0.25, 6.25, 0.5),
    y = log10(y$counts + 1),  
    type='h', 
    lwd=10, 
    lend=2,
    main=paste0(
      genes[lig], 
      " (n = ",
      sum(cells_per_sample$n_cells),
      ")"
    ),
    xlab = "Log10(expr + 1)",
    ylab = "Log10(# cells + 1)", 
    col="darkslategray", 
  )
}
dev.off()





