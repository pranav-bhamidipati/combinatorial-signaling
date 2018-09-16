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
	
## Read in FACS annotations, metadata, and # raw counts per cell 
FACS_annots <- read.csv(
	'FACS/annotations_FACS.csv', 
	stringsAsFactors=F, 
	header=T)

FACS_meta <- read.csv(
	'FACS/metadata_FACS.csv', 
	stringsAsFactors=F, 
	header=T)

counts_per_cell <- read.table(
	'FACS/counts_per_cell.txt', 
	stringsAsFactors=F, 
	header=T)

## Read in the list of data sets
my_filenames <- list.files(
  'FACS/normalized', 
  pattern = ".*\\.txt")

## Define the total set of genes to be studied
all_genes <- BMP_genes_annotated$idx[
  grepl(pattern = "ligand", x = BMP_genes_annotated$role) | 
    grepl(pattern = "receptor", x = BMP_genes_annotated$role)
  ]

## Initialize the master data list and counts_per_cell list
master_dat <- as.list(my_filenames)
names(master_dat) <- my_filenames

## Extract required data
for (i in 1:length(my_filenames)) {
  # Select data set for this iteration
  my_f <- my_filenames[i]
  
  # Read in data
  dat <- read.table(
    paste0('FACS/normalized/',my_f), 
    header = T, 
    stringsAsFactors = F
  )
	
  # Select data for the genes we want
  dat <- dat[dat$gene %in% all_genes,]
  dat$tissue <- sub("_FACS_norm.txt","",my_f)
	
  # Add the data to the master data set list
  master_dat[[my_f]] <- dat
}

# Knit master data set list into one table
master_dat <- do.call("rbind.data.frame",master_dat)

## Write master_dat to file
write.table(
  x = master_dat, 
  file = "FACS/all_tissues_BMP_ligs_recs/BMP_lig_rec_dat.txt", 
  row.names = FALSE,
  col.names = TRUE
)

## If needed, read in this data
master_dat <- read.table(
	"FACS/all_tissues_BMP_ligs_recs/BMP_lig_rec_dat.txt", 
	stringsAsFactors=F, 
	header=T)

counts_per_cell <- read.table(
	"FACS/counts_per_cell.txt", 
	stringsAsFactors=F, 
	header=T)

	
	
	
## After adding all the needed data...

## Define list of tissues
my_tissues = unlist(strsplit(my_filenames, '_FACS_norm.txt'))

## Define list of cells
all_cells <- sort(unique(master_dat$cell))

## Calculate total # cells per sample
cells_per_sample <- as.data.frame(
	table(counts_per_cell$tissue),
	stringsAsFactors = FALSE
)
colnames(cells_per_sample) <- c('sample','n_cells')

## Restrict data to specific genes of interest
my_genes = c(BMP_ligands, BMP_receptors)
my_filter = 10
my_dat <- master_dat[master_dat$expr > my_filter,]
my_dat <- my_dat[my_dat$gene %in% my_genes,]


## Compute metrics on genes of interest
# Descrptive statistics
#Pooled tissues
BMP_stats <- aggregate(
  my_dat$expr,
  by = list(my_dat$gene),
  FUN = function(x) {
    c(
			length(x),        #Determines # of cells ON after thresholding
      min(x),
      mean(x),
      median(x),
      max(x)
    )
  }
)
BMP_stats <- cbind (BMP_stats[1],as.data.frame(BMP_stats[,2]))
BMP_stats$Group.1 <- genes[BMP_stats$Group.1]
colnames(BMP_stats) <- c(
	'gene', 
  'n_cells_ON',
  'min',
  'mean',
  'median',
  'max')
BMP_stats$role <- BMP_genes_annotated$role[
	match(
      BMP_stats$gene,
      BMP_genes_annotated$gene
  )
]
BMP_stats <- BMP_stats[c(
	1,
	ncol(BMP_stats),
	2:(ncol(BMP_stats) - 1)
	)]
BMP_stats <- BMP_stats[
  na.omit(
    match(
      BMP_genes_annotated$gene,
      BMP_stats$gene
    )
  ),
  ]
rownames(BMP_stats) <- 1:nrow(BMP_stats)

# Save to file
write.table(
  x = BMP_stats[order(BMP_stats$role, BMP_stats$gene),], 
  file = "FACS/all_tissues_BMP_ligs_recs/BMP_stats.txt", 
  row.names = FALSE,
  col.names = TRUE
)

# Read from file
BMP_stats <- read.table(
	file = "FACS/all_tissues_BMP_ligs_recs/BMP_stats.txt",
	header = TRUE,
	stringsAsFactors = FALSE
)


#Per tissue	
BMP_tissue_stats <- lapply(
	X = unique(my_tissues),
	FUN = function(tissue) {
		tissue_rows = grepl(tissue, my_dat$tissue)
		my_stats = aggregate(
			my_dat$expr[tissue_rows],
			by = list(my_dat$gene[tissue_rows]),
			FUN = function(x) {
				c(
					length(x),        #Determines # of cells ON **after thresholding**
					min(x),
					mean(x),
					median(x),
					max(x)
				)
			}
		)
		my_stats <- cbind (my_stats[1],as.data.frame(my_stats[,2]))
		my_stats$Group.1 <- genes[my_stats$Group.1]
		colnames(my_stats) <- c('gene', 
														 'n_cells_ON',
														 'min',
														 'mean',
														 'median',
														 'max')
		my_stats$role <- BMP_genes_annotated$role[
			match(
					my_stats$gene,
					BMP_genes_annotated$gene
			)
		]
		my_stats <- my_stats[c(
			1,
			ncol(my_stats),
			2:(ncol(my_stats) - 1)
			)]
		rownames(my_stats) <- 1:nrow(my_stats)
		my_stats <- my_stats[
			na.omit(
				match(
					BMP_genes_annotated$gene,
					my_stats$gene
				)
			),
		]
	}
)
names(BMP_tissue_stats) <- unique(my_tissues)

## Compute P/A matrices in each cell for genes of interest 
# BMP receptors
dat <- my_dat
BMP_rec_pa <- sapply(
	X = BMP_receptors,
	FUN = function(g) {
		g.cells = dat$cell[dat$gene == g]
		as.integer(all_cells %in% g.cells)
	}
)
BMP_rec_pa <- as.data.frame(BMP_rec_pa)
colnames(BMP_rec_pa) <- genes[BMP_receptors]
rownames(BMP_rec_pa) <- all_cells

# BMP ligands
dat <- my_dat,
BMP_lig_pa <- sapply(
	X = BMP_ligands,
	FUN = function(g) {
		g.cells = dat$cell[dat$gene == g]
		as.integer(all_cells %in% g.cells)
	}
)
BMP_lig_pa <- as.data.frame(BMP_lig_pa)
colnames(BMP_lig_pa) <- genes[BMP_ligands]
rownames(BMP_lig_pa) <- all_cells

## BMP words for each cell
BMP_rec_words <- do.call(paste0, BMP_rec_pa)
BMP_lig_words <- do.call(paste0, BMP_lig_pa)






### ***Plots***

## Plot expression distributions
## Pooled samples
# BMP receptors
pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_rec_expr.pdf")
par(mfrow=c(3, 3))
for (rec in BMP_receptors) {
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
		xlab = "Log10(expr + 1)",
		ylab = "Log10(# cells + 1)",
	)
}
dev.off()

# BMP ligands
pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_1.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[1:9]) {
  x = my_dat$expr[my_dat$gene == lig]
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
		xlab = "Log10(expr + 1)",
		ylab = "Log10(# cells + 1)",
	)
}
dev.off()

pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_2.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[10:13]) {
  x = my_dat$expr[my_dat$gene == lig]
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
		xlab = "Log10(expr + 1)",
		ylab = "Log10(# cells + 1)",
	)
}
dev.off()

## Each sample
# BMP receptors and ligands
## Plot expression histograms for each rec/lig in each tissue
for (tissue in unique(my_tissues)) {
  
	writeLines(tissue)
  tissue_rows = grepl(tissue, my_dat$tissue)
  tissue_samples = grepl(tissue, cells_per_sample$sample)
  tissue_path = paste0(
    "FACS",
		"/each_tissue_BMP_ligs_recs/",
		tissue,
		'/'
	)
	
  # BMP Receptors
  pdf(
		paste0(
			tissue_path, 
			tissue,
			'_BMP_rec_expr', 
			".pdf"
		)
	)
  par(mfrow=c(3, 3))
  for (g in BMP_receptors) {
    x = my_dat$expr[my_dat$gene == g & tissue_rows]
    x = c(
			x, 
			numeric(sum(cells_per_sample$n_cells[tissue_samples]) - length(x))
    )
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
				tissue,
        ': ',
        genes[g], 
				" (n = ",
				sum(cells_per_sample$n_cells[tissue_samples]),
				")"
			),
			xlim=c(0, 6), 
			xlab = "Log10(expr + 1)",
			ylab = "Log10(# cells + 1)",
		)
  }
  dev.off()
  
  # BMP Ligands
  # 1
  pdf(
		paste0(
			tissue_path, 
			tissue,
			'_BMP_lig_expr_1', 
			".pdf"
		)
	)
  par(mfrow=c(3, 3))
  for (g in BMP_ligands[1:9]) {
    x = my_dat$expr[my_dat$gene == g & tissue_rows]
    x = c(
			x, 
			numeric(sum(cells_per_sample$n_cells[tissue_samples]) - length(x))
    )
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
				tissue,
        ': ',
        genes[g], 
				" (n = ",
				sum(cells_per_sample$n_cells[tissue_samples]),
				")"
			),
			xlim=c(0, 6), 
			xlab = "Log10(expr + 1)",
			ylab = "Log10(# cells + 1)",
		)
  }
  dev.off()
  
	# BMP Ligands
  # 2
  pdf(
		paste0(
			tissue_path, 
			tissue,
			'_BMP_lig_expr_2', 
			".pdf"
		)
	)
  par(mfrow=c(3, 3))
  for (g in BMP_ligands[10:13]) {
    x = my_dat$expr[my_dat$gene == g & tissue_rows]
    x = c(
			x, 
			numeric(sum(cells_per_sample$n_cells[tissue_samples]) - length(x))
    )
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
				tissue,
        ': ',
        genes[g], 
				" (n = ",
				sum(cells_per_sample$n_cells[tissue_samples]),
				")"
			),
			xlim=c(0, 6), 
			xlab = "Log10(expr + 1)",
			ylab = "Log10(# cells + 1)",
		)
  }
  dev.off()

}


## Plot distribution of # genes ON per cell
## Each tissue separate
# BMP receptors
# 1
pdf("FACS/each_tissue_BMP_ligs_recs/each_tissue_n_BMP_receptors_1.pdf")
par(mfrow=c(3, 3))
for (tissue in unique(my_tissues)[1:9]) {
  dat_n_cells = sum(
		cells_per_sample$n_cells[
			grep(tissue, cells_per_sample$sample)
		]
	)
  dat <- my_dat[
    grep(
			pattern = tissue,
      x = my_dat$tissue,
      fixed = TRUE
		),
  ]
	genes_ON_per_cell <- aggregate(
    dat$gene, 
    by=list(dat$cell), 
    function(x) sum(x %in% BMP_receptors)
  )
  freqs = table(genes_ON_per_cell$x)
	freqs[1] = dat_n_cells - sum(freqs[2:length(freqs)])	
  
	barplot(
    freqs, 
    main = paste0(
      tissue,
      " (n = ", 
      dat_n_cells,
      ")"
    ),
    xlab = "# BMP receptors ON",
    ylab = "# cells",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

# 2
pdf("FACS/each_tissue_BMP_ligs_recs/each_tissue_n_BMP_receptors_2.pdf")
par(mfrow=c(3, 3))
for (tissue in unique(my_tissues)[10:18]) {
  dat_n_cells = sum(
		cells_per_sample$n_cells[
			grep(tissue, cells_per_sample$sample)
		]
	)
  dat <- my_dat[
    grep(
			pattern = tissue,
      x = my_dat$tissue,
      fixed = TRUE
		),
  ]
	genes_ON_per_cell <- aggregate(
    dat$gene, 
    by=list(dat$cell), 
    function(x) sum(x %in% BMP_receptors)
  )
  freqs = table(genes_ON_per_cell$x)
	freqs[1] = dat_n_cells - sum(freqs[2:length(freqs)])	
  
	barplot(
    freqs, 
    main = paste0(
      tissue,
      " (n = ", 
      dat_n_cells,
      ")"
    ),
    xlab = "# BMP receptors ON",
    ylab = "# cells",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

# BMP ligands
# 1
pdf("FACS/each_tissue_BMP_ligs_recs/each_tissue_n_BMP_ligands_1.pdf")
par(mfrow=c(3, 3))
for (tissue in unique(my_tissues)[1:9]) {
  dat_n_cells = sum(
		cells_per_sample$n_cells[
			grep(tissue, cells_per_sample$sample)
		]
	)
  dat <- my_dat[
    grep(
			pattern = tissue,
      x = my_dat$tissue,
      fixed = TRUE
		),
  ]
	genes_ON_per_cell <- aggregate(
    dat$gene, 
    by=list(dat$cell), 
    function(x) sum(x %in% BMP_ligands)
  )
  freqs = table(genes_ON_per_cell$x)
	freqs[1] = dat_n_cells - sum(freqs[2:length(freqs)])	
  
	barplot(
    freqs, 
    main = paste0(
      tissue,
      " (n = ", 
      dat_n_cells,
      ")"
    ),
    xlab = "# BMP ligands ON",
    ylab = "# cells",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

# 2
pdf("FACS/each_tissue_BMP_ligs_recs/each_tissue_n_BMP_ligands_2.pdf")
par(mfrow=c(3, 3))
for (tissue in unique(my_tissues)[10:18]) {
  dat_n_cells = sum(
		cells_per_sample$n_cells[
			grep(tissue, cells_per_sample$sample)
		]
	)
  dat <- my_dat[
    grep(
			pattern = tissue,
      x = my_dat$tissue,
      fixed = TRUE
		),
  ]
	genes_ON_per_cell <- aggregate(
    dat$gene, 
    by=list(dat$cell), 
    function(x) sum(x %in% BMP_ligands)
  )
  freqs = table(genes_ON_per_cell$x)
	freqs[1] = dat_n_cells - sum(freqs[2:length(freqs)])	
  
	barplot(
    freqs, 
    main = paste0(
      tissue,
      " (n = ", 
      dat_n_cells,
      ")"
    ),
    xlab = "# BMP ligands ON",
    ylab = "# cells",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

# All tissues pooled
#BMP receptors
dat <- my_dat
dat_n_cells <- nrow(counts_per_cell)
genes_ON_per_cell <- aggregate(
	dat$gene, 
	by=list(dat$cell), 
	function(x) sum(x %in% BMP_receptors)
)
freqs = table(genes_ON_per_cell$x)
freqs[1] = dat_n_cells - sum(freqs[2:length(freqs)])	

pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_n_BMP_receptors.pdf")
barplot(
	freqs, 
	main = paste0(
		"All tissues pooled (n = ", 
		sum(cells_per_sample$n_cells),
		")"
	),
	xlab = "# BMP receptors ON",
	ylab = "# cells",
	horiz = F,
	col = "darkslategray"
)
dev.off()

#BMP ligands
dat <- my_dat
dat_n_cells <- nrow(counts_per_cell)
genes_ON_per_cell <- aggregate(
	dat$gene, 
	by=list(dat$cell), 
	function(x) sum(x %in% BMP_ligands)
)
freqs = table(genes_ON_per_cell$x)
freqs[1] = dat_n_cells - sum(freqs[2:length(freqs)])	

pdf("FACS/all_tissues_BMP_ligs_recs/all_tissues_n_BMP_ligands.pdf")
barplot(
	freqs, 
	main = paste0(
		"All tissues pooled (n = ", 
		sum(cells_per_sample$n_cells),
		")"
	),
	xlab = "# BMP ligands ON",
	ylab = "# cells",
	horiz = F,
	col = "darkslategray"
)
dev.off()






















