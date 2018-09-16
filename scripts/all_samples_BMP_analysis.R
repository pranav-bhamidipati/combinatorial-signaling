## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
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
	
	
# Notch_genes <- read.table('./Notch_genes.txt', stringsAsFactors = F)
# Notch_genes <- tolower(Notch_genes$V1)
# Notch_genes_idx <- which(genes %in% Notch_genes)

# Notch_genes_annotated <- read.table('./Notch_genes_annotated.txt', 
# stringsAsFactors = F,
# header = T,
# fill = T)

## Read in the list of data sets
dat_filenames <- list.files(
  'droplet/normalized/', 
  pattern = ".*\\.txt")

## Define the total set of genes to be studied
all_genes <- BMP_genes_annotated$idx[
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

cells_per_sample <- data.frame(
  sample=dat_filenames,
  n_cells = integer(length(dat_filenames))
)

## Extract required data
for (i in 1:length(dat_filenames)) {
  # Select data set for this iteration
  my_dat <- dat_filenames[i]
  
  # Read in data
  dat <- read.table(
    paste0('droplet/normalized/',my_dat), 
    header = T, 
    stringsAsFactors = F
  )
  
  colnames(dat) <- c('gene', 'cell','expr')
  dat$tissue <- rep(my_dat, nrow(dat))
  dat$tissue.cell <- paste0(dat$tissue,'.',dat$cell)
  
  # Keep track of # total cells per sample
  cells_per_sample$n_cells[i] = length(unique(dat$cell))
  
  # Select data for the genes we want
  dat <- dat[dat$gene %in% all_genes,]
  
  # Add the data to the master data set
  master_dat <- rbind.data.frame(master_dat, dat)
}

## Write master_dat and cells_per_sample to file
write.table(
  x = master_dat, 
  file = "droplet/all_tissues_BMP_ligs_recs/BMP_lig_rec_dat.txt", 
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  x = cells_per_sample, 
  file = "droplet/cells_per_sample.txt", 
  row.names = FALSE,
  col.names = TRUE
)

## After adding all the needed data...
## Define list of tissues
my_tissues = c(
  "Bladder","Bladder","Bladder",
  "Heart",
  "Kidney","Kidney","Kidney",
  "Liver","Liver","Liver",
  "Lung","Lung","Lung","Lung",
  "Mammary","Mammary",
  "Marrow","Marrow",
  "Muscle","Muscle",
  "Spleen","Spleen",
  "Thymus",
  "Tongue","Tongue","Tongue",
  "Trachea","Trachea"
)

## Define list of cells
all_cells <- sort(unique(master_dat$tissue.cell))

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
    c(sum(x > 10),
      min(x[x > 10]),
      mean(x[x > 10]),
      median(x[x > 10]),
      max(x)
    )
  }
)
BMP_stats <- cbind (BMP_stats[1],as.data.frame(BMP_stats[,2]))
BMP_stats$Group.1 <- genes[BMP_stats$Group.1]
colnames(BMP_stats) <- c('gene', 
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
rownames(BMP_stats) <- 1:nrow(BMP_stats)
BMP_stats <- BMP_stats[
  na.omit(
    match(
      BMP_genes_annotated$gene,
      BMP_stats$gene
    )
  ),
  ]

# Save to file
write.table(
  x = BMP_stats[order(BMP_stats$role, BMP_stats$gene),], 
  file = "droplet/all_tissues_BMP_ligs_recs/BMP_stats.txt", 
  row.names = FALSE,
  col.names = TRUE
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
				c(sum(x > 10),
					min(x[x > 10]),
					mean(x[x > 10]),
					median(x[x > 10]),
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
		g.cells = dat$tissue.cell[dat$gene == g]
		as.integer(all_cells %in% g.cells)
	}
)
BMP_rec_pa <- as.data.frame(BMP_rec_pa)
colnames(BMP_rec_pa) <- genes[BMP_receptors]
rownames(BMP_rec_pa) <- all_cells

# BMP ligands
dat <- my_dat[my_dat$expr > 10,]
BMP_lig_pa <- sapply(
	X = BMP_ligands,
	FUN = function(g) {
		g.cells = dat$tissue.cell[dat$gene == g]
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

## Plot histograms of expression for each BMP receptor	
pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_rec_expr_1.pdf")
par(mfrow=c(3, 3))
for (rec in BMP_receptors[1:9]) {
  x = my_dat$expr[my_dat$gene == rec]
  x = c(x, numeric(sum(cells_per_sample$n_cells) - length(x)))
  x = log10(x + 1)
  y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
  plot(
    x = seq(-0.25, 6.25, 0.5),
    y = log10(y$counts + 1), 
    # log="y", 
    type='h', 
    lwd=10, 
    lend=2,
    main=paste0(
      genes[rec], 
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

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_rec_expr_2.pdf")
par(mfrow=c(3, 3))
for (rec in BMP_receptors[10:18]) {
  x = my_dat$expr[my_dat$gene == rec]
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
      genes[rec], 
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

## Histograms for each ligand
pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_1.pdf")
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

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_2.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[10:18]) {
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

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_3.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[19:27]) {
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

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_lig_expr_4.pdf")
par(mfrow=c(3, 3))
for (lig in BMP_ligands[28:29]) {
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

## Plot expression histograms for each rec/lig in each tissue
for (tissue in unique(my_tissues)) {
  
  tissue_rows = grepl(tissue, my_dat$tissue)
  tissue_samples = grepl(tissue, cells_per_sample$sample)
  
  # BMP Receptors
  # 1
  pdf(paste0(
    "droplet/each_tissue_BMP_ligs_recs/",
    tissue,
    '/',
    tissue,
    '_BMP_rec_expr_1',
    ".pdf"))
  par(mfrow=c(3, 3))
  for (rec in BMP_receptors[1:9]) {
    x = my_dat$expr[my_dat$gene == rec & tissue_rows]
    x = c(x, 
          numeric(
            sum(cells_per_sample$n_cells[tissue_samples]) - length(x)
          )
    )
    x = log10(x + 1)
    y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
    plot(
      x = seq(-0.25, 6.25, 0.5),
      y = log10(y$counts + 1), 
      type='h', 
      lwd=10, 
      lend=2,
      main=paste0(
        tissue,
        ':',
        genes[rec],
        ' (n = ',
        sum(cells_per_sample$n_cells[tissue_samples]),
        ')'
      ),
      xlab = "Log10(expr + 1)",
      ylab = "Log10(# cells + 1)", 
      col="darkslategray"
    )
  }
  dev.off()
  
  # 2
  pdf(paste0(
    "droplet/each_tissue_BMP_ligs_recs/",
    tissue,
    '/',
    tissue,
    '_BMP_rec_expr_2',
    ".pdf"))
  par(mfrow=c(3, 3))
  for (rec in BMP_receptors[10:18]) {
    x = my_dat$expr[my_dat$gene == rec & tissue_rows]
    x = c(x, 
          numeric(
            sum(cells_per_sample$n_cells[tissue_samples]) - length(x)
          )
    )
    x = log10(x + 1)
    y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
    plot(
      x = seq(-0.25, 6.25, 0.5),
      y = log10(y$counts + 1), 
      type='h', 
      lwd=10, 
      lend=2,
      main=paste0(
        tissue,
        ':',
        genes[rec],
        ' (n = ',
        sum(cells_per_sample$n_cells[tissue_samples]),
        ')'
      ),
      xlab = "Log10(expr + 1)",
      ylab = "Log10(# cells + 1)", 
      col="darkslategray"
    )
  }
  dev.off()
  
  # BMP Ligands
  # 1
  pdf(paste0(
    "droplet/each_tissue_BMP_ligs_recs/",
    tissue,
    '/',
    tissue,
    '_BMP_lig_expr_1',
    ".pdf"))
  par(mfrow=c(3, 3))
  for (lig in BMP_ligands[1:9]) {
    x = my_dat$expr[my_dat$gene == lig & tissue_rows]
    x = c(x, 
          numeric(
            sum(cells_per_sample$n_cells[tissue_samples]) - length(x)
          )
    )
    x = log10(x + 1)
    y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
    plot(
      x = seq(-0.25, 6.25, 0.5),
      y = log10(y$counts + 1), 
      type='h', 
      lwd=10, 
      lend=2,
      main=paste0(
        tissue,
        ':',
        genes[lig],
        ' (n = ',
        sum(cells_per_sample$n_cells[tissue_samples]),
        ')'
      ),
      xlab = "Log10(expr + 1)",
      ylab = "Log10(# cells + 1)", 
      col="darkslategray"
    )
  }
  dev.off()
  
  # 2
  pdf(paste0(
    "droplet/each_tissue_BMP_ligs_recs/",
    tissue,
    '/',
    tissue,
    '_BMP_lig_expr_2',
    ".pdf"))
  par(mfrow=c(3, 3))
  for (lig in BMP_ligands[10:18]) {
    x = my_dat$expr[my_dat$gene == lig & tissue_rows]
    x = c(x, 
          numeric(
            sum(cells_per_sample$n_cells[tissue_samples]) - length(x)
          )
    )
    x = log10(x + 1)
    y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
    plot(
      x = seq(-0.25, 6.25, 0.5),
      y = log10(y$counts + 1), 
      type='h', 
      lwd=10, 
      lend=2,
      main=paste0(
        tissue,
        ':',
        genes[lig],
        ' (n = ',
        sum(cells_per_sample$n_cells[tissue_samples]),
        ')'
      ),
      xlab = "Log10(expr + 1)",
      ylab = "Log10(# cells + 1)", 
      col="darkslategray"
    )
  }
  dev.off()
  
  # 3
  pdf(paste0(
    "droplet/each_tissue_BMP_ligs_recs/",
    tissue,
    '/',
    tissue,
    '_BMP_lig_expr_3',
    ".pdf"))
  par(mfrow=c(3, 3))
  for (lig in BMP_ligands[19:27]) {
    x = my_dat$expr[my_dat$gene == lig & tissue_rows]
    x = c(x, 
          numeric(
            sum(cells_per_sample$n_cells[tissue_samples]) - length(x)
          )
    )
    x = log10(x + 1)
    y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
    plot(
      x = seq(-0.25, 6.25, 0.5),
      y = log10(y$counts + 1), 
      type='h', 
      lwd=10, 
      lend=2,
      main=paste0(
        tissue,
        ':',
        genes[lig],
        ' (n = ',
        sum(cells_per_sample$n_cells[tissue_samples]),
        ')'
      ),
      xlab = "Log10(expr + 1)",
      ylab = "Log10(# cells + 1)", 
      col="darkslategray"
    )
  }
  dev.off()
  
  # 4
  pdf(paste0(
    "droplet/each_tissue_BMP_ligs_recs/",
    tissue,
    '/',
    tissue,
    '_BMP_lig_expr_4',
    ".pdf"))
  par(mfrow=c(3, 3))
  for (lig in BMP_ligands[28:29]) {
    x = my_dat$expr[my_dat$gene == lig & tissue_rows]
    x = c(x, 
          numeric(
            sum(cells_per_sample$n_cells[tissue_samples]) - length(x)
          )
    )
    x = log10(x + 1)
    y = hist(x, breaks=seq(-0.5, 6.5, 0.5), plot=F)
    plot(
      x = seq(-0.25, 6.25, 0.5),
      y = log10(y$counts + 1), 
      type='h', 
      lwd=10, 
      lend=2,
      main=paste0(
        tissue,
        ':',
        genes[lig],
        ' (n = ',
        sum(cells_per_sample$n_cells[tissue_samples]),
        ')'
      ),
      xlab = "Log10(expr + 1)",
      ylab = "Log10(# cells + 1)", 
      col="darkslategray"
    )
  }
  dev.off()
  
}


## Plot barplots of genes expressed per cell in each tissue

# Each tissue separate

#BMP receptors
pdf("droplet/each_tissue_BMP_ligs_recs/each_tissue_n_BMP_receptors.pdf")
par(mfrow=c(3, 3))
for (tissue in unique(my_tissues)) {
  dat_n_cells = sum(
		cells_per_sample$n_cells[
			grep(tissue, cells_per_sample$sample)
		]
	)
  dat <- my_dat[
    grepl(pattern = tissue,
          x = my_dat$tissue,
          fixed = TRUE),
    ]
	dat <- dat[dat$expr > 10,] #Filter expression 
	BMP_receptors_per_cell <- aggregate(
    dat$gene, 
    by=list(dat$tissue.cell), 
    function(x) sum(x %in% BMP_receptors)
  )
  freqs = table(BMP_receptors_per_cell$x)
  
  barplot(
    freqs, 
    main = paste0(
      tissue,
      " (n = ", 
      dat_n_cells,
      ")"
    ),
    xlab = "# BMP receptors",
    ylab = "# cells",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

#BMP ligands
pdf("droplet/each_tissue_BMP_ligs_recs/each_tissue_n_BMP_ligands.pdf")
par(mfrow=c(4, 3))
for (tissue in unique(my_tissues)) {
  dat_n_cells = sum(
		cells_per_sample$n_cells[
			grep(tissue, cells_per_sample$sample)
		]
	)
  dat <- my_dat[
    grepl(pattern = tissue,
          x = my_dat$tissue,
          fixed = TRUE),
    ]
	dat <- dat[dat$expr > 10,] #Filter expression 
	BMP_ligands_per_cell <- aggregate(
    dat$gene, 
    by=list(dat$tissue.cell), 
    function(x) sum(x %in% BMP_ligands)
  )
  freqs = table(BMP_ligands_per_cell$x)
  
  barplot(
    freqs, 
    main = paste0(
      tissue,
      " (n = ", 
      dat_n_cells,
      ")"
    ),
    xlab = "# BMP ligands",
    ylab = "# cells",
    horiz = F,
		col = "darkslategray"
  )
}
dev.off()

# All tissues pooled
#BMP receptors
dat <- my_dat
BMP_receptors_per_cell <- aggregate(
	dat$gene, 
	by=list(dat$tissue.cell), 
	function(x) sum(x %in% BMP_receptors)
)
freqs = table(BMP_receptors_per_cell$x)

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_n_BMP_receptors.pdf")
barplot(
	freqs, 
	main = paste0(
		"All tissues pooled (n = ", 
		sum(cells_per_sample$n_cells),
		")"
	),
	xlab = "# BMP receptors",
	ylab = "# cells",
	horiz = F,
	col = "darkslategray"
)
dev.off()

#BMP ligands
dat <- my_dat
BMP_ligands_per_cell <- aggregate(
	dat$gene, 
	by=list(dat$tissue.cell), 
	function(x) sum(x %in% BMP_ligands)
)
freqs = table(BMP_ligands_per_cell$x)

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_n_BMP_ligands.pdf")
barplot(
	freqs, 
	main = paste0(
		"All tissues pooled (n = ", 
		sum(cells_per_sample$n_cells),
		")"
	),
	xlab = "# BMP ligands",
	ylab = "# cells",
	horiz = F,
	col = "darkslategray"
)
dev.off()

## Make a box and whisker plot of BMP expression by # receptors expressed
dat <- my_dat
BMP_receptors_per_cell <- aggregate(
	dat$gene, 
	by=list(dat$tissue.cell), 
	function(x) sum(x %in% BMP_receptors)
)
dat$n_BMP_rec <- BMP_receptors_per_cell$x[
	match(
		dat$tissue.cell, 
		BMP_receptors_per_cell$Group.1
	)
]
dat <- dat[dat$gene %in% BMP_receptors,]

pdf("droplet/test.pdf")
boxplot(
	log10(expr + 1) ~ n_BMP_rec, 
	data = dat,
  xlab = "# BMP receptors", 
	ylab = "Log(expr + 1)",
  main = "Pooled BMP Receptor Expression by # Receptors Expressed"
)
dev.off()







## Backup code for histogram formatting 
par(mfrow=c(3, 3))
colnames <- my_genes[1:9]
for (i in 1:9) {
  hist(my_dat$expr[my_dat$gene == my_genes[i],], 
       # xlim=c(0, 3500), 
       # breaks=seq(0, 3500, 100), 
       main=colnames[i], 
       probability=TRUE, 
       col="gray", 
       border="white")
  d <- density(crime[,i])
  lines(d, col="red")
}

pdf("droplet/all_tissues_BMP_ligs_recs/all_tissues_BMP_rec_expr_2.pdf")
par(mfrow=c(3, 3))
for (rec in BMP_receptors[10:18]) {
  x = log10(my_dat$expr[my_dat$gene == rec])
  hist(
    x = x,
    xlab = "Log10(norm. expr.)",
    ylab = "Proportion of cells",
    # xlim=c(0, 3500), 
    breaks=seq(0, 6, 0.5), 
    main=paste0(
      genes[rec], 
      " (n = ",
      length(x),")"
    ),
    probability=TRUE, 
    col="gray", 
    border="white"
  )
  if (length(x) > 1) {
    d <- density(x)
    lines(d, col="red")
  }
}
dev.off()
