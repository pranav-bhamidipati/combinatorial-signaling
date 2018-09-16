
rm(list = ls())

################################ 1

## Get required parameters
my_tissue = "Marrow"
my_in_file <- paste0("FACS/normalized/",my_tissue,"_FACS_norm.txt")
dist_fun = 'euclidean'

## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('data.table')

## Read in necessary files
genes <- read.table(file = "genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_receptors <- read.table('./BMP_Receptors.txt', stringsAsFactors = F)
BMP_receptors <- which(genes %in% tolower(unlist(BMP_receptors)))
BMP_ligands <- read.table('./BMP_Ligands.txt', stringsAsFactors = F)
BMP_ligands <- which(genes %in% tolower(unlist(BMP_ligands)))
	
counts_per_cell <- read.table(
	"FACS/counts_per_cell.txt", 
	stringsAsFactors=F, 
	header=T)
	
## Calculate total # cells per sample
cells_per_sample <- as.data.frame(
	table(counts_per_cell$tissue),
	stringsAsFactors = FALSE
)
colnames(cells_per_sample) <- c('sample','n_cells')

## Read in data set
E.dat <- read.table(
	file = my_in_file,
	header = TRUE,
	stringsAsFactors = FALSE
)

## Convert sparse, long matrix to dense, wide matrix
# Define list of all genes in the sample
genes = unique(E.dat$gene)

# For each cell, get list of missing genes
missing_genes <- aggregate(
	x = E.dat$gene,
	by = list(E.dat$cell),
	FUN = function(g.list) genes[!(genes %in% g.list)]
)
cells = missing_genes[[1]]
missing_genes = missing_genes[[2]]
names(missing_genes) <- cells

# Get all missing rows
missing_rows = mapply(
	function(cell, g.list) {
		data.frame(
			gene=g.list,
			cell=cell,
			expr=0,
			stringsAsFactors = FALSE
		)
	},
	cells,
	missing_genes,
	SIMPLIFY = FALSE
)
missing_rows <- as.data.frame(data.table::rbindlist(missing_rows))

# Add missing rows to data matrix
E.dat <- rbind(E.dat, missing_rows)

# Cast matrix to wide format
E.dat <- dcast(E.dat, gene ~ cell)
rownames(E.dat) <- E.dat$gene
E.dat <- E.dat[-1]

## Construct kernel
kernFun = function(d,sigma) {
	exp(
		- ( d / sigma ) ^ 2
	)
}

## Compute distances between cols (cells)
D <- dist(
	x = t(E.dat),
	method = dist_fun
)

pdf(paste0("FACS/",my_tissue,"_smoothed_n_BMP_receptors.pdf"))
par(mfrow=c(3,2))
for (sigma in c(1E4, 2E4, 5E4, 1E5, 2E5, 5E5)) {

	## Generate weights from distances using kernel
	W <- kernFun(
		d = as.matrix(D),
		sigma = sigma
		# sigma = mean(D)
	)
	
	## Compute smoothed expression
	writeLines(paste0('Smoothing with sigma = ',sigma))
	E <- as.matrix (E.dat)
	E.dimnames <- dimnames(E)
	writeLines(paste0('.'))
	E = (E %*% W) %*% diag(1 / colSums(W))
	writeLines(paste0('..'))
	E <- as.data.frame(E, stringsAsFactors = FALSE)
	dimnames(E) <- E.dimnames
	E$gene <- rownames(E)
	writeLines(paste0('...'))
	
	## Convert dense, wide matrix back to sparse, long matrix
	writeLines(paste0('Reformatting to long format'))
	E <- reshape2::melt(
		data = E, 
		id.vars = 'gene',
		value.name = 'expr'
	)
	writeLines(paste0('.'))	
	E$gene <- as.integer(E$gene)
	E$variable <- as.character(E$variable)
	colnames(E) <- c('gene','cell','expr')
	writeLines(paste0('..'))
	
	### Calculate # genes ON per cell after smoothing
	## Restrict data to specific genes of interest
	my_genes = c(BMP_ligands, BMP_receptors)
	my_filter = 10
	E <- E[E$gene %in% my_genes,]
	E <- E[E$expr > my_filter,]
	writeLines(paste0('...'))
	
	## Plot distribution of # genes ON per cell
	#BMP receptors
	writeLines(paste0('Plotting distribution for sigma = ',sigma))
	dat_n_cells <- sum(counts_per_cell$tissue == my_tissue)													
	genes_ON_per_cell <- aggregate(
		E$gene, 
		by=list(E$cell), 
		function(x) sum(x %in% BMP_receptors)
	)
	freqs = table(genes_ON_per_cell$x)

	barplot(
		freqs, 
		main = paste0(
			my_tissue,
			" smoothed (sigma=",
			formatC(sigma, format = "e", digits = 1),
			")"
		),
		xlab = "# BMP receptors ON",
		ylab = "# cells",
		horiz = F,
		col = "darkslategray"
	)
	
}
dev.off()
	
	
	
################################ end 1

rm(list = ls())

################################ 2

## Get required parameters
my_tissue = "Brain_Neurons"
my_in_file <- paste0("FACS/normalized/",my_tissue,"_FACS_norm.txt")
dist_fun = 'euclidean'

## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('data.table')

## Read in necessary files
genes <- read.table(file = "genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_receptors <- read.table('./BMP_Receptors.txt', stringsAsFactors = F)
BMP_receptors <- which(genes %in% tolower(unlist(BMP_receptors)))
BMP_ligands <- read.table('./BMP_Ligands.txt', stringsAsFactors = F)
BMP_ligands <- which(genes %in% tolower(unlist(BMP_ligands)))
	
counts_per_cell <- read.table(
	"FACS/counts_per_cell.txt", 
	stringsAsFactors=F, 
	header=T)
	
## Calculate total # cells per sample
cells_per_sample <- as.data.frame(
	table(counts_per_cell$tissue),
	stringsAsFactors = FALSE
)
colnames(cells_per_sample) <- c('sample','n_cells')

## Read in data set
E.dat <- read.table(
	file = my_in_file,
	header = TRUE,
	stringsAsFactors = FALSE
)

## Convert sparse, long matrix to dense, wide matrix
# Define list of all genes in the sample
genes = unique(E.dat$gene)

# For each cell, get list of missing genes
missing_genes <- aggregate(
	x = E.dat$gene,
	by = list(E.dat$cell),
	FUN = function(g.list) genes[!(genes %in% g.list)]
)
cells = missing_genes[[1]]
missing_genes = missing_genes[[2]]
names(missing_genes) <- cells

# Get all missing rows
missing_rows = mapply(
	function(cell, g.list) {
		data.frame(
			gene=g.list,
			cell=cell,
			expr=0,
			stringsAsFactors = FALSE
		)
	},
	cells,
	missing_genes,
	SIMPLIFY = FALSE
)
missing_rows <- as.data.frame(data.table::rbindlist(missing_rows))

# Add missing rows to data matrix
E.dat <- rbind(E.dat, missing_rows)

# Cast matrix to wide format
E.dat <- dcast(E.dat, gene ~ cell)
rownames(E.dat) <- E.dat$gene
E.dat <- E.dat[-1]

## Construct kernel
kernFun = function(d,sigma) {
	exp(
		- ( d / sigma ) ^ 2
	)
}

## Compute distances between cols (cells)
D <- dist(
	x = t(E.dat),
	method = dist_fun
)

pdf(paste0("FACS/",my_tissue,"_smoothed_n_BMP_receptors.pdf"))
par(mfrow=c(3,2))
for (sigma in c(1E4, 2E4, 5E4, 1E5, 2E5, 5E5)) {

	## Generate weights from distances using kernel
	W <- kernFun(
		d = as.matrix(D),
		sigma = sigma
		# sigma = mean(D)
	)
	
	## Compute smoothed expression
	writeLines(paste0('Smoothing with sigma = ',sigma))
	E <- as.matrix (E.dat)
	E.dimnames <- dimnames(E)
	writeLines(paste0('.'))
	E = (E %*% W) %*% diag(1 / colSums(W))
	writeLines(paste0('..'))
	E <- as.data.frame(E, stringsAsFactors = FALSE)
	dimnames(E) <- E.dimnames
	E$gene <- rownames(E)
	writeLines(paste0('...'))
	
	## Convert dense, wide matrix back to sparse, long matrix
	writeLines(paste0('Reformatting to long format'))
	E <- reshape2::melt(
		data = E, 
		id.vars = 'gene',
		value.name = 'expr'
	)
	writeLines(paste0('.'))	
	E$gene <- as.integer(E$gene)
	E$variable <- as.character(E$variable)
	colnames(E) <- c('gene','cell','expr')
	writeLines(paste0('..'))
	
	### Calculate # genes ON per cell after smoothing
	## Restrict data to specific genes of interest
	my_genes = c(BMP_ligands, BMP_receptors)
	my_filter = 10
	E <- E[E$gene %in% my_genes,]
	E <- E[E$expr > my_filter,]
	writeLines(paste0('...'))
	
	## Plot distribution of # genes ON per cell
	#BMP receptors
	writeLines(paste0('Plotting distribution for sigma = ',sigma))
	dat_n_cells <- sum(counts_per_cell$tissue == my_tissue)													
	genes_ON_per_cell <- aggregate(
		E$gene, 
		by=list(E$cell), 
		function(x) sum(x %in% BMP_receptors)
	)
	freqs = table(genes_ON_per_cell$x)

	barplot(
		freqs, 
		main = paste0(
			my_tissue,
			" smoothed (sigma=",
			formatC(sigma, format = "e", digits = 1),
			")"
		),
		xlab = "# BMP receptors ON",
		ylab = "# cells",
		horiz = F,
		col = "darkslategray"
	)
	
}
dev.off()
	
	

################################ end 2

rm(list = ls())

################################ 3

## Get required parameters
my_tissue = "Fat"
my_in_file <- paste0("FACS/normalized/",my_tissue,"_FACS_norm.txt")
dist_fun = 'euclidean'

## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('data.table')

## Read in necessary files
genes <- read.table(file = "genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_receptors <- read.table('./BMP_Receptors.txt', stringsAsFactors = F)
BMP_receptors <- which(genes %in% tolower(unlist(BMP_receptors)))
BMP_ligands <- read.table('./BMP_Ligands.txt', stringsAsFactors = F)
BMP_ligands <- which(genes %in% tolower(unlist(BMP_ligands)))
	
counts_per_cell <- read.table(
	"FACS/counts_per_cell.txt", 
	stringsAsFactors=F, 
	header=T)
	
## Calculate total # cells per sample
cells_per_sample <- as.data.frame(
	table(counts_per_cell$tissue),
	stringsAsFactors = FALSE
)
colnames(cells_per_sample) <- c('sample','n_cells')

## Read in data set
E.dat <- read.table(
	file = my_in_file,
	header = TRUE,
	stringsAsFactors = FALSE
)

## Convert sparse, long matrix to dense, wide matrix
# Define list of all genes in the sample
genes = unique(E.dat$gene)

# For each cell, get list of missing genes
missing_genes <- aggregate(
	x = E.dat$gene,
	by = list(E.dat$cell),
	FUN = function(g.list) genes[!(genes %in% g.list)]
)
cells = missing_genes[[1]]
missing_genes = missing_genes[[2]]
names(missing_genes) <- cells

# Get all missing rows
missing_rows = mapply(
	function(cell, g.list) {
		data.frame(
			gene=g.list,
			cell=cell,
			expr=0,
			stringsAsFactors = FALSE
		)
	},
	cells,
	missing_genes,
	SIMPLIFY = FALSE
)
missing_rows <- as.data.frame(data.table::rbindlist(missing_rows))

# Add missing rows to data matrix
E.dat <- rbind(E.dat, missing_rows)

# Cast matrix to wide format
E.dat <- dcast(E.dat, gene ~ cell)
rownames(E.dat) <- E.dat$gene
E.dat <- E.dat[-1]

## Construct kernel
kernFun = function(d,sigma) {
	exp(
		- ( d / sigma ) ^ 2
	)
}

## Compute distances between cols (cells)
D <- dist(
	x = t(E.dat),
	method = dist_fun
)

pdf(paste0("FACS/",my_tissue,"_smoothed_n_BMP_receptors.pdf"))
par(mfrow=c(3,2))
for (sigma in c(1E4, 2E4, 5E4, 1E5, 2E5, 5E5)) {

	## Generate weights from distances using kernel
	W <- kernFun(
		d = as.matrix(D),
		sigma = sigma
		# sigma = mean(D)
	)
	
	## Compute smoothed expression
	writeLines(paste0('Smoothing with sigma = ',sigma))
	E <- as.matrix (E.dat)
	E.dimnames <- dimnames(E)
	writeLines(paste0('.'))
	E = (E %*% W) %*% diag(1 / colSums(W))
	writeLines(paste0('..'))
	E <- as.data.frame(E, stringsAsFactors = FALSE)
	dimnames(E) <- E.dimnames
	E$gene <- rownames(E)
	writeLines(paste0('...'))
	
	## Convert dense, wide matrix back to sparse, long matrix
	writeLines(paste0('Reformatting to long format'))
	E <- reshape2::melt(
		data = E, 
		id.vars = 'gene',
		value.name = 'expr'
	)
	writeLines(paste0('.'))	
	E$gene <- as.integer(E$gene)
	E$variable <- as.character(E$variable)
	colnames(E) <- c('gene','cell','expr')
	writeLines(paste0('..'))
	
	### Calculate # genes ON per cell after smoothing
	## Restrict data to specific genes of interest
	my_genes = c(BMP_ligands, BMP_receptors)
	my_filter = 10
	E <- E[E$gene %in% my_genes,]
	E <- E[E$expr > my_filter,]
	writeLines(paste0('...'))
	
	## Plot distribution of # genes ON per cell
	#BMP receptors
	writeLines(paste0('Plotting distribution for sigma = ',sigma))
	dat_n_cells <- sum(counts_per_cell$tissue == my_tissue)													
	genes_ON_per_cell <- aggregate(
		E$gene, 
		by=list(E$cell), 
		function(x) sum(x %in% BMP_receptors)
	)
	freqs = table(genes_ON_per_cell$x)

	barplot(
		freqs, 
		main = paste0(
			my_tissue,
			" smoothed (sigma=",
			formatC(sigma, format = "e", digits = 1),
			")"
		),
		xlab = "# BMP receptors ON",
		ylab = "# cells",
		horiz = F,
		col = "darkslategray"
	)
	
}
dev.off()
	
	

################################ end 3

rm(list = ls())

################################ 4

## Get required parameters
my_tissue = "Heart"
my_in_file <- paste0("FACS/normalized/",my_tissue,"_FACS_norm.txt")
dist_fun = 'euclidean'

## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('data.table')

## Read in necessary files
genes <- read.table(file = "genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_receptors <- read.table('./BMP_Receptors.txt', stringsAsFactors = F)
BMP_receptors <- which(genes %in% tolower(unlist(BMP_receptors)))
BMP_ligands <- read.table('./BMP_Ligands.txt', stringsAsFactors = F)
BMP_ligands <- which(genes %in% tolower(unlist(BMP_ligands)))
	
counts_per_cell <- read.table(
	"FACS/counts_per_cell.txt", 
	stringsAsFactors=F, 
	header=T)
	
## Calculate total # cells per sample
cells_per_sample <- as.data.frame(
	table(counts_per_cell$tissue),
	stringsAsFactors = FALSE
)
colnames(cells_per_sample) <- c('sample','n_cells')

## Read in data set
E.dat <- read.table(
	file = my_in_file,
	header = TRUE,
	stringsAsFactors = FALSE
)

## Convert sparse, long matrix to dense, wide matrix
# Define list of all genes in the sample
genes = unique(E.dat$gene)

# For each cell, get list of missing genes
missing_genes <- aggregate(
	x = E.dat$gene,
	by = list(E.dat$cell),
	FUN = function(g.list) genes[!(genes %in% g.list)]
)
cells = missing_genes[[1]]
missing_genes = missing_genes[[2]]
names(missing_genes) <- cells

# Get all missing rows
missing_rows = mapply(
	function(cell, g.list) {
		data.frame(
			gene=g.list,
			cell=cell,
			expr=0,
			stringsAsFactors = FALSE
		)
	},
	cells,
	missing_genes,
	SIMPLIFY = FALSE
)
missing_rows <- as.data.frame(data.table::rbindlist(missing_rows))

# Add missing rows to data matrix
E.dat <- rbind(E.dat, missing_rows)

# Cast matrix to wide format
E.dat <- dcast(E.dat, gene ~ cell)
rownames(E.dat) <- E.dat$gene
E.dat <- E.dat[-1]

## Construct kernel
kernFun = function(d,sigma) {
	exp(
		- ( d / sigma ) ^ 2
	)
}

## Compute distances between cols (cells)
D <- dist(
	x = t(E.dat),
	method = dist_fun
)

pdf(paste0("FACS/",my_tissue,"_smoothed_n_BMP_receptors.pdf"))
par(mfrow=c(3,2))
for (sigma in c(1E4, 2E4, 5E4, 1E5, 2E5, 5E5)) {

	## Generate weights from distances using kernel
	W <- kernFun(
		d = as.matrix(D),
		sigma = sigma
		# sigma = mean(D)
	)
	
	## Compute smoothed expression
	writeLines(paste0('Smoothing with sigma = ',sigma))
	E <- as.matrix (E.dat)
	E.dimnames <- dimnames(E)
	writeLines(paste0('.'))
	E = (E %*% W) %*% diag(1 / colSums(W))
	writeLines(paste0('..'))
	E <- as.data.frame(E, stringsAsFactors = FALSE)
	dimnames(E) <- E.dimnames
	E$gene <- rownames(E)
	writeLines(paste0('...'))
	
	## Convert dense, wide matrix back to sparse, long matrix
	writeLines(paste0('Reformatting to long format'))
	E <- reshape2::melt(
		data = E, 
		id.vars = 'gene',
		value.name = 'expr'
	)
	writeLines(paste0('.'))	
	E$gene <- as.integer(E$gene)
	E$variable <- as.character(E$variable)
	colnames(E) <- c('gene','cell','expr')
	writeLines(paste0('..'))
	
	### Calculate # genes ON per cell after smoothing
	## Restrict data to specific genes of interest
	my_genes = c(BMP_ligands, BMP_receptors)
	my_filter = 10
	E <- E[E$gene %in% my_genes,]
	E <- E[E$expr > my_filter,]
	writeLines(paste0('...'))
	
	## Plot distribution of # genes ON per cell
	#BMP receptors
	writeLines(paste0('Plotting distribution for sigma = ',sigma))
	dat_n_cells <- sum(counts_per_cell$tissue == my_tissue)													
	genes_ON_per_cell <- aggregate(
		E$gene, 
		by=list(E$cell), 
		function(x) sum(x %in% BMP_receptors)
	)
	freqs = table(genes_ON_per_cell$x)

	barplot(
		freqs, 
		main = paste0(
			my_tissue,
			" smoothed (sigma=",
			formatC(sigma, format = "e", digits = 1),
			")"
		),
		xlab = "# BMP receptors ON",
		ylab = "# cells",
		horiz = F,
		col = "darkslategray"
	)
	
}
dev.off()
	
	

	
	
	
	
	
	
	
	
	
	
	