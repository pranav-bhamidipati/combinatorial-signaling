## Do t-SNE stuff

## Get required parameters
my_tissue = "Liver"
my_in_file <- paste0("FACS/normalized/",my_tissue,"_FACS_norm.txt")
dist_fun = 'euclidean'
my_threshold = 10

## Set up environment
# setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('data.table')
usePackage('Rtsne')
usePackage('scales')

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

## Isolate data to be used later
E.dat.BMP <- E.dat[E.dat$gene %in% BMP_receptors,]
all_cells <- unique(E.dat$cell)

## Convert sparse, long matrix to dense, wide matrix
# Define list of all genes in the sample
genes_in_sample = unique(E.dat$gene)

# For each cell, get list of missing genes
missing_genes <- aggregate(
	x = E.dat$gene,
	by = list(E.dat$cell),
	FUN = function(g.list) genes_in_sample[!(genes_in_sample %in% g.list)]
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


########### **PLOTS** ###########

## Plot distribution of distances between cells before smoothing
pdf(paste0("FACS/",my_tissue,"_dists_unsmoothed.pdf"))
x = D
y = hist(x, breaks=seq(0, 1E6, 2E4), plot=F)
plot(
	x = seq(1E4, 9.9E5, 2E4),
	# y = log10(y$counts + 1), 
	y = y$counts,
	# log="y", 
	type='h', 
	lwd=10, 
	lend=2,
	main=paste0(
		"Cell-cell distance, unsmoothed", 
		" (n = ",
		sum(cells_per_sample$n_cells),
		")"
	),
	xlab = "distance",
	ylab = "# cell pairs", 
	col="darkslategray", 
)
dev.off()

## Plot distribution of distances between cells after smoothing
pdf(paste0("FACS/",my_tissue,"_dists_smoothed.pdf"))
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
	# E <- as.data.frame(E, stringsAsFactors = FALSE)
	# dimnames(E) <- E.dimnames
	# E$gene <- rownames(E)
	# writeLines(paste0('...'))
	
	# ## Convert dense, wide matrix back to sparse, long matrix
	# writeLines(paste0('Reformatting to long format'))
	# E <- reshape2::melt(
		# data = E, 
		# id.vars = 'gene',
		# value.name = 'expr'
	# )
	# writeLines(paste0('.'))	
	# E$gene <- as.integer(E$gene)
	# E$variable <- as.character(E$variable)
	# colnames(E) <- c('gene','cell','expr')
	# writeLines(paste0('..'))
	
	# ### Calculate # genes ON per cell after smoothing
	# ## Restrict data to specific genes of interest
	# my_genes = c(BMP_ligands, BMP_receptors)
	# my_filter = 10
	# E <- E[E$gene %in% my_genes,]
	# E <- E[E$expr > my_filter,]
	# writeLines(paste0('...'))
	
	# ## Plot distribution of # genes ON per cell
	# #BMP receptors
	# writeLines(paste0('Plotting distribution for sigma = ',sigma))
	# dat_n_cells <- sum(counts_per_cell$tissue == my_tissue)													
	# genes_ON_per_cell <- aggregate(
		# E$gene, 
		# by=list(E$cell), 
		# function(x) sum(x %in% BMP_receptors)
	# )
	# freqs = table(genes_ON_per_cell$x)

	## Calculate distances between cells after smoothing
	## Compute distances between cols (cells)
	D.smooth <- dist(
		x = t(E),
		method = dist_fun
	)

	## Plot distribution of distances 
	writeLines(paste0('Plotting distribution for sigma = ',sigma))
	dat_n_cells <- sum(counts_per_cell$tissue == my_tissue)													

  x = D.smooth
  # x = c(x, numeric(sum(cells_per_sample$n_cells) - length(x)))
  # x = log10(x + 1)
  y = hist(x, breaks=seq(0, 1E6, 2E4), plot=F)
  plot(
    x = seq(1E4, 9.9E5, 2E4),
    # y = log10(y$counts + 1), 
		y = y$counts,
    # log="y", 
    type='h', 
    lwd=10, 
    lend=2,
    main=paste0(
      "Cell-cell distance ", 
      " (n = ",
      sum(cells_per_sample$n_cells),
      "), sigma = (",
			sigma,
			")"
    ),
    xlab = "distance",
    ylab = "# cell pairs", 
    col="darkslategray", 
  )
	
	
	
	# barplot(
		# freqs, 
		# main = paste0(
			# my_tissue,
			# " smoothed (sigma=",
			# formatC(sigma, format = "e", digits = 1),
			# ")"
		# ),
		# xlab = "# BMP receptors ON",
		# ylab = "# cells",
		# horiz = F,
		# col = "darkslategray"
	# )
	
}
dev.off()
	
	

## Plot 2D t-SNE of expression with different sigmas
## Define parameters
closeness_threshold = 5E4

## Define color scheme 
colors = rainbow(3)
names(colors) = c("self", "close", "far")

## Make plots
pdf(paste0("FACS/",my_tissue,"_smoothed_tSNE_panel.pdf"))
par(mfrow=c(3, 3), oma=c(0,0,2,0))

# Unsmoothed
dists = D
cell_labels = dists[1:(attr(dists, "Size") - 1)] > closeness_threshold
cell_labels = as.integer(as.logical(cell_labels))
cell_labels = c("close","far")[cell_labels + 1]
cell_labels = c('self', cell_labels)
tsne <- Rtsne(
	t(E.dat), 
	dims = 2, 
	perplexity=30, 
	verbose=TRUE, 
	max_iter = 500
)

plot(
	rbind(tsne$Y[-1,], tsne$Y[1,]), #*Plot first cell last so it's on top!*
	pch=c(rep(1, nrow(tsne$Y) - 1), 18),
	cex=c(rep(1, nrow(tsne$Y) - 1), 2),
	col=colors[c(cell_labels[-1], cell_labels[1])],
	lwd=1.5,
	axes=F,
	xlab="",
	ylab="",
	main=paste0(
		"Unsmoothed"
	)
)
	
# Smoothed
	
for (sigma in c(1E4, 2E4, 5E4, 1E5, 2E5, 5E5, 1E6, 2E6)) {
	
	## Generate weights from distances using kernel
	W <- kernFun(
		d = as.matrix(D),
		sigma = sigma
	)
	
	## Compute smoothed expression
	writeLines(paste0('Smoothing with sigma = ',sigma))
	E <- as.matrix (E.dat)
	E.dimnames <- dimnames(E)
	writeLines(paste0('.'))
	E = (E %*% W) %*% diag(1 / colSums(W))
	writeLines(paste0('..'))
	## Find distances between first cell and all other cells
	D.cell = numeric(1)
	for (i in 2:dim(E)[2]) {
		D.cell[i] = dist(t(E[,c(1,i)]), method = dist_fun)
	}
	
	## Compute t-SNE and plot result
	dists = D.cell
	cell_labels = dists > closeness_threshold
	cell_labels = as.integer(as.logical(cell_labels))
	cell_labels = c("close","far")[cell_labels + 1]
	cell_labels = c('self', cell_labels[-1])
	writeLines(paste0('Computing t-SNE'))
	tsne <- Rtsne(
		t(E), 
		dims = 2, 
		perplexity=30, 
		verbose=TRUE, 
		max_iter = 500
	)
	
	writeLines(paste0('Plotting t-SNE'))
	
	plot(
		rbind(tsne$Y[-1,], tsne$Y[1,]), #*Plot first cell last so it's on top!*
		pch=c(rep(1, nrow(tsne$Y) - 1), 18),
		cex=c(rep(1, nrow(tsne$Y) - 1), 2),
		col=colors[c(cell_labels[-1], cell_labels[1])],
		lwd=1.5,
		axes=F,
		xlab="",
		ylab="",
		main=paste0(
			"Sigma = ",
			formatC(sigma, format = "e", digits = 0)
		)
	)
		
}

title("Liver t-SNE; R=BMPR2, G=BMPR1a, B=ACVR2a", outer=TRUE)

dev.off()
	
	
## BMP words for receptors of interest
my_receptors <- c(
	'bmpr2', 
	'bmpr1a',
	'acvr2a'
)

## Get expression for receptors of interest
my_rec_expr <- lapply(
	X = which(genes %in% my_receptors),
	FUN = function(g) {
		ON_cells = match(E.dat.BMP$cell[E.dat.BMP$gene == g], all_cells)
		g.expr = numeric(length(all_cells))
		g.expr[ON_cells] = E.dat.BMP$expr[E.dat.BMP$gene == g]
		return(g.expr)
	}
)
my_rec_expr <- as.data.frame(my_rec_expr)
colnames(my_rec_expr) <- my_receptors
rownames(my_rec_expr) <- all_cells

## Convert log of expression to RGB color values
my_rec_clr = sapply(
	X = my_rec_expr,
	FUN = function(col) (log(col + 1) / max(log(col + 1))) * 255
)
my_rec_clr <- as.data.frame(my_rec_clr)
my_rec_clr = do.call(rgb, c(list(my_rec_clr), list(maxColorValue=255), list(names=all_cells)))
my_rec_clr = my_rec_clr[order(names(my_rec_clr))]

## Plot 2D t-SNE of expression using RGB color-codes
# Unsmoothed
dists = D
tsne <- Rtsne(
	t(E.dat), 
	dims = 2, 
	perplexity=30, 
	verbose=TRUE, 
	max_iter = 500
)

pdf(paste0("FACS/",my_tissue,"_smoothed_tSNE_RGB_panel.pdf"))
par(mfrow=c(3,3), oma=c(0,0,2,0))

plot(
	tsne$Y, 
	col=alpha(my_rec_clr, 0.5),
	lwd=1,
	pch=16,
	cex=1,
	axes=F,
	xlab="",
	ylab="",
	main="Unsmoothed"
)

# Smoothed!!

E.genes <- rownames(E.dat)
E.cells <- colnames(E.dat)

for (sigma in c(1E4, 2E4, 5E4, 1E5, 2E5, 5E5, 1E6, 2E6)) {

	## Generate weights from distances using kernel
	W <- kernFun(
		d = as.matrix(D),
		sigma = sigma
	)
	
	## Compute smoothed expression
	writeLines(paste0('Smoothing with sigma = ',sigma))
	E <- as.matrix (E.dat)
	writeLines(paste0('.'))
	E = (E %*% W) %*% diag(1 / colSums(W))
	writeLines(paste0('..'))
	## Get expression of desired genes (my_receptors)
	my_rec_expr = E[which(E.genes %in% which(genes %in% my_receptors)), ]
	my_rec_expr <- as.data.frame(t(my_rec_expr))
	colnames(my_rec_expr) <- my_receptors
	rownames(my_rec_expr) <- all_cells
	## Convert log of expression to RGB color values
	my_rec_clr = sapply(
		X = my_rec_expr,
		FUN = function(col) (log(col + 1) / max(log(col + 1))) * 255
	)
	my_rec_clr <- as.data.frame(my_rec_clr)
	my_rec_clr = do.call(
		rgb, 
		c(
			list(my_rec_clr), 
			list(names=all_cells),
			list(maxColorValue=255)
		)
	)
	my_rec_clr = my_rec_clr[order(names(my_rec_clr))]
	
	## Compute t-SNE and plot result
	writeLines(paste0('Computing t-SNE'))
	tsne <- Rtsne(
		t(E), 
		dims = 2, 
		perplexity=30, 
		verbose=TRUE, 
		max_iter = 500
	)

	writeLines(paste0('Plotting t-SNE'))

	plot(
		tsne$Y, 
		col=alpha(my_rec_clr, 0.5),
		lwd=1,
		pch=16,
		cex=1,
		axes=F,
		xlab="",
		ylab="",
		main=paste0(
			"Sigma = ",
			formatC(sigma, format = "e", digits = 0)
		)
	)
}

title("Liver t-SNE \nR = BMPR2, G = BMPR1a, B = ACVR2a", outer=T)
dev.off()
	
	
	
	
	
	