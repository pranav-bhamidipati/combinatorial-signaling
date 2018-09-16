#### Filename:			FACS_smoothing.R
#### Author:				Pranav Bhamidipati
#### Description:		Use kernel smoothing to average each 
####									cell's expression with that of its 
####									neighbors in gene expression space.
#### Arguments:			read_path
####									-	path to read expression matrix
####									- must be sparse, long matrix with header
####								  	and without row names
####								  - must be space-delimited
####								  - must have a '.txt' extension
####								write_path
####									- path to write smoothed expression
####									- writes sparse, long matrix with header
####										and without row names
####								  - must have a '.txt' extension
####								  - output will be written as space-delimited
####								dist_fun
####									- used by dist() to calculate the 
####										distance between cells (must be one of 
####										"euclidean", "maximum", "manhattan", 
####										"canberra", "binary" or "minkowski"). 
####								sigma 
####									- parameter applied to kernel function:
####											exp( - ( d / sigma ) ^ 2 )
####										where d = distance between cells i and j
####									- if unsure of an appropriate sigma, the 
####										mean of the distances from dist() may 
####										be a good starting point
#### Notes:					After smoothing, many zero expression 
####									values become non-zeros due to averaging, 
####									and the size of the output file may be 
####									much larger than the input. I've included
####									a threshold option (commented out)
####									so that if the threshold is known,
####									it would save a lot of storage to apply
####									it before saving the smoothed data.


########################## RUN THIS BEFORE RUNNING FULL SCRIPT

## Set up environment
if (!is.element('reshape2', installed.packages()[,1]))
  install.packages('reshape2', dep = TRUE)
require('reshape2', character.only = TRUE)
if (!is.element('data.table', installed.packages()[,1]))
  install.packages('data.table', dep = TRUE)
require('data.table', character.only = TRUE)

########################## ENTER PARAMETERS, THEN SOURCE FILE

## Get required parameters
read_path <- "FACS/normalized/Liver_FACS_norm.txt"
write_path <- "FACS/smoothed/Liver_FACS_smoothed.txt"

dist_fun = 'euclidean'
sigma = 1.0 * 10^4

## Remove comment ('#') to enforce threshold, otherwise 
##		uses default value of expr_threshold = 0
# expr_threshold = 10


########################## SCRIPT

## Read in necessary files
genes <- read.table(file = "genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_receptors <- read.table('./BMP_Receptors.txt', stringsAsFactors = F)
BMP_receptors <- which(genes %in% tolower(unlist(BMP_receptors)))
# BMP_ligands <- read.table('./BMP_Ligands.txt', stringsAsFactors = F)
# BMP_ligands <- which(genes %in% tolower(unlist(BMP_ligands)))
	
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
dat <- read.table(
	file = read_path,
	header = TRUE,
	stringsAsFactors = FALSE
)

## Convert sparse, long matrix to dense, wide matrix
# Define list of all genes in the sample
genes_in_sample = unique(dat$gene)

# For each cell, get list of missing genes
missing_genes <- aggregate(
	x = dat$gene,
	by = list(dat$cell),
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
dat <- rbind(dat, missing_rows)

# Cast matrix to wide format
dat <- dcast(dat, gene ~ cell)
rownames(dat) <- dat$gene
dat <- dat[-1]

## Compute distances between cols (cells)
D <- dist(
	x = t(dat),
	method = dist_fun
)

## Construct kernel
kernFun = function(d,sigma) {
	exp(
		- ( d / sigma ) ^ 2
	)
}

## Generate weights from distances using kernel
W <- kernFun(
	d = as.matrix(D),
	sigma = sigma
)
lol = hist(W)
lol$breaks
lol$counts

## Compute smoothed expression
E <- as.matrix (dat)
E.dimnames <- dimnames(E)

E = (E %*% W) %*% diag(1 / colSums(W))

E <- as.data.frame(E, stringsAsFactors = FALSE)
dimnames(E) <- E.dimnames
E$gene <- rownames(E)

## Convert dense, wide matrix back to sparse, long matrix
E <- reshape2::melt(
	data = E, 
	id.vars = 'gene',
	value.name = 'expr'
)

## Apply threshold if desired. If not, remove zero values
expr_threshold = mget('expr_threshold', ifnotfound = list(0))[[1]]
E <- E[E$expr > expr_threshold,]
E$gene <- as.integer(E$gene)
E$variable <- as.character(E$variable)
colnames(E) <- c('gene','cell','expr')

## Save smoothed data to file
write.table(
	x = E,
	file = write_path,
	row.names = FALSE
)












