## Set up environment
setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('heatmap3')
usePackage('colorRamps')
usePackage('lattice') 

## Read in droplet data set (marrow)
dat_filenames <- list.files(
  'droplet/normalized/', 
  pattern = "Marrow.*\\.txt")
dat <- lapply(
  dat_filenames,
  function(x) {
    read.table(
      paste0('droplet/normalized/',x), 
      header = T, 
      stringsAsFactors = F
    )
  }
)
names(dat) <- dat_filenames

## Melt data into one df
dat_df <- mapply(
  function(x, y) {
    cbind(x, tissue=rep(y, nrow(x)))
  },
  dat,
  names(dat),
  SIMPLIFY = F
)
dat_df = do.call(rbind.data.frame, dat_df)
colnames(dat_df) <- c('gene', 'cell','expr','tissue')
rownames(dat_df) <- c()
dat_df$tissue.cell <- paste0(dat_df$tissue,'.',dat_df$cell)

## Read in master gene list and signaling gene list
genes <- read.table(file = "droplet/genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)

BMP_genes <- read.table('./BMP_genes.txt', stringsAsFactors = F)
BMP_genes <- tolower(BMP_genes$V1)
BMP_genes_idx <- which(genes %in% BMP_genes)

Notch_genes <- read.table('./Notch_genes.txt', stringsAsFactors = F)
Notch_genes <- tolower(Notch_genes$V1)
Notch_genes_idx <- which(genes %in% Notch_genes)

BMP_genes_annotated <- read.table('./BMP_genes_annotated.txt', 
                                stringsAsFactors = F,
                                header = T)

Notch_genes_annotated <- read.table('./Notch_genes_annotated.txt', 
                                    stringsAsFactors = F,
                                    header = T,
                                    fill = T)

## Compute metrics on genes
all_expr <- dat_df[,c('gene','expr')]
all_expr <- all_expr[order(all_expr[,1]),]
rownames(all_expr) <- c()
colnames(all_expr) <- c('gene','expr')
means = aggregate(x = all_expr$expr, by = list(all_expr$gene), FUN = mean)
colnames(means) <- c('gene','mean_expr')
vars = aggregate(x = all_expr$expr, by = list(all_expr$gene), FUN = var)
colnames(vars) <- c('gene','var_expr')
cvs = aggregate(
  x = all_expr$expr, 
  by = list(all_expr$gene), 
  FUN = function(x) var(x) / mean(x) * 100
  )
colnames(cvs) <- c('gene','CV_expr')
sds = aggregate(x = all_expr$expr, by = list(all_expr$gene), FUN = sd)
colnames(sds) <- c('gene', 'sd_expr')

#BMP
BMP_means <- means[means$gene %in% BMP_genes_idx,]
BMP_vars  <- vars[vars$gene %in% BMP_genes_idx,]
BMP_cvs   <- cvs[cvs$gene %in% BMP_genes_idx,]
BMP_sds   <- sds[sds$gene %in% BMP_genes_idx,]

#Notch
Notch_means <- means[means$gene %in% Notch_genes_idx,]
Notch_vars  <- vars[vars$gene %in% Notch_genes_idx,]
Notch_cvs   <- cvs[cvs$gene %in% Notch_genes_idx,]
Notch_sds   <- sds[sds$gene %in% Notch_genes_idx,]

## Clean genes by removing those with NA values or outlier CVs
nrow(BMP_cvs)                             # 41 total genes
sum(!is.na(BMP_cvs$CV_expr))              # 2 NA values
sum(BMP_cvs$CV_expr < 13000, na.rm = T)   # 6 CV outliers
BMP_clean_idx <- BMP_cvs$gene[which(BMP_cvs$CV_expr < 13000)]
table(BMP_clean_idx)

nrow(Notch_cvs)                             # 42 total genes
sum(!is.na(Notch_cvs$CV_expr))              # 2 NA values
sum(Notch_cvs$CV_expr < 15000, na.rm = T)   # 11 CV outliers
Notch_clean_idx <- Notch_cvs$gene[which(Notch_cvs$CV_expr < 15000)]

## Get data for cleaned gene sets

# Only cleaned data
BMP_clean_dat <- dat_df[dat_df$gene %in% BMP_clean_idx,]
Notch_clean_dat <- dat_df[dat_df$gene %in% Notch_clean_idx,]

# All data (clean and noisy)
BMP_dat <- dat_df[dat_df$gene %in% BMP_genes_idx,]
Notch_dat <- dat_df[dat_df$gene %in% Notch_genes_idx,]

## Analyze how many signaling receptors/ligands expressed per cell
# BMP
BMP_ligands <- BMP_genes_annotated$idx[
  grepl(pattern = "ligand", x = BMP_genes_annotated$role)
  ]
BMP_receptors <- BMP_genes_annotated$idx[
  grepl(pattern = "receptor", x = BMP_genes_annotated$role)
  ]
BMP_smads <- BMP_genes_annotated$idx[
  grepl(pattern = "transducer", x = BMP_genes_annotated$role)
  ]

BMP_ligands_per_cell <- aggregate(
  dat_df$gene, 
  by=list(dat_df$tissue.cell), 
  function(x) sum(x %in% BMP_ligands)
  )
table(BMP_ligands_per_cell[,2])
  
BMP_receptors_per_cell <- aggregate(
  dat_df$gene, 
  by=list(dat_df$tissue.cell), 
  function(x) sum(x %in% BMP_receptors)
)
table(BMP_receptors_per_cell[,2])

BMP_smads_per_cell <- aggregate(
  dat_df$gene, 
  by=list(dat_df$tissue.cell), 
  function(x) sum(x %in% BMP_smads)
)
table(BMP_smads_per_cell[,2])

# Notch
Notch_ligands <- Notch_genes_annotated$idx[
  grepl(pattern = "ligand", x = Notch_genes_annotated$role)
  ]
Notch_receptors <- Notch_genes_annotated$idx[
  grepl(pattern = "receptor", x = Notch_genes_annotated$role)
  ]

Notch_ligands_per_cell <- aggregate(
  dat_df$gene, 
  by=list(dat_df$tissue.cell), 
  function(x) sum(x %in% Notch_ligands)
)
table(Notch_ligands_per_cell[,2])

Notch_receptors_per_cell <- aggregate(
  dat_df$gene, 
  by=list(dat_df$tissue.cell), 
  function(x) sum(x %in% Notch_receptors)
)
table(Notch_receptors_per_cell[,2])

## Plots

# BMP
pdf(file = "droplet/Marrow/barplots-histograms/BMP_ligands_count.pdf")
barplot(
  table(BMP_ligands_per_cell$x), 
  main = paste0("Marrow (n = 4112)\n# BMP ligands expressed per cell"),
  xlab = "# ligands",
  ylab = "# cells",
  horiz = F
)
dev.off()

pdf(file = "droplet/Marrow/barplots-histograms/BMP_receptors_count.pdf")
barplot(
  table(BMP_receptors_per_cell$x), 
  main = paste0("Marrow (n = 4112)\n# BMP receptors expressed per cell"),
  xlab = "# receptors",
  ylab = "# cells",
  horiz = F
)
dev.off()

pdf(file = "droplet/Marrow/barplots-histograms/BMP_smads_count.pdf")
barplot(
  table(BMP_smads_per_cell$x), 
  main = paste0("Marrow (n = 4112)\n# BMP Smads expressed per cell"),
  xlab = "# smads",
  ylab = "# cells",
  horiz = F
)
dev.off()

# Notch
pdf(file = "droplet/Marrow/barplots-histograms/Notch_ligands_count.pdf")
barplot(
  table(Notch_ligands_per_cell$x), 
  main = paste0("Marrow (n = 4112)\n# Notch ligands expressed per cell"),
  xlab = "# ligands",
  ylab = "# cells",
  horiz = F
)
dev.off()

pdf(file = "droplet/Marrow/barplots-histograms/Notch_receptors_count.pdf")
barplot(
  table(Notch_receptors_per_cell$x), 
  main = paste0("Marrow (n = 4112)\n# Notch receptors expressed per cell"),
  xlab = "# receptors",
  ylab = "# cells",
  horiz = F
)
dev.off()

## Make any data frams needed for heatmaps

#BMP ligands and receptors
used_BMP_ligs <- BMP_ligands[BMP_ligands %in% BMP_dat$gene]
BMP_ligs_dat <- BMP_dat[BMP_dat$gene %in% used_BMP_ligs,]

used_BMP_recs <- BMP_receptors[BMP_receptors %in% BMP_dat$gene]
BMP_recs_dat <- BMP_dat[BMP_dat$gene %in% used_BMP_recs,]

#Most important Notch genes (specified by Rachael)
Notch_genes_imp <- Notch_genes_annotated$idx[
  nchar(Notch_genes_annotated$role) > 0 &
    Notch_genes_annotated$role != "protease"
  ]
used_Notch_genes <- Notch_genes_imp[Notch_genes_imp %in% Notch_dat$gene]
Notch_imp_dat <- Notch_dat[Notch_dat$gene %in% Notch_genes_imp,]

## Initialize matrices without any skipped gene x cell combinations

#BMP
BMP_mtx <- data.frame(
  gene=rep(BMP_clean_idx,
           length(unique(BMP_clean_dat$tissue.cell))),
  tissue.cell=rep(unique(BMP_clean_dat$tissue.cell),
                  each=length(BMP_clean_idx)),
  expr=numeric(length(BMP_clean_idx) * length(unique(BMP_clean_dat$tissue.cell))),
  stringsAsFactors = F
)

BMP_ligs_mtx <- data.frame(
  gene=rep(used_BMP_ligs,
           length(unique(BMP_ligs_dat$tissue.cell))),
  tissue.cell=rep(unique(BMP_ligs_dat$tissue.cell),
                  each=length(used_BMP_ligs)),
  expr=numeric(length(used_BMP_ligs) * length(unique(BMP_ligs_dat$tissue.cell))),
  stringsAsFactors = F
)

BMP_recs_mtx <- data.frame(
  gene=rep(used_BMP_recs,
           length(unique(BMP_recs_dat$tissue.cell))),
  tissue.cell=rep(unique(BMP_recs_dat$tissue.cell),
                  each=length(used_BMP_recs)),
  expr=numeric(length(used_BMP_recs) * length(unique(BMP_recs_dat$tissue.cell))),
  stringsAsFactors = F
)

#Notch
Notch_mtx <- data.frame(
  gene=rep(Notch_clean_idx,
           length(unique(Notch_dat$tissue.cell))),
  tissue.cell=rep(unique(Notch_dat$tissue.cell),
                  each=length(Notch_clean_idx)),
  expr=numeric(length(Notch_clean_idx) * length(unique(Notch_dat$tissue.cell))),
  stringsAsFactors = F
)

Notch_imp_mtx <- data.frame(
  gene=rep(used_Notch_genes,
           length(unique(Notch_imp_dat$tissue.cell))),
  tissue.cell=rep(unique(Notch_imp_dat$tissue.cell),
                  each=length(used_Notch_genes)),
  expr=numeric(length(used_Notch_genes) * length(unique(Notch_imp_dat$tissue.cell))),
  stringsAsFactors = F
)

## Fill matrices with signaling expression values

#BMP
BMP_mtx_rows = do.call(paste0, BMP_mtx[c('gene','tissue.cell')])
BMP_clean_dat_rows = do.call(paste0, BMP_clean_dat[c('gene','tissue.cell')])
matching_rows = match(BMP_clean_dat_rows, BMP_mtx_rows)
BMP_mtx$expr[matching_rows] = BMP_clean_dat$expr

BMP_mtx = acast(BMP_mtx, gene ~ tissue.cell, value.var = 'expr')
rownames(BMP_mtx) = genes[as.integer(rownames(BMP_mtx))]
BMP_log_mtx = log10(BMP_mtx+1)

#BMP ligands
BMP_ligs_mtx_rows = do.call(paste0, BMP_ligs_mtx[c('gene','tissue.cell')])
BMP_ligs_dat_rows = do.call(paste0, BMP_ligs_dat[c('gene','tissue.cell')])
matching_rows = match(BMP_ligs_dat_rows, BMP_ligs_mtx_rows)
BMP_ligs_mtx$expr[matching_rows] = BMP_ligs_dat$expr

BMP_ligs_mtx = acast(BMP_ligs_mtx, gene ~ tissue.cell, value.var = 'expr')
rownames(BMP_ligs_mtx) = genes[as.integer(rownames(BMP_ligs_mtx))]
BMP_ligs_log_mtx = log10(BMP_ligs_mtx+1)

#BMP receptors
BMP_recs_mtx_rows = do.call(paste0, BMP_recs_mtx[c('gene','tissue.cell')])
BMP_recs_dat_rows = do.call(paste0, BMP_recs_dat[c('gene','tissue.cell')])
matching_rows = match(BMP_recs_dat_rows, BMP_recs_mtx_rows)
BMP_recs_mtx$expr[matching_rows] = BMP_recs_dat$expr

BMP_recs_mtx = acast(BMP_recs_mtx, gene ~ tissue.cell, value.var = 'expr')
rownames(BMP_recs_mtx) = genes[as.integer(rownames(BMP_recs_mtx))]
BMP_recs_log_mtx = log10(BMP_recs_mtx+1)

#Notch
Notch_mtx_rows = do.call(paste0, Notch_mtx[c('gene','tissue.cell')])
Notch_dat_rows = do.call(paste0, Notch_dat[c('gene','tissue.cell')])
matching_rows = match(Notch_dat_rows, Notch_mtx_rows)
Notch_mtx$expr[matching_rows] = Notch_dat$expr

Notch_mtx = acast(Notch_mtx, gene ~ tissue.cell, value.var = 'expr')
rownames(Notch_mtx) = genes[as.integer(rownames(Notch_mtx))]
Notch_log_mtx = log10(Notch_mtx+1)

#Notch genes marked most important by Rachael
Notch_imp_mtx_rows = do.call(paste0, Notch_imp_mtx[c('gene','tissue.cell')])
Notch_imp_dat_rows = do.call(paste0, Notch_imp_dat[c('gene','tissue.cell')])
matching_rows = match(Notch_imp_dat_rows, Notch_imp_mtx_rows)
Notch_imp_mtx$expr[matching_rows] = Notch_imp_dat$expr

Notch_imp_mtx = acast(Notch_imp_mtx, gene ~ tissue.cell, value.var = 'expr')
rownames(Notch_imp_mtx) = genes[as.integer(rownames(Notch_imp_mtx))]
Notch_imp_log_mtx = log10(Notch_imp_mtx+1)

## Clustering and heatmaps 

# Color palettes

my_palette <- colorRamps::matlab.like(300)

# my_palette <- c(
#   colorRampPalette(c("dark blue", "green"))(n = 20),
#   colorRampPalette(c("green", "red"))(n = 179)
#   )

# my_palette <- colorRamps::matlab.like2(299)

# BMP genes
#Euclidean distance between cols
BMP_dist_euc <- dist(x = t(BMP_mtx),
                     method = "euclidean"
)
BMP_hc_euc <- hclust(d = BMP_dist_euc, 
                     method = "complete"
)
BMP_dend_euc <- as.dendrogram(BMP_hc_euc)

#Euclidean distance with log10(x+1) transformation
BMP_log_dist_euc <- dist(x = t(BMP_log_mtx),
                         method = "euclidean"
)
BMP_log_hc_euc <- hclust(d = BMP_log_dist_euc, 
                         method = "complete"
)
BMP_log_dend_euc <- as.dendrogram(BMP_log_hc_euc)

#Manhattan distance between cols 
BMP_dist_man <- dist(x = t(BMP_mtx),
                     method = "manhattan"
)
BMP_hc_man <- hclust(d = BMP_dist_man, 
                     method = "complete"
)
BMP_dend_man <- as.dendrogram(BMP_hc_man)

#Manhattan distance with log10(x+1) transformation
BMP_log_dist_man <- dist(x = t(BMP_log_mtx),
                         method = "manhattan"
)
BMP_log_hc_man <- hclust(d = BMP_log_dist_man, 
                         method = "complete"
)
BMP_log_dend_man <- as.dendrogram(BMP_log_hc_man)

#Correlation 
BMP_cor <- cor(t(BMP_mtx))
BMP_cor_hc <- hclust(d = as.dist(BMP_cor), method = 'complete')
BMP_cor_dend <- as.dendrogram(BMP_cor_hc)

#Plots
pdf(file="droplet/Marrow/heatmaps/BMP_euclidean.pdf")
BMP_hm <- heatmap3(
  x = BMP_mtx, 
  main = 'Marrow BMP (Euclidean)',
  # Rowv=NA,
  Colv=BMP_dend_euc,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/BMP_log10_euclidean.pdf")
BMP_hm <- heatmap3(
  x = BMP_log_mtx, 
  main = 'Marrow BMP (Log10-Euclidean)',
  # Rowv=NA,
  Colv=BMP_log_dend_euc,
  labCol = NA,
  col = my_palette, 
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/BMP_manhattan.pdf")
BMP_hm <- heatmap3(
  x = BMP_mtx, 
  main = 'Marrow BMP (Manhattan)',
  # Rowv=NA,
  Colv=BMP_dend_man,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/BMP_log10_manhattan.pdf")
BMP_hm <- heatmap3(
  x = BMP_log_mtx, 
  main = 'Marrow BMP (Log10-Manhattan)',
  # Rowv=NA,
  Colv=BMP_log_dend_man,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/BMP_correlation.pdf")
BMP_hm <- heatmap3(
  x = BMP_cor, 
  main = 'Marrow BMP (Correlation)',
  distfun = NA,
  scale = "none",
  showRowDendro = F,
  # Rowv=NA,
  Colv=BMP_cor_dend,
  Rowv=BMP_cor_dend,
  # labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

# BMP ligands

#Euclidean distance between cols
BMP_ligs_dist <- dist(x = t(BMP_ligs_mtx),
                      method = "euclidean"
)
BMP_ligs_hc <- hclust(d = BMP_ligs_dist, 
                      method = "complete"
)
BMP_ligs_dend <- as.dendrogram(BMP_ligs_hc)

#Euclidean distance with log10(x+1) transformation
BMP_ligs_log_dist <- dist(x = t(BMP_ligs_log_mtx),
                          method = "euclidean"
)
BMP_ligs_log_hc <- hclust(d = BMP_ligs_log_dist, 
                          method = "complete"
)
BMP_ligs_log_dend <- as.dendrogram(BMP_ligs_log_hc)

#Plots
pdf(file="droplet/Marrow/heatmaps/BMP_ligands.pdf")
BMP_hm <- heatmap3(
  x = BMP_ligs_mtx, 
  main = 'Marrow\nCells (n = 2600) vs BMP ligands',
  Rowv=NA,
  Colv=BMP_ligs_dend,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/BMP_ligands_log.pdf")
BMP_hm <- heatmap3(
  x = BMP_ligs_log_mtx, 
  main = 'Marrow\nCells (n = 2600) vs BMP ligands',
  Rowv=NA,
  Colv=BMP_ligs_log_dend,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

# BMP receptors

#Euclidean distance between cols
BMP_recs_dist <- dist(x = t(BMP_recs_mtx),
                      method = "euclidean"
)
BMP_recs_hc <- hclust(d = BMP_recs_dist, 
                      method = "complete"
)
BMP_recs_dend <- as.dendrogram(BMP_recs_hc)

#Euclidean distance with log10(x+1) transformation
BMP_recs_log_dist <- dist(x = t(BMP_recs_log_mtx),
                          method = "euclidean"
)
BMP_recs_log_hc <- hclust(d = BMP_recs_log_dist, 
                          method = "complete"
)
BMP_recs_log_dend <- as.dendrogram(BMP_recs_log_hc)

#Plots
pdf(file="droplet/Marrow/heatmaps/BMP_receptors.pdf")
BMP_hm <- heatmap3(
  x = BMP_recs_mtx, 
  main = 'Marrow\nCells (n = 1588) vs BMP receptors',
  Rowv=NA,
  Colv=BMP_recs_dend,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/BMP_receptors_log.pdf")
BMP_hm <- heatmap3(
  x = BMP_recs_log_mtx, 
  main = 'Marrow\nCells (n = 1588) vs BMP receptors',
  Rowv=NA,
  Colv=BMP_recs_log_dend,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

# Notch genes
#Euclidean distance between cols
Notch_dist_euc <- dist(x = t(Notch_mtx),
                       method = "euclidean"
)
Notch_hc_euc <- hclust(d = Notch_dist_euc, 
                       method = "complete"
)
Notch_dend_euc <- as.dendrogram(Notch_hc_euc)

#Euclidean distance with log10(x+1) transformation
Notch_log_dist_euc <- dist(x = t(Notch_log_mtx),
                           method = "euclidean"
)
Notch_log_hc_euc <- hclust(d = Notch_log_dist_euc, 
                           method = "complete"
)
Notch_log_dend_euc <- as.dendrogram(Notch_log_hc_euc)

#Manhattan distance between cols 
Notch_dist_man <- dist(x = t(Notch_mtx),
                       method = "manhattan"
)
Notch_hc_man <- hclust(d = Notch_dist_man, 
                       method = "complete"
)
Notch_dend_man <- as.dendrogram(Notch_hc_man)

#Manhattan distance with log10(x+1) transformation
Notch_log_dist_man <- dist(x = t(Notch_log_mtx),
                           method = "manhattan"
)
Notch_log_hc_man <- hclust(d = Notch_log_dist_man, 
                           method = "complete"
)
Notch_log_dend_man <- as.dendrogram(Notch_log_hc_man)

#Correlation 
Notch_cor <- cor(t(Notch_mtx))
Notch_cor_hc <- hclust(d = as.dist(Notch_cor), method = 'complete')
Notch_cor_dend <- as.dendrogram(Notch_cor_hc)

#Plots
pdf(file="droplet/Marrow/heatmaps/Notch_euclidean.pdf")
Notch_hm <- heatmap3(
  x = Notch_mtx, 
  main = 'Marrow Notch (Euclidean)',
  # Rowv=NA,
  Colv=Notch_dend_euc,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/Notch_log10_euclidean.pdf")
Notch_hm <- heatmap3(
  x = log10(Notch_log_mtx + 1), 
  main = 'Marrow Notch (Log10-Euclidean)',
  # Rowv=NA,
  Colv=Notch_log_dend_euc,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/Notch_manhattan.pdf")
Notch_hm <- heatmap3(
  x = Notch_mtx, 
  main = 'Marrow Notch (Manhattan)',
  # Rowv=NA,
  Colv=Notch_dend_man,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/Notch_log10_manhattan.pdf")
Notch_hm <- heatmap3(
  x = log10(Notch_log_mtx + 1), 
  main = 'Marrow Notch (Log10-Manhattan)',
  # Rowv=NA,
  Colv=Notch_log_dend_man,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/Notch_correlation.pdf")
Notch_hm <- heatmap3(
  x = Notch_cor, 
  main = 'Marrow Notch (Correlation)',
  distfun = NA,
  scale = "none",
  showRowDendro = F,
  # Rowv=NA,
  Colv=Notch_cor_dend,
  Rowv=Notch_cor_dend,
  # labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

# Most important Notch genes

#Euclidean distance between cols
Notch_imp_dist <- dist(x = t(Notch_imp_mtx),
                      method = "euclidean"
)
Notch_imp_hc <- hclust(d = Notch_imp_dist, 
                      method = "complete"
)
Notch_imp_dend <- as.dendrogram(Notch_imp_hc)

#Euclidean distance with log10(x+1) transformation
Notch_imp_log_dist <- dist(x = t(Notch_imp_log_mtx),
                          method = "euclidean"
)
Notch_imp_log_hc <- hclust(d = Notch_imp_log_dist, 
                          method = "complete"
)
Notch_imp_log_dend <- as.dendrogram(Notch_imp_log_hc)

#Plots
pdf(file="droplet/Marrow/heatmaps/Notch_imp_genes.pdf")
BMP_hm <- heatmap3(
  x = Notch_imp_mtx, 
  main = 'Marrow\nCells (n = 2160) vs Notch genes',
  Rowv=NA,
  Colv=Notch_imp_dend,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/Marrow/heatmaps/Notch_imp_genes_log.pdf")
BMP_hm <- heatmap3(
  x = Notch_imp_log_mtx, 
  main = 'Marrow\nCells (n = 2160) vs Notch genes',
  Rowv=NA,
  Colv=Notch_imp_log_dend,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()








