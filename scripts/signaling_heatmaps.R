# Set up environment
setwd('~/Caltech/Thomson_Summer/')
source('./scripts/droplet_lib.R')
usePackage('reshape2')
usePackage('heatmap3')

# Read in droplet data sets (heart, marrow, muscle, thymus)
dat_filenames <- list.files(
  'droplet/normalized/', 
  pattern = "\\.txt")[c(4, 17, 19, 23)]
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

# Melt data into one df
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

# Read in master gene list and signaling gene list
genes <- read.table(file = "droplet/genes.txt", stringsAsFactors = F)
genes <- tolower(genes$V1)
sig_genes <- read.table(file = "all_signaling_genes.txt", stringsAsFactors = F)
sig_genes <- sig_genes$x
sig_genes_idx <- which(genes %in% sig_genes)
BMP_genes <- read.table('./BMP_genes.txt', stringsAsFactors = F)
BMP_genes <- tolower(BMP_genes$V1)
BMP_genes_idx <- which(genes %in% BMP_genes)
Notch_genes <- read.table('./Notch_genes.txt', stringsAsFactors = F)
Notch_genes <- tolower(Notch_genes$V1)
Notch_genes_idx <- which(genes %in% Notch_genes)

# Select samples of genes for study 
all_expr <- dat_df[,c('gene','expr')]
all_expr <- all_expr[order(all_expr[,1]),]
rownames(all_expr) <- c()
colnames(all_expr) <- c('gene','expr')
means = aggregate(x = all_expr$expr, by = list(all_expr$gene), FUN = mean)
colnames(means) <- c('gene','mean_expr')
vars = aggregate(x = all_expr$expr, by = list(all_expr$gene), FUN = var)
colnames(vars) <- c('gene','var_expr')

#Random sample of genes with above-median expression (to avoid excessive zero values)
random_high_pass_sample <- sample(
  means[means$mean_expr > quantile(means$mean_expr, 0.5),1], 
  size = 500, 
  replace = F)
#Top 1000 and 500 genes by variance
top_1000_vars_sample <- vars[order(vars$var_expr, decreasing = T)[1:1000],1]
top_500_vars_sample <- vars[order(vars$var_expr, decreasing = T)[1:500],1]
#Signaling genes
signaling_sample <- sig_genes_idx[sig_genes_idx %in% all_expr$gene]
BMP_sample <- BMP_genes_idx[BMP_genes_idx %in% all_expr$gene]
Notch_sample <- Notch_genes_idx[Notch_genes_idx %in% all_expr$gene]

# Generate data sets based on gene samples
#Top variance
top_1000_vars_dat <- dat_df[dat_df$gene %in% top_1000_vars_sample,]
top_500_vars_dat <- dat_df[dat_df$gene %in% top_500_vars_sample,]
#Random high-pass
rhp_dat <- dat_df[dat_df$gene %in% random_high_pass_sample,]
#Signaling
sig_dat <- dat_df[dat_df$gene %in% signaling_sample,]
BMP_dat <- dat_df[dat_df$gene %in% BMP_sample,]
Notch_dat <- dat_df[dat_df$gene %in% Notch_sample,]

# Initialize matrices without any skipped gene x cell combinations
#Top variance
top_1000_vars_mtx <- data.frame(
  gene=rep(top_1000_vars_sample,
           length(unique(top_1000_vars_dat$tissue.cell))),
  tissue.cell=rep(unique(top_1000_vars_dat$tissue.cell),
           each=length(top_1000_vars_sample)),
  expr=numeric(length(top_1000_vars_sample) * length(unique(top_1000_vars_dat$tissue.cell))),
  stringsAsFactors = F
)

top_500_vars_mtx <- data.frame(
  gene=rep(top_500_vars_sample,
           length(unique(top_500_vars_dat$tissue.cell))),
  tissue.cell=rep(unique(top_500_vars_dat$tissue.cell),
                  each=length(top_500_vars_sample)),
  expr=numeric(length(top_500_vars_sample) * length(unique(top_500_vars_dat$tissue.cell))),
  stringsAsFactors = F
)

#Random high-pass
rhp_mtx <- data.frame(
  gene=rep(random_high_pass_sample,
           length(unique(rhp_dat$tissue.cell))),
  tissue.cell=rep(unique(rhp_dat$tissue.cell),
           each=length(random_high_pass_sample)),
  expr=numeric(length(random_high_pass_sample) * 
                 length(unique(rhp_dat$tissue.cell))),
  stringsAsFactors = F
)

#Signaling
sig_mtx <- data.frame(
  gene=rep(signaling_sample,
           length(unique(sig_dat$tissue.cell))),
  tissue.cell=rep(unique(sig_dat$tissue.cell),
           each=length(signaling_sample)),
  expr=numeric(length(signaling_sample) * length(unique(sig_dat$tissue.cell))),
  stringsAsFactors = F
)

BMP_mtx <- data.frame(
  gene=rep(BMP_sample,
           length(unique(BMP_dat$tissue.cell))),
  tissue.cell=rep(unique(BMP_dat$tissue.cell),
                  each=length(BMP_sample)),
  expr=numeric(length(BMP_sample) * length(unique(BMP_dat$tissue.cell))),
  stringsAsFactors = F
)

Notch_mtx <- data.frame(
  gene=rep(Notch_sample,
           length(unique(Notch_dat$tissue.cell))),
  tissue.cell=rep(unique(Notch_dat$tissue.cell),
                  each=length(Notch_sample)),
  expr=numeric(length(Notch_sample) * length(unique(Notch_dat$tissue.cell))),
  stringsAsFactors = F
)

# Fill matrices with signaling expression values
#Top variance
top_1000_vars_mtx_rows = do.call(paste0, top_1000_vars_mtx[c('gene','tissue.cell')])
top_1000_vars_dat_rows = do.call(paste0, top_1000_vars_dat[c('gene','tissue.cell')])
matching_rows = match(top_1000_vars_dat_rows, top_1000_vars_mtx_rows)
top_1000_vars_mtx$expr[matching_rows] = top_1000_vars_dat$expr

top_1000_vars_mtx = acast(top_1000_vars_mtx, gene ~ tissue.cell, value.var = 'expr')
top_1000_vars_log_mtx = log10(top_1000_vars_mtx+1)

top_500_vars_mtx_rows = do.call(paste0, top_500_vars_mtx[c('gene','tissue.cell')])
top_500_vars_dat_rows = do.call(paste0, top_500_vars_dat[c('gene','tissue.cell')])
matching_rows = match(top_500_vars_dat_rows, top_500_vars_mtx_rows)
top_500_vars_mtx$expr[matching_rows] = top_500_vars_dat$expr

top_500_vars_mtx = acast(top_500_vars_mtx, gene ~ tissue.cell, value.var = 'expr')
top_500_vars_log_mtx = log10(top_500_vars_mtx+1)

#Random high-pass
rhp_mtx_rows = do.call(paste0, rhp_mtx[c('gene','tissue.cell')])
rhp_dat_rows = do.call(paste0, rhp_dat[c('gene','tissue.cell')])
matching_rows = match(rhp_dat_rows, rhp_mtx_rows)
rhp_mtx$expr[matching_rows] = rhp_dat$expr

rhp_mtx = acast(rhp_mtx, gene ~ tissue.cell, value.var = 'expr')
rhp_log_mtx = log10(rhp_mtx+1)

#Signaling
sig_mtx_rows = do.call(paste0, sig_mtx[c('gene','tissue.cell')])
sig_dat_rows = do.call(paste0, sig_dat[c('gene','tissue.cell')])
matching_rows = match(sig_dat_rows, sig_mtx_rows)
sig_mtx$expr[matching_rows] = sig_dat$expr

sig_mtx = acast(sig_mtx, gene ~ tissue.cell, value.var = 'expr')
sig_log_mtx = log10(sig_mtx+1)

#BMP
BMP_mtx_rows = do.call(paste0, BMP_mtx[c('gene','tissue.cell')])
BMP_dat_rows = do.call(paste0, BMP_dat[c('gene','tissue.cell')])
matching_rows = match(BMP_dat_rows, BMP_mtx_rows)
BMP_mtx$expr[matching_rows] = BMP_dat$expr

BMP_mtx = acast(BMP_mtx, gene ~ tissue.cell, value.var = 'expr')
BMP_log_mtx = log10(BMP_mtx+1)

#Notch
Notch_mtx_rows = do.call(paste0, Notch_mtx[c('gene','tissue.cell')])
Notch_dat_rows = do.call(paste0, Notch_dat[c('gene','tissue.cell')])
matching_rows = match(Notch_dat_rows, Notch_mtx_rows)
Notch_mtx$expr[matching_rows] = Notch_dat$expr

Notch_mtx = acast(Notch_mtx, gene ~ tissue.cell, value.var = 'expr')
Notch_log_mtx = log10(Notch_mtx+1)

## Clustering and heatmaps
# Top 1000 (too big, kept crashing)
# Top 500 
#Euclidean distance between cols
top_500_dist_euc <- dist(x = t(top_500_vars_mtx),
                         method = "euclidean"
)
top_500_hc_euc <- hclust(d = top_500_dist_euc, 
                           method = "complete"
)
top_500_dend_euc <- as.dendrogram(top_500_hc_euc)

#Euclidean distance with log10(x+1) transformation
top_500_log_dist_euc <- dist(x = t(top_500_vars_log_mtx),
                         method = "euclidean"
)
top_500_log_hc_euc <- hclust(d = top_500_log_dist_euc, 
                         method = "complete"
)
top_500_log_dend_euc <- as.dendrogram(top_500_log_hc_euc)

#Manhattan distance between cols 
top_500_dist_man <- dist(x = t(top_500_vars_mtx),
                     method = "manhattan"
                     )
top_500_hc_man <- hclust(d = top_500_dist_man, 
                         method = "complete"
                         )
top_500_dend_man <- as.dendrogram(top_500_hc_man)

#Manhattan distance with log10(x+1) transformation
top_500_log_dist_man <- dist(x = t(top_500_vars_log_mtx),
                             method = "manhattan"
)
top_500_log_hc_man <- hclust(d = top_500_log_dist_man, 
                             method = "complete"
)
top_500_log_dend_man <- as.dendrogram(top_500_log_hc_man)

#Plots
pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/top_500_euclidean.pdf")
top_500_vars_hm <- heatmap3(
  x = top_500_vars_mtx, 
  Rowv=NA,
  Colv=top_500_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/top_500_log10_euclidean.pdf")
top_500_vars_hm <- heatmap3(
  x = log10(top_500_vars_log_mtx + 1), 
  Rowv=NA,
  Colv=top_500_log_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/top_500_manhattan.pdf")
top_500_vars_hm <- heatmap3(
  x = top_500_vars_mtx, 
  Rowv=NA,
  Colv=top_500_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/top_500_log10_manhattan.pdf")
top_500_vars_hm <- heatmap3(
  x = log10(top_500_vars_log_mtx + 1), 
  Rowv=NA,
  Colv=top_500_log_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

# Random high-pass 
#Euclidean distance between cols
rhp_dist_euc <- dist(x = t(rhp_mtx),
                         method = "euclidean"
)
rhp_hc_euc <- hclust(d = rhp_dist_euc, 
                         method = "complete"
)
rhp_dend_euc <- as.dendrogram(rhp_hc_euc)

#Euclidean distance with log10(x+1) transformation
rhp_log_dist_euc <- dist(x = t(rhp_log_mtx),
                             method = "euclidean"
)
rhp_log_hc_euc <- hclust(d = rhp_log_dist_euc, 
                             method = "complete"
)
rhp_log_dend_euc <- as.dendrogram(rhp_log_hc_euc)

#Manhattan distance between cols 
rhp_dist_man <- dist(x = t(rhp_mtx),
                         method = "manhattan"
)
rhp_hc_man <- hclust(d = rhp_dist_man, 
                         method = "complete"
)
rhp_dend_man <- as.dendrogram(rhp_hc_man)

#Manhattan distance with log10(x+1) transformation
rhp_log_dist_man <- dist(x = t(rhp_log_mtx),
                             method = "manhattan"
)
rhp_log_hc_man <- hclust(d = rhp_log_dist_man, 
                             method = "complete"
)
rhp_log_dend_man <- as.dendrogram(rhp_log_hc_man)

#Plots
pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/rhp_euclidean.pdf")
rhp_hm <- heatmap3(
  x = rhp_mtx, 
  Rowv=NA,
  Colv=rhp_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/rhp_log10_euclidean.pdf")
rhp_hm <- heatmap3(
  x = log10(rhp_log_mtx + 1), 
  Rowv=NA,
  Colv=rhp_log_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/rhp_manhattan.pdf")
rhp_hm <- heatmap3(
  x = rhp_mtx, 
  Rowv=NA,
  Colv=rhp_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/rhp_log10_manhattan.pdf")
rhp_hm <- heatmap3(
  x = log10(rhp_log_mtx + 1), 
  Rowv=NA,
  Colv=rhp_log_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

# Signaling genes
#Euclidean distance between cols
sig_dist_euc <- dist(x = t(sig_mtx),
                     method = "euclidean"
)
sig_hc_euc <- hclust(d = sig_dist_euc, 
                     method = "complete"
)
sig_dend_euc <- as.dendrogram(sig_hc_euc)

#Euclidean distance with log10(x+1) transformation
sig_log_dist_euc <- dist(x = t(sig_log_mtx),
                         method = "euclidean"
)
sig_log_hc_euc <- hclust(d = sig_log_dist_euc, 
                         method = "complete"
)
sig_log_dend_euc <- as.dendrogram(sig_log_hc_euc)

#Manhattan distance between cols 
sig_dist_man <- dist(x = t(sig_mtx),
                     method = "manhattan"
)
sig_hc_man <- hclust(d = sig_dist_man, 
                     method = "complete"
)
sig_dend_man <- as.dendrogram(sig_hc_man)

#Manhattan distance with log10(x+1) transformation
sig_log_dist_man <- dist(x = t(sig_log_mtx),
                         method = "manhattan"
)
sig_log_hc_man <- hclust(d = sig_log_dist_man, 
                         method = "complete"
)
sig_log_dend_man <- as.dendrogram(sig_log_hc_man)

#Plots
pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/sig_euclidean.pdf")
sig_hm <- heatmap3(
  x = sig_mtx, 
  Rowv=NA,
  Colv=sig_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/sig_log10_euclidean.pdf")
sig_hm <- heatmap3(
  x = log10(sig_log_mtx + 1), 
  Rowv=NA,
  Colv=sig_log_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/sig_manhattan.pdf")
sig_hm <- heatmap3(
  x = sig_mtx, 
  Rowv=NA,
  Colv=sig_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/sig_log10_manhattan.pdf")
sig_hm <- heatmap3(
  x = log10(sig_log_mtx + 1), 
  Rowv=NA,
  Colv=sig_log_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

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

#Plots
pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/BMP_euclidean.pdf")
BMP_hm <- heatmap3(
  x = BMP_mtx, 
  Rowv=NA,
  Colv=BMP_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/BMP_log10_euclidean.pdf")
BMP_hm <- heatmap3(
  x = log10(BMP_log_mtx + 1), 
  Rowv=NA,
  Colv=BMP_log_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/BMP_manhattan.pdf")
BMP_hm <- heatmap3(
  x = BMP_mtx, 
  Rowv=NA,
  Colv=BMP_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/BMP_log10_manhattan.pdf")
BMP_hm <- heatmap3(
  x = log10(BMP_log_mtx + 1), 
  Rowv=NA,
  Colv=BMP_log_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
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

#Plots
pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/Notch_euclidean.pdf")
Notch_hm <- heatmap3(
  x = Notch_mtx, 
  Rowv=NA,
  Colv=Notch_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/Notch_log10_euclidean.pdf")
Notch_hm <- heatmap3(
  x = log10(Notch_log_mtx + 1), 
  Rowv=NA,
  Colv=Notch_log_dend_euc,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/Notch_manhattan.pdf")
Notch_hm <- heatmap3(
  x = Notch_mtx, 
  Rowv=NA,
  Colv=Notch_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/Notch_log10_manhattan.pdf")
Notch_hm <- heatmap3(
  x = log10(Notch_log_mtx + 1), 
  Rowv=NA,
  Colv=Notch_log_dend_man,
  labCol = NA,
  # col = cm.colors(256), 
  useRaster = T
)
dev.off()

## Plots with altered colors
#BMP

my_palette <- c(colorRampPalette(c("white", "blue", "purple", "black"))(n = 299))

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/BMP_euclidean_colors.pdf")
BMP_hm <- heatmap3(
  x = BMP_mtx, 
  Rowv=NA,
  Colv=BMP_dend_euc,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/BMP_manhattan_colors.pdf")
BMP_hm <- heatmap3(
  x = BMP_mtx, 
  Rowv=NA,
  Colv=BMP_dend_man,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

#Notch
pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/Notch_euclidean_colors.pdf")
Notch_hm <- heatmap3(
  x = Notch_mtx, 
  Rowv=NA,
  Colv=Notch_dend_euc,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()

pdf(file="droplet/HeartMarrowMuscleSpleen/heatmaps/Notch_manhattan_colors.pdf")
Notch_hm <- heatmap3(
  x = Notch_mtx, 
  Rowv=NA,
  Colv=Notch_dend_man,
  labCol = NA,
  col = my_palette,
  useRaster = T
)
dev.off()



