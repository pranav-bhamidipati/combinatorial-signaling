source('~/Caltech/Thomson_Summer/scripts/droplet_lib.R')
setwd('~/Caltech/Thomson_Summer/')

# Load all droplet data sets
dat_filenames <- list.files('droplet/', pattern = "\\.mtx")
dat <- sapply(dat_filenames,
              function(x) readAllMtx('./droplet/', x, fixed = T))
dat <- lapply(dat, function(x) as.data.frame(summary(x)))
names(dat) <- dat_filenames

# #Load four droplet data sets (heart, marrow, muscle, thymus) and make into dfs
# dat_filenames <- list.files('droplet/', pattern = "\\.mtx")[c(4, 17, 19, 23)]
# dat <- sapply(dat_filenames, 
#               function(x) readAllMtx('./droplet/', x, fixed = T))
# dat <- lapply(dat, function(x) as.data.frame(summary(x)))
# names(dat) <- dat_filenames

## ------------------------------------------------------------------------
# Find sum of reads in each cell
dat_nreads <- lapply(
  X = dat,
  FUN = function(df) {
    aggregate(x = df$x, 
              by = list(cell=df$j), 
              FUN=sum)
  })

## ------------------------------------------------------------------------
# Normalize data to # reads per 10^6 transcripts in cell
dat_norm <- mapply(
  function(df, nreads) {
    df$x <- df$x / rep(nreads$x, table(df$j)) * 10^6
    return(df)
  },
  dat,
  dat_nreads,
  SIMPLIFY = F
)

# # Normalize using f(x) = log2(x + 1)
# # NOTE: Matt said we may have to do this, but try without it first

## ------------------------------------------------------------------------
# Write normalized data to files

for (n in 1:length(dat_norm)) {
  fname = names(dat_norm)[n]
  fname = paste0(
    './droplet/normalized/',
    substr(fname, 1, nchar(fname) - 4),
    '_norm.txt')
  print(fname)
  write.table(
    x = dat_norm[[n]],
    file = fname,
    row.names = F
  )
}
