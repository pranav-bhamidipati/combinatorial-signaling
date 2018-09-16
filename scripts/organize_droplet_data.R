# Rename the data matrix files from "matrix.mtx" to their tissue-specific 
# name and put them directly in the droplet data folder.
path = '~/Caltech/Thomson_Summer/'
setwd(dir = path)
mtx_from <- list.files(path = './droplet/', pattern = '\\.mtx$', recursive = T, full.names = T)
mtx_to <- list.files(path = './droplet/', full.names = T)
mtx_to <- paste0(mtx_to,'.mtx')
file.rename(mtx_from, mtx_to)
#### I then manually deleted the folders that used to hold the droplet data, 
#### along with the text files containing the gene lists (all identical) and barcodes