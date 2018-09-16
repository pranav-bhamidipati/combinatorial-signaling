source('~/Caltech/Thomson_Summer/scripts/droplet_lib.R')
setwd('~/Caltech/Thomson_Summer/')

## ------------------------------------------------------------------------
# Define lists of signaling genes
#BMP
BMP_genes <- "BMP2
BMP3
BMP4
BMP5
BMP6
BMP7
BMP8a
GDF2
BMP10
BMP15
GDF1
GDF3
GDF5
GDF6
GDF7
MSTN
GDF9
GDF10
GDF11
GDF15
Nodal
INHBA
INHBB
INHBC
INHBE
INHA
TGFb1
TGFb2
TGFb3
Smad1
Smad2
Smad3
Smad5
Smad9
Smad4
Smad6
Smad7
ACVRL1
ACVR1
BMPR1A
ACVR1B
TGFBR1
BMPR1B
ACVR1C
ACVR2A
ACVR2b
BMPR2
TGFBR2
TGFBR3
ENG
CFC1
RGMa
RGMB
HFE2
Bambi
CHRD
FST
LEFTY1
NOG
CHRDL1
CHRDL2
GREM1
GREM2
CER1
NBL1
LEFTY2
Id1"
BMP_genes <- strsplit(
  x = BMP_genes,
  split = '\n'
)
BMP_genes <- BMP_genes[[1]]

#Notch
Notch_genes <- "DLL1
DLL3
DLL4
DTX1
JAG1
JAG2
ADAM10
PSEN1
PSEN2
PEN2
NOTCH1
NOTCH2
NOTCH3
NOTCH4
HES1
HES5
HES6
HEY1
HEYL
CDKN1A
CFLAR
IL2RA
NFKB1
CCND1
ERBB2
FOSL1
PPARG
FOS
NFKB2
NR4A2
STAT6
CD44
CHUK
IFNG
IL17B
KRT1
LOR
MAP2K7
PDPK1
PTCRA
MFNG
RFNG
LFNG
DLK1
DLK2
NOTCH2NL"
Notch_genes <- strsplit(
  x = Notch_genes,
  split = '\n'
)
Notch_genes <- Notch_genes[[1]]
# ***Remove NOTCH2NL for mouse data because it's a human-specific gene!***
Notch_genes <- Notch_genes[-length(Notch_genes)]

#Wnt
Wnt_genes <- "Myc
CCND1
TCF1
PPARD
MMP7
Axin2
CD44
CTNNB1
Wnt
LRP5
LRP6
Frizzled
TAK1
NLK
TAZ
Snail1
FoxO
LEF1"
Wnt_genes <- strsplit(
  x = Wnt_genes,
  split = '\n'
)
Wnt_genes <- Wnt_genes[[1]]

#Shh
Shh_genes <- "CCND1
CCNE2
Myc
HHIP
FoxM1
Ptc
Smo
Gpc3
Gpc4
Gpc6
Cdo
Boc
Gas1
Kif7
Sufu
PKA
CK1
GSK3
bTRCP
Shh
Dhh
Ihh
Gli1
Gli2
Gli3"
Shh_genes <- strsplit(
  x = Shh_genes,
  split = '\n'
)
Shh_genes <- Shh_genes[[1]]

#Make a list of all of them
sig_genes <- list(
  BMP_genes, 
  Notch_genes, 
  Wnt_genes, 
  Shh_genes
  )
#Make lowercase for easy string matching
sig_genes <- sapply(sig_genes, tolower)

# Save lists of signaling genes
write.table(x = sig_genes[[1]], file = './BMP_genes.txt', row.names = F, col.names = F)
write.table(x = sig_genes[[2]], file = './Notch_genes.txt', row.names = F, col.names = F)
write.table(x = sig_genes[[3]], file = './Wnt_genes.txt', row.names = F, col.names = F)
write.table(x = sig_genes[[4]], file = './Shh_genes.txt', row.names = F, col.names = F)
write.table(x = sort(unique(unlist(sig_genes))), 
            file = './all_signaling_genes.txt', row.names = F, col.names = F)

## ------------------------------------------------------------------------
# Load the list of genes
genes <- read.table(
  file = "C:/Users/Pranav/Documents/Caltech/Thomson_Summer/droplet/genes.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
genes <- genes$V1
#Make all chars lowercase for string matching
genes <- sapply(genes, tolower)

## ------------------------------------------------------------------------
#Search for each gene name in the genes list and store index in another list
#No match is stored as index -1
sig_genes_idx <- sapply(sig_genes, 
                        match, 
                        table=genes, 
                        nomatch=-1)

# #If this worked perfectly, there should be zero no-matches, and there should be 3 numbers repeated (the 3 multiply listed genes)
table(unlist(sig_genes_idx) > 0)    #15 no-matches

#Identify which genes were not found in the genes list
unlist(sig_genes)[unlist(sig_genes_idx) < 0]

## ------------------------------------------------------------------------
#Remove unmapped genes from signaling genes list
sig_genes <- mapply(
  FUN = function(genes,indices) {
    genes[indices>0]
  },
  sig_genes,
  sig_genes_idx
)
sig_genes_idx <- sapply(
  sig_genes_idx,
  function(x) x[which(x>0)]
)
