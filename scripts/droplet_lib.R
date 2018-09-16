## ------------------------------------------------------------------------
# Load a package. Install it first if it hasn't been installed.

usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

## ------------------------------------------------------------------------
# Read .mtx files in a given directory into the environment. 
# Parameter 'pattern' allows you to search for files (using regular  
# expressions by default or explicit matching if fixed = TRUE)

readAllMtx <- function(path, pattern = '\\.mtx$', fixed = FALSE) {
  
  #Load required packages
  if (!is.element('Matrix', installed.packages()[,1])) {
    install.packages(p, dep = TRUE)
  }
  require('Matrix', character.only = TRUE)
  
  # Find all matching filenames
  if (fixed) {
    files <- list.files(path = path)
    files <- grep(pattern = pattern, x = files, value = TRUE, fixed = TRUE)
  } else {
    files <- list.files(path = path, pattern = pattern)
  }
  if (substr(path, nchar(path), nchar(path)) != '/') path = paste0(path, '/')
  files <- paste0(path, files)
  
  # Make a list and fill it with the matrix objects from readMM()
  mtx_list <- sapply(X = files, FUN = readMM)
  
  # Return output
  return(mtx_list)
}

