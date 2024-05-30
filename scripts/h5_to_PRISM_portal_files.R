library(tidyverse)
library(useful)
library(magrittr)
library(rhdf5)


read_dataset <- function(file = "data/depmap_datasets.h5", dataset, rownames_depmap_ids = TRUE) {
  X <- h5read(file, name = paste0(dataset, "/mat"))
  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"))
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"))
  colnames(X) <- column_meta$column_name
  if(rownames_depmap_ids){
    rownames(X) <- row_meta$ModelID
  }else{
    rownames(X) <- row_meta$CCLEName
  }
  X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
  X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE] 
  return(X)
}



datasets <- h5ls("data/depmap_datasets.h5") 
datasets <- unique(substr(setdiff(datasets$group, "/"),2,100))


# write continuous features -----


dataset = "GeneEffect"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/xpr.csv")

dataset = "Expression"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/ge.csv")

dataset = "ExpressionGeneSetEnrichment"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/gse.csv")

dataset = "Repurposing"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/rep.csv")


datasets
dataset = "RNAi"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/shrna.csv")


dataset = "RPPA"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/rppa.csv")

dataset = "GDSC.Proteomics"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/gdsc_prot.csv")


dataset = "CopyNumber"
X <- read_dataset("data/depmap_datasets.h5", dataset, rownames_depmap_ids = FALSE)
write.csv(X, "data/cna.csv")


# write discrete features ----

# TBD

# write random forest matrices ----

# TBD

