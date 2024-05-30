library(tidyverse)
library(useful)
library(magrittr)
library(rhdf5)
library(taigr)

# make sure to change the taigaclient path below.
options(taigaclient.path=path.expand("/opt/miniconda3/envs/taigapy/bin/taigaclient"))

# release path
release = 'public-24q2-356f'
version = '42'


# RELEASE FILES ------

# THIS STILL NEDS TO BE UPDATED !!! 
map <- load.from.taiga(data.name= release , data.version = version, data.file = 'Model') %>% 
  dplyr::distinct(CCLEName, ModelID)

# Expression
OmicsExpressionGeneSetEnrichment <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsExpressionGeneSetEnrichment')
OmicsExpressionProteinCodingGenesTPMLogp1 <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsExpressionProteinCodingGenesTPMLogp1')
prc_ids <- word(colnames(OmicsExpressionProteinCodingGenesTPMLogp1), -1)
gene_ids <- word(colnames(OmicsExpressionProteinCodingGenesTPMLogp1))
colnames(OmicsExpressionProteinCodingGenesTPMLogp1) %<>% word() 

# Copy Number
# OmicsArmLevelCNA <- load.from.taiga(data.name='internal-24q2-3719', data.version=63, data.file='OmicsArmLevelCNA')
OmicsCNGene <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsCNGene') 
OmicsCNGene <- log2(OmicsCNGene + 1)
flag <- (word(colnames(OmicsCNGene), -1) %in% prc_ids) | (word(colnames(OmicsCNGene)) %in% gene_ids)
OmicsCNGene <- OmicsCNGene[, flag]
colnames(OmicsCNGene) %<>% word()

# Fusion
OmicsFusionFiltered <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsFusionFiltered') %>% 
  dplyr::distinct(ModelID, FusionName,
                  LeftGene, RightGene) %>% 
  dplyr::filter((word(LeftGene) %in% gene_ids) | (word(RightGene) %in% gene_ids)) %>% 
  dplyr::distinct(ModelID, FusionName) %>% 
  dplyr::mutate(dummy = 1) %>%
  reshape2::acast(ModelID ~ FusionName, value.var= "dummy", fill = 0) 

# LoH
OmicsLoH <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsLoH')
flag <- (word(colnames(OmicsLoH), -1) %in% prc_ids) | (word(colnames(OmicsLoH)) %in% gene_ids)
OmicsLoH <- OmicsLoH[, flag]
colnames(OmicsLoH) %<>% word()

# Signatures
OmicsSignatures <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsSignatures')

# Mutations
OmicsSomaticMutationsMatrixDamaging <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsSomaticMutationsMatrixDamaging')
OmicsSomaticMutationsMatrixHotspot <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsSomaticMutationsMatrixHotspot')
flag <- (word(colnames(OmicsSomaticMutationsMatrixHotspot), -1) %in% prc_ids) | (word(colnames(OmicsSomaticMutationsMatrixHotspot)) %in% gene_ids)
OmicsSomaticMutationsMatrixHotspot <- OmicsSomaticMutationsMatrixHotspot[, flag]
flag <- (word(colnames(OmicsSomaticMutationsMatrixDamaging), -1) %in% prc_ids) | (word(colnames(OmicsSomaticMutationsMatrixDamaging)) %in% gene_ids)
OmicsSomaticMutationsMatrixDamaging <- OmicsSomaticMutationsMatrixDamaging[, flag]
colnames(OmicsSomaticMutationsMatrixDamaging) %<>% word()
colnames(OmicsSomaticMutationsMatrixHotspot) %<>% word()
OmicsSomaticMutationsMatrixDamaging <- (OmicsSomaticMutationsMatrixDamaging > 0 ) + 0
OmicsSomaticMutationsMatrixHotspot <- (OmicsSomaticMutationsMatrixHotspot > 0 ) + 0
OmicsSomaticMutations <- load.from.taiga(data.name= release, data.version=version, data.file='OmicsSomaticMutations') %>% 
  dplyr::filter(Hotspot) %>% 
  dplyr::distinct(ModelID, HugoSymbol, ProteinChange) %>% 
  dplyr::filter((HugoSymbol %in% gene_ids), !is.na(ProteinChange)) %>%
  dplyr::mutate(cn = paste0(HugoSymbol, ".", ProteinChange),
                dummy = 1) %>% 
  reshape2::acast(ModelID ~ cn, value.var = "dummy", fill = 0) 
  


# # PRISM OncRef 
# PRISMOncologyReferenceAUCMatrix <- load.from.taiga(data.name='prism-oncology-reference-set-24q2-f31b', data.version=4, data.file='PRISMOncologyReferenceAUCMatrix')
# PRISMOncologyReferenceCompoundList <- load.from.taiga(data.name='prism-oncology-reference-set-24q2-f31b', data.version=4, data.file='PRISMOncologyReferenceCompoundList')
# colnames(PRISMOncologyReferenceAUCMatrix) <- tibble( SampleID = colnames(PRISMOncologyReferenceAUCMatrix))  %>% 
#   dplyr::left_join(dplyr::distinct(PRISMOncologyReferenceCompoundList, SampleID, CompoundName)) %>%
#   .$CompoundName

# CRISPR 
CRISPRGeneEffect <- load.from.taiga(data.name= release, data.version=version, data.file='CRISPRGeneEffect')
flag = (word(colnames(CRISPRGeneEffect)) %in% gene_ids) | (word(colnames(CRISPRGeneEffect), -1) %in% prc_ids)
CRISPRGeneEffect <- CRISPRGeneEffect[, flag]
colnames(CRISPRGeneEffect) %<>% word()


# LEGACY DATASETS  ------


# PRISM Repurposing
Repurposing.23Q2.Extended.Primary.Compound.List <- load.from.taiga(data.name='repurposing-23q2-a803', data.version=3, data.file='Repurposing_23Q2_Extended_Primary_Compound_List')
Repurposing.23Q2.Extended.Primary.Data.Matrix <- load.from.taiga(data.name='repurposing-23q2-a803', data.version=3, data.file='Repurposing_23Q2_Extended_Primary_Data_Matrix')
REPURPOSING.PRIMARY <- Repurposing.23Q2.Extended.Primary.Compound.List %>% 
  dplyr::inner_join(tibble(IDs = rownames(Repurposing.23Q2.Extended.Primary.Data.Matrix),
                          n = apply(Repurposing.23Q2.Extended.Primary.Data.Matrix, 1, function(x) sum(is.finite(x))),
                          ns = apply(Repurposing.23Q2.Extended.Primary.Data.Matrix, 1, function(x) sum(x < log2(.25), na.rm = T)))) %>% 
  dplyr::group_by(Drug.Name) %>% 
  dplyr::top_n(1, n) %>%
  dplyr::top_n(1, ns) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(Repurposing.23Q2.Extended.Primary.Data.Matrix %>% 
                     reshape2::melt() %>% 
                     dplyr::rename(IDs = Var1)) %>%
  reshape2::acast(Var2 ~ Drug.Name, value.var = "value",
                  fun.aggregate = function(x) median(x, na.rm = T))
  
# RNAi
RNAi <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=13, data.file='gene_effect')
flag = (word(colnames(RNAi)) %in% gene_ids) | (word(colnames(RNAi), -1) %in% prc_ids)
RNAi <- RNAi[, flag]
colnames(RNAi) %<>% word()
rownames(RNAi) <- tibble(CCLEName = rownames(RNAi)) %>%
  dplyr::left_join(map) %>%
  .$ModelID

# RPPA
CCLE.RPPA.20180123 <- load.from.taiga(data.name='depmap-rppa-1b43', data.version=1, data.file='CCLE_RPPA_20180123')
CCLE.RPPA.Ab.info.20180123 <- load.from.taiga(data.name='depmap-rppa-1b43', data.version=1, data.file='CCLE_RPPA_Ab_info_20180123')
colnames(CCLE.RPPA.20180123) <- tibble(Antibody_Name = colnames(CCLE.RPPA.20180123)) %>% 
  dplyr::left_join(CCLE.RPPA.Ab.info.20180123 %>% 
                     dplyr::distinct(Antibody_Name, Target_Genes)) %>% 
  dplyr::mutate(cn = paste0(Antibody_Name, "_", Target_Genes)) %>% 
  .$cn
rownames(CCLE.RPPA.20180123) <- tibble(CCLEName = rownames(CCLE.RPPA.20180123)) %>%
  dplyr::left_join(map) %>%
  .$ModelID

# # SangerMethylation
# MethylationExpressionImpactSanger <- load.from.taiga(data.name='internal-beta-features-7130', data.version=31, data.file='MethylationExpressionImpactSanger')
# flag <- (word(colnames(MethylationExpressionImpactSanger), -1) %in% prc_ids) | (word(colnames(MethylationExpressionImpactSanger)) %in% gene_ids)
# MethylationExpressionImpactSanger <- MethylationExpressionImpactSanger[, flag]
# colnames(MethylationExpressionImpactSanger) %<>% word()

# SangerProteomics
sangerProt <- load.from.taiga(data.name='sanger-proteomics-be1c', data.version=2, data.file='sangerProt')



# WRITE THE MATRICES INTO A H5 FILE---------


h5createFile("data/depmap_datasets.h5")


X = OmicsExpressionGeneSetEnrichment
name = "ExpressionGeneSetEnrichment"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))


X = OmicsExpressionProteinCodingGenesTPMLogp1
name = "Expression"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))



X = OmicsCNGene
name = "CopyNumber"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))


X = OmicsFusionFiltered
name = "Fusion"
X <- X[, colSums(X, na.rm = T) > 2 ]
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))


X = OmicsLoH
name = "LoH"
X <- X[, colSums(X, na.rm = T) > 2 ]
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))

X = OmicsSignatures
name = "Signatures"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))

X1 <- OmicsSomaticMutationsMatrixDamaging
colnames(X1) %<>% paste0(., ".DAM")
X2 <- OmicsSomaticMutationsMatrixHotspot 
colnames(X2) %<>% paste0(., ".HS")
X3 <- OmicsSomaticMutations 
X <- reshape2::melt(X1) %>% 
  dplyr::filter(value == 1) %>% 
  dplyr::bind_rows(reshape2::melt(X2) %>% 
                     dplyr::filter(value == 1)) %>%
  dplyr::bind_rows(reshape2::melt(X3) %>% 
                     dplyr::filter(value == 1)) %>%
  reshape2::acast(Var1 ~ Var2, value.var = "value", fill = 0) 
rm(X1, X2, X3)
X <- X[,colSums(X, na.rm = T) > 2] 
name = "Mutation"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))


X <- CRISPRGeneEffect
name = "GeneEffect"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))



X <- REPURPOSING.PRIMARY
name = "Repurposing"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))


X <- RNAi 
X <- X[!is.na(rownames(X)),]
name = "RNAi"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))



X <- CCLE.RPPA.20180123
name = "RPPA"
X <- X[!is.na(rownames(X)),]
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))



X <- sangerProt
name = "GDSC.Proteomics"
!c(any(duplicated(colnames(X))) , any(is.na(colnames(X))) , any(duplicated(rownames(X))) , any(is.na(rownames(X))))
h5createGroup("data/depmap_datasets.h5", name)
row_meta <- tibble(ix = 1:dim(X)[1], ModelID = rownames(X)) %>% 
  dplyr::left_join(map)
column_meta <- tibble(ix = 1:dim(X)[2], column_name = colnames(X))   
h5write(row_meta, "data/depmap_datasets.h5", paste0(name, "/row_meta"))
h5write(column_meta, "data/depmap_datasets.h5",  paste0(name, "/column_meta"))
h5createDataset("data/depmap_datasets.h5",  paste0(name, "/mat"), dim(X),
                storage.mode = "double", chunk=c(dim(X)[1],1), level=9)
h5write(X, file= "data/depmap_datasets.h5",
        name= paste0(name, "/mat"), index=list(1:dim(X)[1],1:dim(X)[2]))

h5ls("data/depmap_datasets.h5")
