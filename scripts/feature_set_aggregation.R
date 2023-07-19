library(tidyverse)
library(taigr)
library(magrittr)
library(useful)

Expression <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='OmicsExpressionProteinCodingGenesTPMLogp1')
colnames(Expression) %<>% word()

Expression <- Expression[, apply(Expression, 2, function(x) var(x, na.rm = T)) > .1]

CRISPRGeneEffect <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='CRISPRGeneEffect')
colnames(CRISPRGeneEffect) %<>% word()
CRISPRGeneEffect <- CRISPRGeneEffect[, apply(CRISPRGeneEffect, 2,  function(x) var(x, na.rm = T)) > 0.01]


MatrixDamaging <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='OmicsSomaticMutationsMatrixDamaging')
MatrixHotspot <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='OmicsSomaticMutationsMatrixHotspot')
colnames(MatrixDamaging) %<>% word()
colnames(MatrixHotspot) %<>% word()
MUT <- MatrixDamaging %>%
  reshape2::melt() %>%
  dplyr::filter(value > 0) %>% 
  dplyr::mutate(Var2 = paste0(Var2, "_", "Dmg")) %>% 
  dplyr::bind_rows(  
    MatrixHotspot %>%
      reshape2::melt() %>%
      dplyr::filter(value > 0) %>% 
      dplyr::mutate(Var2 = paste0(Var2, "_", "Hs"))) %>%
  dplyr::mutate(value = 1) %>% 
  reshape2::acast(Var1 ~ Var2 , value.var  ="value", fill = 0)
MUT <- MUT[,apply(MUT, 2,  function(x) sum(x, na.rm = T)) > 9]


compound_list <- data.table::fread("~/Desktop/PREP/RELEASE FILES_23Q2/FINAL/Repurposing_23Q2_Extended_Primary_Compound_List.csv")
compound_matrix <- data.table::fread("~/Desktop/PREP/RELEASE FILES_23Q2/FINAL/Repurposing_23Q2_Extended_Primary_Data_Matrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix()
rownames(compound_matrix) <- tibble(IDs = rownames(compound_matrix) ) %>%
  dplyr::left_join(compound_list %>%
                     dplyr::distinct(IDs, Drug.Name,repurposing_target) %>% 
                     dplyr::group_by(IDs) %>% 
                     dplyr::top_n(1, Drug.Name)) %>%
  .$Drug.Name
compound_matrix %<>% t()
compound_matrix <- compound_matrix[,apply(compound_matrix, 2,  function(x) sum(x < log2(.25), na.rm = T)) > 4]

compound_matrix %<>% 
  reshape2::melt() %>%
  dplyr::filter(is.finite(value)) %>% 
  reshape2::acast(Var1 ~ Var2, value.var = "value", fun.aggregate = median) 


OmicsSignatures <- load.from.taiga(data.name='internal-beta-features-7130', data.version=31, data.file='OmicsSignatures')

DriverGeneFunctionalAlterations <- load.from.taiga(data.name='internal-beta-features-7130', data.version=31, data.file='DriverGeneFunctionalAlterations')
DriverGeneFunctionalAlterations %<>%
  tidyr::pivot_longer(2:1016, names_to = "cn", values_to = "value") %>% 
  dplyr::filter(value %in% c("complete_loss", "enhanced/new_function", "partial_loss")) %>% 
  dplyr::mutate(dummy = 1) %>%
  reshape2::acast(DepMap_ID ~ cn + value, value.var = "dummy", fill = 0) 
DriverGeneFunctionalAlterations <- DriverGeneFunctionalAlterations[, apply(DriverGeneFunctionalAlterations, 2, function(x) sum(x, na.rm = T)) > 4]


Model <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='Model') %>% 
  dplyr::mutate(dummy = 1) %>% 
  reshape2::acast(ModelID ~ DepmapModelType, value.var = "dummy", fill = 0) 
Model <- Model[, apply(Model, 2, function(x) sum(x, na.rm = T)) > 9] 


OmicsCNGene <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='OmicsCNGene')
colnames(OmicsCNGene) %<>% word()
OmicsCNGene <- OmicsCNGene[, apply(OmicsCNGene, 2, function(x) var(x, na.rm = T)) > 0.05] 

MethylationExpressionImpactBroad <- load.from.taiga(data.name='internal-beta-features-7130', data.version=31, data.file='MethylationExpressionImpactBroad')
colnames(MethylationExpressionImpactBroad) %<>% word()
MethylationExpressionImpactBroad <- MethylationExpressionImpactBroad[, apply(MethylationExpressionImpactBroad, 2, function(x) var(x, na.rm = T)) > 0.01] 


MethylationExpressionImpactSanger <- load.from.taiga(data.name='internal-beta-features-7130', data.version=31, data.file='MethylationExpressionImpactSanger')
colnames(MethylationExpressionImpactSanger) %<>% word()
MethylationExpressionImpactSanger <- MethylationExpressionImpactSanger[, apply(MethylationExpressionImpactSanger, 2, function(x) var(x, na.rm = T)) > 0.01] 


webster.dict.w.manual.names <- load.from.taiga(data.name='webster-dictionary-with-manual-names-d259', data.version=4, data.file='webster_dict_w_manual_names') %>%
  column_to_rownames("DepMap_ID") %>%
  as.matrix() 


map <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='Model') %>%
  dplyr::distinct(CCLEName, ModelID)


gene.effect <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=13, data.file='gene_effect')
colnames(gene.effect) %<>% word()
rownames(gene.effect) <- tibble(CCLEName = rownames(gene.effect)) %>%
  dplyr::left_join(map) %>%
  .$ModelID

gene.effect <- gene.effect[, apply(gene.effect, 2,  function(x) var(x, na.rm = T)) > 0.01]



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




MOLTEN <- reshape2::melt(Expression) %>%
  dplyr::mutate(type = "EXP") %>% 
  dplyr::filter(is.finite(value)) %>% 
  dplyr::bind_rows(reshape2::melt(CRISPRGeneEffect) %>%
                     dplyr::mutate(type = "XPR")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) %>%
  dplyr::bind_rows(reshape2::melt(MUT) %>%
                     dplyr::mutate(type = "MUT") %>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2))) %>%
  dplyr::bind_rows(reshape2::melt(compound_matrix) %>%
                     dplyr::mutate(type = "REP")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) %>%
  dplyr::bind_rows(reshape2::melt(OmicsSignatures) %>%
                     dplyr::mutate(type = "OMICS.SIG") %>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2))) %>%
  dplyr::bind_rows(reshape2::melt(DriverGeneFunctionalAlterations) %>%
                     dplyr::mutate(type = "ALT")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) %>%
  dplyr::bind_rows(reshape2::melt(Model) %>%
                     dplyr::mutate(type = "LIN")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) %>%
  dplyr::bind_rows(reshape2::melt(OmicsCNGene) %>%
                     dplyr::mutate(type = "CN")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) %>%
  dplyr::bind_rows(reshape2::melt(MethylationExpressionImpactBroad) %>%
                     dplyr::mutate(type = "MET.BROAD")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) %>%
  dplyr::bind_rows(reshape2::melt(MethylationExpressionImpactSanger) %>%
                     dplyr::mutate(type = "MET.GDSC") %>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2))) %>%
  dplyr::bind_rows(reshape2::melt(webster.dict.w.manual.names) %>%
                     dplyr::mutate(type = "WEBSTER")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) )  %>%
  dplyr::bind_rows(reshape2::melt(CCLE.RPPA.20180123) %>%
                     dplyr::mutate(type = "RPPA")%>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2))) %>%
  dplyr::bind_rows(reshape2::melt(gene.effect) %>%
                     dplyr::mutate(type = "RNAi") %>% 
                     dplyr::filter(is.finite(value),
                                   !is.na(Var1),
                                   !is.na(Var2)) ) 


MATRIX <- MOLTEN %>%
  reshape2::acast(Var1 ~ type + Var2, value.var = "value")



MATRIX %>% 
  write.csv("data/BIOMARKER_FEATURES_INTERNAL.csv")

