
library(taigr)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggthemes)
library(scales)

source("scripts/base_functions.R")

CRISPRGeneEffect <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='CRISPRGeneEffect')
colnames(CRISPRGeneEffect) %<>% word()

results <- robust_linear_model(CRISPRGeneEffect, CRISPRGeneEffect[, "CTNNB1"])


p1 =results %>%  
  as.tibble() %>% 
  dplyr::filter(x != "CTNNB1") %>% 
  dplyr::mutate(q.val = p.adjust(p.val, method = "BH"),
                q.val.rob = p.adjust(p.val.rob, method = "BH")) %>% 
  dplyr::arrange(rho) %>% dplyr::mutate(ix = 1:n()) %>% 
  dplyr::arrange(-rho) %>% dplyr::mutate(ix = pmin(ix, 1:n())) %>% 
  ggplot(aes(x = rho, y = -log10(q.val))) +
  geom_point(aes(color = log2(ns)), shape = 1, size = 2) +
  geom_text_repel(aes(label = ifelse(ix < 6, x, NA))) +
  scale_color_viridis_b() + theme_bw()

p2 =results %>%  
  as.tibble() %>% 
  dplyr::filter(x != "CTNNB1") %>% 
  dplyr::mutate(q.val = p.adjust(p.val, method = "BH"),
                q.val.rob = p.adjust(p.val.rob, method = "BH")) %>% 
  dplyr::arrange(rho) %>% dplyr::mutate(ix = 1:n()) %>% 
  dplyr::arrange(-rho) %>% dplyr::mutate(ix = pmin(ix, 1:n())) %>% 
  ggplot(aes(x = rho, y = -log10(q.val.rob))) +
  geom_point(aes(color = log2(ns)), shape = 1, size = 2) +
  geom_text_repel(aes(label = ifelse(ix < 6, x, NA))) +
  scale_color_viridis_b() + theme_bw()

genes <- results %>%  
  as.tibble() %>% 
  dplyr::filter(x != "CTNNB1") %>% 
  dplyr::mutate(q.val = p.adjust(p.val, method = "BH"),
                q.val.rob = p.adjust(p.val.rob, method = "BH")) %>% 
  dplyr::arrange(rho) %>% dplyr::mutate(ix = 1:n()) %>% 
  dplyr::arrange(-rho) %>% dplyr::mutate(ix = pmin(ix, 1:n())) %>%
  dplyr::filter(ix < 4) %>% 
  .$x  %>% unique


p3 = CRISPRGeneEffect[, genes] %>% 
  reshape2::melt() %>% 
  dplyr::left_join(tibble(Var1 = rownames(CRISPRGeneEffect),
                          CTNNB1 = CRISPRGeneEffect[,"CTNNB1"])) %>% 
  dplyr::filter(is.finite(value), is.finite(CTNNB1)) %>% 
  dplyr::group_by(Var2) %>% dplyr::mutate(n = n(), r = cor(value, CTNNB1)) %>% 
  ggplot() +
  geom_point(aes(x = CTNNB1, y = value,
                 color = Var2), show.legend = F) +
  geom_abline() + 
  facet_wrap(reorder(paste0(Var2, "\nr = ", round(r,2)), -r) ~ .) +
  scale_color_wsj() + 
  theme_bw() + labs(x = "CTNNB1 Gene Effect", y = "Gene Effect")


cowplot::plot_grid(cowplot::plot_grid(p1,p2, nrow = 1),
                   p3, ncol = 1, rel_heights = c(2,3))



res_internal <- univariate_biomarker_table(Y = CRISPRGeneEffect[, "CTNNB1", drop = FALSE], 
                           path = "results/CTNNB1_internal")

res_public <- univariate_biomarker_table(Y = CRISPRGeneEffect[, "CTNNB1", drop = FALSE], 
                                           features = "public", 
                                           path = "results/CTNNB1_public")


# ONC REF # don't run!  ----  

AUC.matrix <- load.from.taiga(data.name='prism-oncref001-portal-files--3044', data.version=4, data.file='AUC_matrix')
IC50.matrix <- load.from.taiga(data.name='prism-oncref001-portal-files--3044', data.version=4, data.file='IC50_matrix')
viability.matrix <- load.from.taiga(data.name='prism-oncref001-portal-files--3044', data.version=4, data.file='viability_matrix')

Condition.annotations <- load.from.taiga(data.name='prism-oncref001-portal-files--3044', data.version=4, data.file='Condition_annotations')
Model.annotations <- load.from.taiga(data.name='prism-oncref001-portal-files--3044', data.version=4, data.file='Model_annotations')
compound.list <- load.from.taiga(data.name='prism-oncref001-portal-files--3044', data.version=4, data.file='compound_list')

V = viability.matrix %>% 
  reshape2::melt() %>% 
  dplyr::left_join(Condition.annotations, by = c("Var1" = "index")) %>%
  dplyr::left_join(Model.annotations, by = c("Var2" = "index")) %>% 
  dplyr::group_by(compound_name, cell_line_name, dose) %>% 
  dplyr::summarise(LFC = median(value, na.rm = T)) %>%
  dplyr::mutate(dose = as.character(dose)) %>% 
  dplyr::bind_rows(AUC.matrix %>% 
                     reshape2::melt() %>%
                     dplyr::rename(LFC = value, compound_name = Var1, cell_line_name = Var2) %>% 
                     dplyr::mutate(dose = "log2(AUC)", LFC = log2(LFC)) ) %>% 
  dplyr::bind_rows(IC50.matrix %>% 
                     reshape2::melt() %>%
                     dplyr::rename(LFC = value, compound_name = Var1, cell_line_name = Var2) %>% 
                     dplyr::mutate(dose = "log2(IC50)", LFC = log2(LFC)) ) %>%
  dplyr::left_join(compound.list %>% 
                     dplyr::distinct(IDs, Drug.Name, Target),
                   by = c("compound_name" = "IDs")) %>% 
  dplyr::mutate(cn = paste0(compound_name, "::", dose,
                            "::", Drug.Name, "::", Target)) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(LFC)) %>% 
  reshape2::acast(cell_line_name ~ cn, value.var = "LFC")


res_onc_public <- univariate_biomarker_table(Y = V, 
                                         features = "public", 
                                         path = "results/onc_ref")

res_onc_public <- data.table::fread("results/onc_ref.csv") %>% 
  dplyr::filter(ns > 9, n > 24) %>% 
  dplyr::group_by(y, feature.set) %>% 
  dplyr::arrange(rho) %>% dplyr::mutate(ix = 1:n()) %>% 
  dplyr::arrange(-rho) %>% dplyr::mutate(iix = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter((ix < 26) | (iix < 26)) %>% 
  dplyr::mutate(IDs = word(y, sep = fixed("::")),
                dose = word(y,2, sep = fixed("::")),
                Drug.Name = word(y, 3, sep = fixed("::")),
                target = word(y, 4, sep = fixed("::"))) %>%
  dplyr::group_by(IDs, x) %>%
  dplyr::mutate(m = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(m > 1)


res_onc_public %>% 
  dplyr::group_by(feature.set, IDs) %>%
  dplyr::top_n(1, m) %>% 
  dplyr::top_n(1, abs(rho)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(feature.set %in% c("REP", "MUT", "EXP", "XPR", "RPPA", "RNAi")) %>% 
  dplyr::mutate(feature = paste0(feature," - m:", m , " - r:", round(rho,2), "  @ ", dose)) %>% 
  dplyr::select(-rho, -dose, -beta, -n, -ns, -q.val, -q.val.rob, -ix, -iix, -x, -y, -m, -p.val, -p.val.rob) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(Drug.Name != "NA") %>% 
  dplyr::group_by(IDs, feature.set) %>%
  tidyr::pivot_wider(names_from = "feature.set", values_from = "feature") %>%
  View
