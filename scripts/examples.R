
library(taigr)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggthemes)
library(scales)

source("scripts/base_functions.R")

CRISPRGeneEffect <- load.from.taiga(data.name='internal-23q2-1e49', data.version=95, data.file='CRISPRGeneEffect')
colnames(CRISPRGeneEffect) %<>% word()
#CRISPRGeneEffect <- CRISPRGeneEffect[, apply(CRISPRGeneEffect, 2,  function(x) var(x, na.rm = T)) > 0.01]

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


