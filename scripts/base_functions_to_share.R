
robust_linear_model <- function(X, y, W = NULL) {
  # W shouldn't have any NAs!!!
  Y = matrix(y, dim(X)[1], dim(X)[2])
  
  # Centering X 
  if(!is.null(W)){
    W = cbind(1, W) 
    H = diag(nrow(W)) - crossprod(t(W), tcrossprod(solve(crossprod(W)), W))
    
    X.mask = is.na(X); y.mask = is.na(y); 
    mask <- X.mask | y.mask
    X[mask] = 0; Y[mask] = 0
    X <- crossprod(H, X); X[mask] = NA
    Y <- crossprod(H, Y); Y[mask] = NA
    d = dim(W)[2] + 1
  } else{
    d = 2
  }
  
  flag = apply(X, 2, function(x) var(x,na.rm = T)) 
  Y <- Y[, flag > 0.001, drop = FALSE]
  X <- X[, flag > 0.001, drop = FALSE] # this threshold is very arbitrary! 
  
  # this step is necessary due to the NA's in X !! 
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  
  # Defining S and computing mean and variances 
  S = X * Y
  n = colSums(!is.na(S))
  muS = colMeans(S, na.rm = T)
  #sigmaS = apply(S, 2, function(x) sd(x, na.rm = T))
  sigmaS = sqrt(rowMeans((t(S) - muS)^2, na.rm = T))
  X[is.na(S)] = NA; Y[is.na(S)] = NA
  varX = colMeans(X^2, na.rm = T)
  varY = colMeans(Y^2, na.rm = T)
  
  # Estimation and inference
  beta = muS / varX
  rho = muS / sqrt(varX * varY)
  p.val <- 2*pt(abs(sqrt(n-d) * muS/sigmaS), 
                df = n-d, lower.tail = FALSE)
  
  
  # Homoskedastic p-values
  p.val.homoskedastic <- 2*pt(sqrt( (n-d) /(1/rho^2-1)),
                              df = n-d, lower.tail = FALSE)
  # Experimental robustness score
  ns = apply(S,2, function(s) sum(cumsum(sort(s/sum(s, na.rm = T), decreasing = T)) < 1, na.rm = T)) + 1
  
  ns2 = apply(X, 2, function(x){
    x = x[is.finite(x)]
    length(x) - sort(table(x), decreasing = TRUE)[1]
    })
  
  return(list(x = names(beta), rho = rho, beta = beta,
              p.val.rob = p.val, 
              p.val = p.val.homoskedastic, 
              n = n,
              ns = pmin(ns, ns2)))
}

read_dataset <- function(file = 'https://assets.clue.io/testing/depmap_datasets.h5', dataset, rownames_depmap_ids = TRUE) {
  require(rhdf5)
  if(word(file, sep = fixed("://")) %in% c("s3", "http", "https")){
    s3 = TRUE
  } else{
    s3 = FALSE
  }
  X <- h5read(file, name = paste0(dataset, "/mat"), s3 = s3)
  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = s3)
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


univariate_biomarker_table <- function(Y, W = NULL,
                                       file = 'https://assets.clue.io/testing/depmap_datasets.h5',
                                       features = c('CopyNumber', 'Expression', 'ExpressionGeneSetEnrichment',
                                                    'Fusion', 'GDSC.Proteomics', 'GeneEffect', 'LoH', 'Mutation',
                                                    'RNAi', 'RPPA', 'Repurposing', 'Signatures'),
                                       homoskedastic = TRUE, n.X.min = 25, ns.min = 3, q.val.max = .2){
  require(tidyverse)
  require(magrittr)
  require(parallel)
  require(parallelly)
  require(rhdf5)

  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  
  RESULTS <- list()
  for(feat in features){
    print(feat)
    X <- read_dataset(file , feat)
    cl = intersect(rownames(X), rownames(Y))
    
    X <- X[cl, ];   W = W[cl, , drop = F];    Y = Y[cl, , drop = F]
    X <- X[ , apply(X, 2, function(x) sum(is.finite(x))) >= n.X.min, drop = FALSE]
    
    if((dim(X)[1] >= n.X.min) & dim(X)[2] > 0){
      f <- function(ix){
        res <- robust_linear_model(X = X, y = Y[,ix], W = W) %>% 
          as.tibble()%>%
          dplyr::rename(feature = x) %>%
          dplyr::mutate(q.val = p.adjust(p.val, method = "BH"),
                        q.val.rob = p.adjust(p.val.rob, method = "BH"),
                        feature.set = feat, 
                        y = colnames(Y)[ix]) %>%
          dplyr::ungroup() %>%
          dplyr::filter(ns >= ns.min)
        
        if(homoskedastic){
          res <- res %>%
            dplyr::filter(q.val <= q.val.max)
        } else{
          res <- res %>%
            dplyr::filter(q.val.rob <= q.val.max)
        }
        
        if(nrow(res) > 0){
          res <- res %>%
            dplyr::arrange(rho) %>% 
            dplyr::mutate(rank = 1:n()) %>% 
            dplyr::arrange(-rho) %>% 
            dplyr::mutate(rank = pmin(rank, 1:n()))
        }
        return(res)
      }
      
      RES <- mclapply(1:dim(Y)[2], f, mc.cores = parallelly::availableCores() )
      RESULTS[[feat]] <- dplyr::bind_rows(RES)
    }
  }
  
  RESULTS <- dplyr::bind_rows(RESULTS)
  return(RESULTS)
}



calculate_fitted_fold_changes <- function(Lower_Limit, Upper_Limit, Slope, Inflection, md, MD, k = 12){
  g <- function(x) Lower_Limit + (Upper_Limit - Lower_Limit) / (1 + (x / Inflection)^-Slope)
  tibble(dose = 2^seq(log2(md), log2(MD), length.out = k),
         fc = g(dose))
}

# write a random_forest_biomarker_table function! 
random_forest <- function (X, y, W = NULL, k = 10, vc = 0.01, lm = 25, nu = 10, 
                           p0 = 0.01) {
  y.clean <- y[is.finite(y)]
  X.clean <- X[, apply(X, 2, function(x) all(is.finite(x)))]
  cl <- sample(dplyr::intersect(rownames(X.clean), names(y.clean)))
  X.clean <- X.clean[cl, ]
  y.clean <- y.clean[cl]
  colnames(X.clean) %<>% make.names()
  N = floor(length(cl)/k)
  yhat_rf <- rep(NA, length(y.clean))
  names(yhat_rf) <- names(y.clean)
  SS = tibble()
  for (kx in 1:k) {
    test <- cl[((kx - 1) * N + 1):(kx * N)]
    if (kx == k) 
      test <- cl[((kx - 1) * N + 1):length(cl)]
    train <- dplyr::setdiff(cl, test)
    X_train <- X.clean[train, , drop = F]
    X_test <- X.clean[test, , drop = F]
    X_train <- X_train[, apply(X_train, 2, var) > vc, drop = F]
    b <- gausscov::f2st(y.clean[train], X_train, lm = lm, 
                        nu = nu, p0 = p0)
    features <- colnames(X_train)[b[[1]][, 2]]
    if (length(features > 0)) {
      X_train <- X_train[, features, drop = F]
      if (!is.null(W)) {
        X_train <- cbind(X_train, W[train, ])
        X_test <- cbind(X_test, W[test, ])
      }
      rf <- ranger::ranger(y = y.clean[train], x = X_train, 
                           importance = "impurity")
      yhat_rf[test] <- predict(rf, data = as.data.frame(X_test[, 
                                                               colnames(X_train), drop = F]))$predictions
      ss <- tibble::tibble(feature = names(rf$variable.importance), 
                           RF.imp = rf$variable.importance/sum(rf$variable.importance), 
                           fold = kx)
      SS %<>% dplyr::bind_rows(ss)
    }
  }
  if (nrow(SS) > 0) {
    RF.importances <- SS %>% dplyr::distinct(feature, RF.imp, 
                                             fold) %>% reshape2::acast(feature ~ fold, value.var = "RF.imp", 
                                                                       fill = 0)
    RF.table <- tibble::tibble(feature = rownames(RF.importances), 
                               RF.imp.mean = RF.importances %>% apply(1, mean), 
                               RF.imp.sd = RF.importances %>% apply(1, function(x) sd(x, 
                                                                                      na.rm = T)), RF.imp.stability = RF.importances %>% 
                                 apply(1, function(x) sum(x > 0)/k)) %>% dplyr::filter(feature != 
                                                                                         "(Intercept)") %>% dplyr::arrange(desc(RF.imp.mean)) %>% 
      dplyr::mutate(rank = 1:n())
    mse <- mean((yhat_rf - y.clean)^2, na.rm = T)
    mse.se <- sqrt(var((yhat_rf - y.clean)^2, na.rm = T))/sqrt(length(y.clean))
    r2 <- 1 - (mse/var(y.clean, na.rm = T))
    ps <- cor(y.clean, yhat_rf, use = "pairwise.complete.obs")
    RF.table %<>% dplyr::mutate(MSE = mse, MSE.se = mse.se, 
                                R2 = r2, PearsonScore = ps)
    return(list(model_table = RF.table, predictions = yhat_rf))
  }
  else {
    return(list(model_table = tibble(), predictions = tibble()))
  }
}


# EXAMPLES -----------

library(tidyverse)
library(taigr)


## make sure to change the taigaclient path below if you are using newer version of taiga 
# options(taigaclient.path=path.expand("/opt/miniconda3/envs/taigapy/bin/taigaclient"))
PRISMOncologyReferenceCompoundList <- load.from.taiga(data.name='prism-oncology-reference-set-24q2-f31b', data.version=4, data.file='PRISMOncologyReferenceCompoundList')
PRISMOncologyReferenceAUCMatrix <- load.from.taiga(data.name='prism-oncology-reference-set-24q2-f31b', data.version=4, data.file='PRISMOncologyReferenceAUCMatrix')



colnames(PRISMOncologyReferenceAUCMatrix) <- tibble(SampleID = colnames(PRISMOncologyReferenceAUCMatrix)) %>% 
  dplyr::left_join(PRISMOncologyReferenceCompoundList %>% 
                     dplyr::distinct(SampleID, CompoundName)) %>% 
  .$CompoundName


Y <- PRISMOncologyReferenceAUCMatrix[, c("MRTX-1133", "ZIFTOMENIB", "WO2022249060 EXAMPLE 42")]

# Let's calculate the biomarkers for only Mutations, Fusions, and Signatures 
# Note the the full list of available featuresets are : 'CopyNumber', 'Expression', 'ExpressionGeneSetEnrichment', 'Fusion', 'GDSC.Proteomics', 'GeneEffect', 'LoH', 'Mutation', 'RNAi', 'RPPA', 'Repurposing', 'Signatures'

# DO NOT RUN! 
# If we don't profile a "file" argument it uses the features from the AWS: 'https://assets.clue.io/testing/depmap_datasets.h5', which can be slow and heavy. 
# Also it sometimes give random errors, I am stil trying to fix it! :) 
results <- univariate_biomarker_table(Y, features = c( 'Mutation', 'Signatures'))


# If you are going to use this function often, I recommend downloading depmap_datasets.h5 file to you local disc and providing the file path: 
results_local <- univariate_biomarker_table(Y, file = "~/Downloads/depmap_datasets.h5", features = c( 'Mutation', "Fusion", 'Signatures'))


# vanity libraries
library(ggrepel)
library(scales)
library(ggthemes)


results_local %>% 
  ggplot(aes(x = rho, y = -log10(q.val), color = feature.set)) +
  geom_point(show.legend = FALSE) +
  geom_text_repel(aes(label = ifelse(rank < 3, feature, NA)),
                  show.legend = FALSE) +
  facet_grid(feature.set ~ y, scales = "free_y") +
  scale_color_wsj() + theme_bw() +
  labs(x = "Correlation", y = "-log10(q)", color = NULL)





# Experimental example - feel free to ignore very much wip. -----

PRISMOncologyReferenceAnalyteMeta <- load.from.taiga(data.name='prism-oncology-reference-set-24q2-f31b', data.version=4, data.file='PRISMOncologyReferenceAnalyteMeta')
PRISMOncologyReferenceDoseResponseParameters <- load.from.taiga(data.name='prism-oncology-reference-set-24q2-f31b', data.version=4, data.file='PRISMOncologyReferenceDoseResponseParameters')


k = 48

test_data <- PRISMOncologyReferenceDoseResponseParameters %>%  
  dplyr::left_join(PRISMOncologyReferenceAnalyteMeta) %>% 
  dplyr::left_join(PRISMOncologyReferenceCompoundList) %>% 
  dplyr::filter(CompoundName %in% c("MRTX-1133", "ZIFTOMENIB",
                                    "WO2022249060 EXAMPLE 42")) %>% 
  dplyr::group_by(CompoundName, CompoundPlate, depmap_id) %>% 
  dplyr::top_n(1, cellset) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(CompoundName, SampleID, CompoundPlate, depmap_id,
                  Lower_Limit, Upper_Limit, Slope, Inflection, md, MD) 


molten_data <- test_data %>%   
  dplyr::group_by(CompoundName, CompoundPlate, SampleID, depmap_id,
                  Lower_Limit, Upper_Limit, Slope, Inflection, md, MD) %>% 
  dplyr::summarize(calculate_fitted_fold_changes(Lower_Limit, Upper_Limit, Slope, Inflection, md, MD, k = k)) %>%
  dplyr::ungroup()  %>% 
  dplyr::mutate(cn = paste0(CompoundName, "::", CompoundPlate, "::", SampleID, "::", dose)) 


test_biomarkers <- molten_data %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "fc") %>%
  univariate_biomarker_table(file = "~/Downloads/depmap_datasets.h5",
                             features =c( 'Fusion', 'Mutation', 'Signatures'), q.val.max = 1, ns.min = 0) %>% 
  dplyr::left_join(molten_data %>% 
                     dplyr::distinct(CompoundName, SampleID, CompoundPlate, dose, cn) ,
                   by = c("y" = "cn")) 

test_biomarkers %>% 
  dplyr::group_by(feature, feature.set, CompoundName, CompoundPlate) %>% 
  dplyr::summarize(e.agg = mean(1/sqrt(p.val) - 1, na.rm = T),
                   nsm = max(ns, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(p.agg = pmin(1, 1/e.agg)) %>%  
  dplyr::group_by(feature.set, CompoundName, CompoundPlate) %>% 
  dplyr::mutate(q = p.adjust(p.agg, method = "BH")) %>% 
  dplyr::filter(nsm > 2) %>% 
  dplyr::top_n(20, -log10(q)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(q < 0.2) %>% 
  dplyr::left_join(test_biomarkers) %>%  
  dplyr::group_by(feature, feature.set, CompoundName) %>%  
  dplyr::mutate(mc = max(abs(rho), na.rm = T)) %>% 
  dplyr::ungroup()  %>% 
  ggplot(aes(x = log10(dose),
             y = rho,
             group = paste0(feature, q),
             color = -log10(q.val),
             color = ns
             )) +
  geom_hline(aes(yintercept = 0), color ="black") + 
  geom_line(lwd = 1) +
  geom_text_repel(aes(label = ifelse(abs(rho) == mc, feature, NA)),
                  color = "red", size = 3) + 
  theme_base() +
  scale_color_viridis_b(option = "B") + 
  facet_grid(feature.set ~ CompoundName,
             scales = "free_x")
  
