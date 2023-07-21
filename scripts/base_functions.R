
robust_linear_model <- function(X, y, W = NULL) {
  # W shouldn't have any NAs!!!
  Y = matrix(y, dim(X)[1], dim(X)[2])
  
  # Centering X 
  if(!is.null(W)){
    #require(splines)
    #if(is.null(d)){
    #  d <- pmin(round((dim(W)[1]-25)/25/dim(W)[2]), 5)
    #}  
    #if(d > 1){
    #  W <- matrix(apply(W, 2, function(x) splines::ns(x, df = d)),
    #              nrow = dim(W)[1])
    #} 
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
  
  
  # Homoskedastic p-values: Only for comparison! 
  p.val.homoskedastic <- 2*pt(sqrt( (n-d) /(1/rho^2-1)),
                              df = n-d, lower.tail = FALSE)
  # Experimental robustness score
  ns = apply(S,2, function(s) sum(cumsum(sort(s/sum(s, na.rm = T), decreasing = T)) < 1, na.rm = T)) + 1
  
  return(list(x = names(beta), rho = rho, beta = beta,
              p.val.rob = p.val, 
              p.val = p.val.homoskedastic, 
              n = n,
              ns = ns))
}



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

univariate_biomarker_table <- function(Y, W = NULL,
                                       path = NULL,
                                       features = "internal",
                                       lab.n = 5){
  require(taigr);
  require(tidyverse)
  require(magrittr)
  library(ggrepel)
  library(ggthemes)
  
  if(!is.matrix(Y)){
    Y = as.matrix(Y)
  }
  
  if(features == "internal"){
    BIOMARKER.FEATURES <- load.from.taiga(data.name='mustafa-biomarker-features--d353', data.version=2, data.file='BIOMARKER_FEATURES_INTERNAL')
  }else{
    BIOMARKER.FEATURES <- load.from.taiga(data.name='mustafa-biomarker-features--d353', data.version=2, data.file='BIOMARKER_FEATURES_PUBLIC')
  }
  cl = intersect(rownames(BIOMARKER.FEATURES), rownames(Y))
  
  BIOMARKER.FEATURES <- BIOMARKER.FEATURES[cl, ]
  BIOMARKER.FEATURES <- BIOMARKER.FEATURES[, apply(BIOMARKER.FEATURES, 2, function(x) sum(is.finite(x))) > 24, drop = FALSE]
  W = W[cl,, drop = F]
  Y = Y[cl, , drop = F]
  
  f <- function(yy){
    res <- robust_linear_model(X = BIOMARKER.FEATURES,  # !!!
                               y = yy,
                               W = W) %>% 
      as.tibble()%>%
      dplyr::mutate(feature.set = word(x, sep = fixed("_")),
                    feature = word(x,2,-1, sep = fixed("_"))) %>% 
      dplyr::group_by(feature.set) %>%
      dplyr::mutate(q.val = p.adjust(p.val, method = "BH"),
                    q.val.rob = p.adjust(p.val.rob, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is.finite(p.val),
                    ns > 2, q.val < .5)
    return(res)
  }
  
  RES = list(); ix= 1
  for(nn in colnames(Y)){
    print(paste0(ix ,"/", dim(Y)[2], " - ", nn)) 
    RES[[nn]] = f(Y[,nn])
    RES[[nn]] %<>% dplyr::mutate(y = nn) 
    ix = ix + 1
  }
  
  RES %<>% dplyr::bind_rows()
  
  
  if(!is.null(path)){
    RES %>%
      write_csv(paste0(path, ".csv"))
    
    ys <- RES$y %>% unique()
    
    pdf(paste0(path, ".pdf"), width = 14, height = 7)
    
    for(yy in ys){
      p = RES%>%
        dplyr::group_by(feature.set,y) %>%
        dplyr::arrange(rho) %>% dplyr::mutate(ix = 1:n()) %>% 
        dplyr::arrange(-rho) %>% dplyr::mutate(ix = pmin(1:n(), ix)) %>% 
        dplyr::arrange(q.val.rob) %>% dplyr::mutate(ix = pmin(1:n(), ix)) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(type2 = ifelse(feature.set %in% c("LIN", "ALT", "MUT", "RPPA", "OMICS.SIG"), 
                                     ifelse(features == "internal", " Lineage/Driver Alt./Mutations/RPPA/Omics Sig.", " Lineage/Mutations/RPPA"), 
                                     ifelse(feature.set %in% c("EXP", "MET.BROAD", "MET.GDSC", "CN"), 
                                            ifelse(features == "internal", " Expression/Copy Number/Methylation", " Expression/Copy Number"), 
                                            ifelse(feature.set %in% c("XPR", "RNAi", "REP", "WEBSTER"), " CRISPR/RNAi/Repurposing/Webster", "Other")))) %>% 
        dplyr::filter(y == yy,
                      ns > 4) %>% 
        ggplot(aes(x = rho, 
                   y = -log10(q.val),
                   color = feature.set)) +
        geom_point(shape = 1,
                   aes(size = cut(ns, ceiling(4*2^((0:5)*log2(length(cl)/4) / 5)),
                                  include.lowest = T)), 
                   show.legend = T) +
        geom_text_repel(aes(label = ifelse(ix <= lab.n, feature, NA)),
                        color = "black",
                        show.legend = F, size = 3, max.overlaps = 25) + 
        facet_wrap(type2 ~ .,
                   scales = "free_y",
                   ncol = 4) +
        coord_cartesian(ylim = c(0, NA),
                        xlim = c(-1,1)) + 
        theme_clean(#base_family = "GillSans",
          base_size = 14) + 
        theme(legend.position = "bottom") + 
        guides(color=guide_legend(nrow=2,byrow=TRUE, 
                                  label.theme = element_text(size = 10),
                                  title.theme =  element_text(size = 10)),
               size=guide_legend(title.position = "top", 
                                 label.theme = element_text(size = 10), 
                                 title.theme =  element_text(size = 12))) + 
        scale_size_manual(values = c(1, 2, 3, 3.5, 4)) + 
        scale_color_viridis_d(option = "H") + 
        labs(x = "Correlation", y = "-log10(q)", color = "", size = "Stability score (ns)",
             title = yy)
      
      
      g <- ggplotGrob(p)
      ax <- g[["grobs"]][g$layout$name == "axis-l-1-1"][[1]]
      pp = p + 
        geom_vline(aes(xintercept=0)) +
        annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
                          xmax=0, xmin=0) +
        theme(axis.text.y = element_blank(), 
              axis.ticks.y=element_blank(),
              axis.line.y =  element_blank()) +
        coord_cartesian(ylim = c(0.1,NA))
      
      print(pp)
      
    }
    dev.off()
    
  }
 
  return(RES)
}

# write a random_forest_biomarker_table function! 

