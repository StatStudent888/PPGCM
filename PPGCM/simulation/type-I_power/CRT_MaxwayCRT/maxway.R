# ============================================================
# Model-X CRT and Maxway CRT implementation
#
# Purpose:
#   This script implements the model-X CRT and Maxway CRT procedures
#   used as competing methods in the simulation studies. The main
#   function Maxway_CRT() fits working models for X | Z and Y | Z,
#   calibrates the fitted X model using low-dimensional information
#   extracted from the Y model, and returns p-values for the model-X
#   and Maxway CRT variants.
#
# Notation:
#   A_label : labeled exposure X
#   y_label : labeled outcome Y
#   Z_label : labeled conditioning variables Z
#   A_add   : additional/unlabeled observations of X
#   Z_add   : corresponding additional/unlabeled observations of Z
#
# Supported working models:
#   - Gaussian lasso
#   - Binomial lasso
#   - Gaussian random forest
#
# Main output:
#   Maxway_CRT() returns p-values for the model-X CRT and the
#   corresponding Maxway CRT procedures. The exact returned components
#   depend on the selected models for X and Y.
#
# Main required packages:
#   glmnet, randomForest, MASS, CompQuadForm
# ============================================================

# ------------------------------------------------------------
# Main function for the model-X CRT and Maxway CRT.
#
# The procedure first fits X | Z and Y | Z, extracts a low-dimensional
# summary of the outcome regression, and then recalibrates the fitted
# X model using this summary. Both the original model-X CRT and the
# calibrated Maxway CRT p-values are returned.
#
# A_add and Z_add provide additional observations for estimating X | Z.
# Z_label_add and y_label_add can optionally provide an independent
# sample for out-of-sample estimation of the Y regression.
# ------------------------------------------------------------
Maxway_CRT <- function(A_label, Z_label, y_label, A_add = NULL, Z_add = NULL, 
                       Z_label_add = NULL, y_label_add = NULL,
                       model_x = 'Gaussian_lasso', model_y = 'Gaussian_lasso', 
                       RF.num.trees = c(100, 100, 50, 20), 
                       lambda.seq = NULL, k = NULL, seed = 1, M = 500) {
  
  # Check that the requested working models are supported.
  if (!model_y %in% c('Gaussian_lasso', 'Binomial_lasso', 'Gaussian_RF')){
    print('Error: please correctly specify the model for Y.')
    return()
  }
  
  if (!model_x %in% c('Gaussian_lasso', 'Binomial_lasso', 'Gaussian_RF')){
    print('Error: please correctly specify the model for X.')
    return()
  }
  # If no additional samples are supplied, use the labeled sample itself.
  # Likewise, use the labeled sample for the auxiliary Y fit by default.
  if (is.null(A_add)){
    A_add <- A_label
    Z_add <- Z_label
  }
  
  if (is.null(Z_label_add)){
    Z_label_add <- Z_label
    y_label_add <- y_label
  }
  
  p <- ncol(Z_label)
  n <- length(A_label)
  N <- length(A_add)
  
  # Choose the dimension of the low-dimensional adjustment set
  # automatically when k is not specified by the user.
  if (is.null(k)){
    if (model_x == 'Gaussian_lasso' | model_x == 'Binomial_lasso'){
      k <- as.integer(min(2 * log(p), p - 1, N^(1/3)))
    }
    if (model_x == 'Gaussian_RF'){
      k <- as.integer(2 * log(p))
    }
    if (model_y == 'Binomial_lasso'){
      k <- as.integer(log(p))
    }
  }
  
  
  # ============================================================
  # Fit the working model for X | Z
  # ============================================================
  # The original Maxway code denotes the exposure X by A.
  
  # Continuous X: lasso regression.
  
  if (model_x == 'Gaussian_lasso'){
    set.seed(seed)
    fit_x <- fit_gaussian(Z_add, A_add, Z_label, A_label, 
                          model = 'linear', lambda.seq = lambda.seq)
    res_x_lab <- fit_x$res_test
    res_x_add <- fit_x$res_train
  }
  
  # Binary X: logistic lasso regression.
  
  if (model_x == 'Binomial_lasso'){
    set.seed(seed)
    fit_x <- fit_binary(Z_add, A_add, Z_label, A_label, 
                        model = 'linear', lambda.seq = lambda.seq)
    res_x_lab <- fit_x$res_test
    res_x_add <- fit_x$res_train
    pred_x_lab <- fit_x$pred_test
    pred_x_add <- fit_x$pred_train
  }
  
  # Continuous X: random-forest regression.
  
  if (model_x == 'Gaussian_RF'){
    set.seed(seed)
    fit_x <- fit_gaussian(Z_add, A_add, Z_label, A_label, 
                          model = 'RF', RF.num.trees = RF.num.trees[1], CV = T)
    res_x_lab <- fit_x$res_test
    res_x_add <- fit_x$res_train
  }
  
  
  
  # ============================================================
  # Fit the working model for Y | Z
  # ============================================================
  # The fitted outcome model is used both to obtain Y residuals and to
  # construct a low-dimensional summary of Z for the Maxway calibration.
  
  if (model_y == 'Gaussian_lasso'){
    
    # Fit Y | Z on the labeled sample and identify the k covariates
    # with the largest fitted coefficients.
    
    set.seed(seed)
    fit_y <- fit_gaussian(Z_label, y_label, Z_add, NULL, 
                          model = 'linear', lambda.seq = lambda.seq)
    res_y_lab <- fit_y$res_train
    beta_y <- fit_y$coef
    gZ_label <- fit_y$pred_train
    gZ_add <- fit_y$pred_test
    
    # Select the k most important Z variables according to the fitted Y model.
    beta_sort <- sort(abs(as.vector(beta_y)), decreasing = T, index.return = T)
    index_gZ <- beta_sort$ix[1:k]
    
    # Augment the selected Z variables with the component of the fitted
    # outcome regression orthogonal to them.
    gZ_label <- orthogonalize(gZ_label, Z_label[,index_gZ])
    gZ_add <- orthogonalize(gZ_add, Z_add[,index_gZ])
    
    # Repeat the outcome-model construction using the optional independent
    # sample, producing an out-of-sample version of the Maxway adjustment.
    
    set.seed(seed)
    fit_y_out <- fit_gaussian(Z_label_add, y_label_add, rbind(Z_add, Z_label), NULL, 
                              model = 'linear', lambda.seq = lambda.seq)
    beta_y_out <- fit_y_out$coef
    gZ_add_out <- fit_y_out$pred_test[1:nrow(Z_add)]
    gZ_label_out <- fit_y_out$pred_test[(nrow(Z_add) + 1):(nrow(Z_add) + nrow(Z_label))]
    res_y_lab_out <- y_label - gZ_label_out
    
    beta_sort_out <- sort(abs(as.vector(beta_y_out)), decreasing = T, index.return = T)
    index_gZ_out <- beta_sort_out$ix[1:k]
    
    gZ_label_out <- orthogonalize(gZ_label_out, Z_label[,index_gZ_out])
    gZ_add_out <- orthogonalize(gZ_add_out, Z_add[,index_gZ_out])
    
  }
  
  # Binary Y: logistic-lasso outcome model.
  if (model_y == 'Binomial_lasso'){
    
    # In sample or CV fit:
    
    set.seed(seed)
    fit_y <- fit_binary(Z_label, y_label, Z_add, NULL, 
                        model = 'linear', lambda.seq = lambda.seq)
    beta_y <- fit_y$coef
    offset_label <- logit(fit_y$pred_train)
    offset_add <- logit(fit_y$pred_test)
    
    beta_sort <- sort(abs(as.vector(beta_y)), decreasing = T, index.return = T)
    index_gZ <- beta_sort$ix[1:k]
    
    gZ_label <- cbind(logit(pred_x_lab), Z_label[,index_gZ]) 
    gZ_add <- cbind(logit(pred_x_add), Z_add[,index_gZ]) 
    
    # Out of sample fit:
    
    set.seed(seed)
    fit_y_out <- fit_binary(Z_label_add, y_label_add, rbind(Z_add, Z_label), NULL, 
                            model = 'linear', lambda.seq = lambda.seq)
    beta_y_out <- fit_y_out$coef
    offset_add_out <- logit(fit_y_out$pred_test[1:nrow(Z_add)])
    offset_label_out <- logit(fit_y_out$pred_test[(nrow(Z_add) + 1):(nrow(Z_add) + nrow(Z_label))])
    
    beta_sort_out <- sort(abs(as.vector(beta_y_out)), decreasing = T, index.return = T)
    index_gZ_out <- beta_sort_out$ix[1:k]
    
    gZ_label_out <- orthogonalize(offset_label_out, Z_label[,index_gZ_out])
    gZ_add_out <- orthogonalize(offset_add_out, Z_add[,index_gZ_out])
    
  }
  
  # Continuous Y: random-forest outcome model.
  if (model_y == 'Gaussian_RF'){
    
    # In sample or CV fit:
    
    set.seed(seed)
    fit_y <- fit_gaussian(Z_label, y_label, Z_add, NULL, 
                          model = 'RF', RF.num.trees = RF.num.trees[2], CV = T)
    res_y_lab <- fit_y$res_train
    gZ_label <- fit_y$pred_train
    gZ_add <- fit_y$pred_test
    
    RF_y_importance <- fit_y$model$importance
    imp_fit <- as.vector(RF_y_importance)
    imp_sort <- sort(imp_fit, decreasing = T, index.return = T)
    index_gZ <- imp_sort$ix[1:k]
    
    gZ_label <- cbind(gZ_label, Z_label[,index_gZ])
    gZ_add <- cbind(gZ_add, Z_add[,index_gZ])
    
    # Out of sample fit:
    
    set.seed(seed)
    fit_y_out <- fit_gaussian(Z_label_add, y_label_add, rbind(Z_add, Z_label), NULL, 
                              model = 'RF', RF.num.trees = RF.num.trees[2], CV = T)
    gZ_add_out <- fit_y_out$pred_test[1:nrow(Z_add)]
    gZ_label_out <- fit_y_out$pred_test[(nrow(Z_add) + 1):(nrow(Z_add) + nrow(Z_label))]
    res_y_lab_out <- y_label - gZ_label_out
    
    RF_y_importance_out <- fit_y_out$model$importance
    imp_fit_out <- as.vector(RF_y_importance_out)
    imp_sort_out <- sort(imp_fit_out, decreasing = T, index.return = T)
    index_gZ_out <- imp_sort_out$ix[1:k]
    
    gZ_label_out <- cbind(gZ_label_out, Z_label[,index_gZ_out])
    gZ_add_out <- cbind(gZ_add_out, Z_add[,index_gZ_out])
    
  }
  
  # ============================================================
  # Calibrate the fitted X | Z model
  # ============================================================
  # Regress the X-model residuals (or probabilities for binary X) on the
  # low-dimensional information extracted from Y | Z. This calibration
  # is the key step distinguishing Maxway CRT from the original model-X CRT.
  
  # Calibration for continuous X fitted by lasso.
  if (model_x == 'Gaussian_lasso'){
    cal_fit <- cal_fit_gaussian(gZ_add, res_x_add, gZ_label, res_x_lab)
    res_x_lab_cal <- cal_fit$res_test
    
    cal_fit_out <- cal_fit_gaussian(gZ_add_out, res_x_add, gZ_label_out, res_x_lab)
    res_x_lab_cal_out <- cal_fit_out$res_test
  }
  
  # Calibration for binary X.
  if (model_x == 'Binomial_lasso'){
    cal_fit <- cal_fit_binary(gZ_add, A_add, pred_x_add, gZ_label, A_label, pred_x_lab)
    pred_x_lab_cal <- cal_fit$pred_test
    
    cal_fit_out <- cal_fit_binary(gZ_add_out, A_add, pred_x_add,
                                  gZ_label_out, A_label, pred_x_lab)
    pred_x_lab_cal_out <- cal_fit_out$pred_test
  }
  
  # Calibration for continuous X fitted by random forest.
  if (model_x == 'Gaussian_RF'){
    cal_fit <- cal_fit_gaussian(gZ_add, res_x_add, gZ_label, res_x_lab, 
                                RF.num.trees = RF.num.trees[3], model = 'RF')
    res_x_lab_cal <- cal_fit$res_test
    
    cal_fit_out <- cal_fit_gaussian(gZ_add_out, res_x_add, gZ_label_out, res_x_lab, 
                                    RF.num.trees = RF.num.trees[3], model = 'RF')
    res_x_lab_cal_out <- cal_fit_out$res_test
  }
  
  
  # ============================================================
  # Compute CRT p-values
  # ============================================================
  # Depending on the types of X and Y, use the corresponding dCRT
  # statistic and return both the original model-X and calibrated
  # Maxway p-values.
  
  # Continuous X and continuous Y with lasso working models.
  if (model_y == 'Gaussian_lasso' & model_x == 'Gaussian_lasso'){
    
    # Original model-X dCRT.
    d0CRT_pvl <- dCRT_res(res_x_lab, res_y_lab, Z_label[,index_gZ], 
                          k = k, d.interaction = F)
    
    dICRT_pvl <- dCRT_res(res_x_lab, res_y_lab, Z_label[,index_gZ], 
                          k = k, d.interaction = T)
    
    # In-sample calibrated Maxway dCRT.
    cald0CRT_pvl <- dCRT_res(res_x_lab_cal, res_y_lab, Z_label[,index_gZ], 
                             k = k, d.interaction = F)
    
    caldICRT_pvl <- dCRT_res(res_x_lab_cal, res_y_lab, Z_label[,index_gZ], 
                             k = k, d.interaction = T)
    
    # Out-of-sample calibrated Maxway dCRT.
    cald0CRT_pvl_out <- dCRT_res(res_x_lab_cal_out, res_y_lab_out, Z_label[,index_gZ_out], 
                                 k = k, d.interaction = F)
    
    caldICRT_pvl_out <- dCRT_res(res_x_lab_cal_out, res_y_lab_out, Z_label[,index_gZ_out], 
                                 k = k, d.interaction = T)
    
    return(list(MX_d0CRT_pvl = d0CRT_pvl, MX_dICRT_pvl = dICRT_pvl, 
                Maxway_d0CRT_pvl = cald0CRT_pvl, Maxway_dICRT_pvl = caldICRT_pvl,
                Maxway_out_d0CRT_pvl = cald0CRT_pvl_out,
                Maxway_out_dICRT_pvl = caldICRT_pvl_out))
  }
  
  # Binary X and continuous Y: use Monte Carlo conditional resampling
  # from the fitted Bernoulli model for X | Z.
  if (model_y == 'Gaussian_lasso' & model_x == 'Binomial_lasso'){
    
    
    cald0CRT_pvl <- dCRT_res_binary(A_label, pred_x_lab_cal, res_y_lab, 
                                    Z_label[,index_gZ], k = 0, d.interaction = F, M = M)
    
    caldICRT_pvl <- dCRT_res_binary(A_label, pred_x_lab_cal, res_y_lab,
                                    Z_label[,index_gZ], k = k, d.interaction = T, M = M)
    
    
    cald0CRT_pvl_out <- dCRT_res_binary(A_label, pred_x_lab_cal_out, res_y_lab_out,
                                        Z_label[,index_gZ_out], k = 0, d.interaction = F,
                                        M = M)
    
    caldICRT_pvl_out <- dCRT_res_binary(A_label, pred_x_lab_cal_out, res_y_lab_out,
                                        Z_label[,index_gZ_out], k = k, d.interaction = T,
                                        M = M)
    
    d0CRT_pvl <- dCRT_res_binary(A_label, pred_x_lab, res_y_lab, Z_label[,index_gZ], 
                                 k = 0, d.interaction = F, M = M)
    
    dICRT_pvl <- dCRT_res_binary(A_label, pred_x_lab, res_y_lab, Z_label[,index_gZ],
                                 k = k, d.interaction = T, M = M)
    
    return(list(Maxway_d0CRT_pvl = cald0CRT_pvl, Maxway_dICRT_pvl = caldICRT_pvl,
                Maxway_out_d0CRT_pvl = cald0CRT_pvl_out, Maxway_out_dICRT_pvl = caldICRT_pvl_out,
                MX_d0CRT_pvl = d0CRT_pvl, MX_dICRT_pvl = dICRT_pvl))
  }
  
  
  
  
  
  if (model_y == 'Gaussian_RF' & model_x == 'Gaussian_RF'){
    
    dCRT_pvl <- dCRT_res_nonp(res_x_lab, res_y_lab, gZ_label,
                              method = 'RF', RF.num.trees = RF.num.trees[4], M = M)

    caldCRT_pvl <- dCRT_res_nonp(res_x_lab_cal, res_y_lab, gZ_label,
                                 method = 'RF', RF.num.trees = RF.num.trees[4], M = M)
    
    caldCRT_pvl_out <- dCRT_res_nonp(res_x_lab_cal_out, res_y_lab_out, gZ_label_out,
                                     method = 'RF', RF.num.trees = RF.num.trees[4], M = M)
    
    return(list(MX_CRT_pvl = dCRT_pvl, Maxway_CRT_pvl = caldCRT_pvl,
                Maxway_out_CRT_pvl = caldCRT_pvl_out))
    
  }
  
  
  # ============================================================
  # Additional case not used in the present simulation study
  # ============================================================
  # Binary X and binary Y are retained for completeness but are not
  # required for the simulations reported in the paper.
  
  if (model_y == 'Binomial_lasso' & model_x == 'Binomial_lasso'){
    
    
    cald0CRT_pvl <- dCRT_binary_Y(A_label, pred_x_lab_cal, y_label, offset_label,
                                  Z_label[,index_gZ], k = 0, d.interaction = F, M = M)
    
    caldICRT_pvl <- dCRT_binary_Y(A_label, pred_x_lab_cal, y_label, offset_label,
                                  Z_label[,index_gZ], k = k, d.interaction = T, M = M)

    d0CRT_pvl <- dCRT_binary_Y(A_label, pred_x_lab, y_label, offset_label,
                               Z_label[,index_gZ], k = 0, d.interaction = F, M = M)
    
    dICRT_pvl <- dCRT_binary_Y(A_label, pred_x_lab, y_label, offset_label,
                               Z_label[,index_gZ], k = k, d.interaction = T, M = M)
    
    return(list(Maxway_d0CRT_pvl = cald0CRT_pvl, Maxway_dICRT_pvl = caldICRT_pvl,
                MX_CRT_pvl = d0CRT_pvl, MX_CRT_pvl = dICRT_pvl))
  }
  
  
}



# ------------------------------------------------------------
# Fit a regression model for a continuous response.
#
# Supported models are lasso and random forest. The function returns
# fitted values and residuals on the training sample and, when supplied,
# predictions and residuals on a test sample.
# ------------------------------------------------------------

fit_gaussian <- function(Z_train, y_train, Z_test = NULL, y_test = NULL, 
                         model = 'linear', RF.num.trees = 100, lambda.seq = NULL,
                         CV = F){
  p <- ncol(Z_train)
  # Lasso regression with cross-validation for the penalty parameter.
  if (model == 'linear'){
    cv_lasso <- cv.glmnet(Z_train, y_train, alpha = 1, family = 'gaussian', 
                          lambda = lambda.seq, dfmax = as.integer(p / 2))
    lamb_opt <- cv_lasso$lambda.min
    fit_lasso <- glmnet(Z_train, y_train, alpha = 1, 
                        lambda = lamb_opt, family = 'gaussian')
    fit_beta <- fit_lasso$beta
    mean_pred_train <- as.vector(predict(fit_lasso, Z_train))
    res_train <- y_train - mean_pred_train
    
    mean_pred_test <- NULL
    res_test <- NULL
    if (! is.null(Z_test)){
      mean_pred_test <- as.vector(predict(fit_lasso, Z_test))
      if (! is.null(y_test)){
        res_test <- y_test - mean_pred_test
      }
    }
    return(list(pred_train = mean_pred_train, pred_test = mean_pred_test,
                res_train = res_train, res_test = res_test, coef = as.vector(fit_beta)))
  }
  
  # Random-forest regression, optionally using internal cross-fitting
  # to obtain out-of-fold training predictions.
  if (model == 'RF'){
    Z_train <- as.data.frame(Z_train)
    Z_test <- as.data.frame(Z_test)
    
    if (CV == F){
      
      fit_rf <- randomForest(x = Z_train, y = y_train, ntree = RF.num.trees)
      mean_pred_train <- predict(fit_rf, Z_train)
      res_train <- y_train - mean_pred_train
      
      mean_pred_test <- NULL
      res_test <- NULL
      if (! is.null(Z_test)){
        mean_pred_test <- as.vector(predict(fit_rf, Z_test))
        if (! is.null(y_test)){
          res_test <- y_test - mean_pred_test
        }
      }
    }
    
    # Five-fold cross-fitting for the training predictions when CV = TRUE.
    if (CV == T){
      n <- length(y_train)
      K <- 5
      resample <- sample(1:n, n)
      mean_pred_train <- rep(0, n)
      
      if (! is.null(Z_test)){
        mean_pred_test <- 0
      }
      for (k in 1:K){
        if (k == K){
          index_set <- ((k - 1) * as.integer(n / K) + 1):n
        }else{
          index_set <- ((k - 1) * as.integer(n / K) + 1):(k * as.integer(n / K))
        }
        test_set <- resample[index_set]
        train_set <- setdiff(1:n, test_set)

        fit_rf <- randomForest(x = Z_train[train_set,], y = y_train[train_set],
                               ntree = RF.num.trees)
        
        mean_pred_train[test_set] <- as.vector(predict(fit_rf, Z_train[test_set,]))
        
        if (! is.null(Z_test)){
          mean_pred_test <- mean_pred_test + as.vector(predict(fit_rf, Z_test))
        }
      }
      res_train <- y_train - mean_pred_train
      res_test <- NULL
      if (! is.null(Z_test)){
        mean_pred_test <- mean_pred_test / K
        if (! is.null(y_test)){
          res_test <- y_test - mean_pred_test
        }
      }
    }
    
    return(list(pred_train = mean_pred_train, pred_test = mean_pred_test,
                res_train = res_train, res_test = res_test, model = fit_rf))
  }
  
}


# ------------------------------------------------------------
# Fit a regression model for a binary response.
#
# The current implementation uses logistic lasso and returns fitted
# probabilities and residuals for the training and optional test samples.
# ------------------------------------------------------------

fit_binary <- function(Z_train, y_train, Z_test = NULL, y_test = NULL, 
                       model = 'linear', RF.num.trees = 100, lambda.seq = NULL,
                       CV = F){
  p <- ncol(Z_train)
  if (model == 'linear'){
    cv_lasso <- cv.glmnet(Z_train, y_train, alpha = 1, family = 'binomial', 
                          lambda = lambda.seq, dfmax = as.integer(p / 2))
    lamb_opt <- cv_lasso$lambda.min
    fit_lasso <- glmnet(Z_train, y_train, alpha = 1, 
                        lambda = lamb_opt, family = 'binomial')
    fit_beta <- fit_lasso$beta
    mean_pred_train <- glogit(as.vector(predict(fit_lasso, Z_train)))
    res_train <- y_train - mean_pred_train
    
    mean_pred_test <- NULL
    res_test <- NULL
    if (! is.null(Z_test)){
      mean_pred_test <- glogit(as.vector(predict(fit_lasso, Z_test)))
      if (! is.null(y_test)){
        res_test <- y_test - mean_pred_test
      }
    }
    return(list(pred_train = mean_pred_train, pred_test = mean_pred_test,
                res_train = res_train, res_test = res_test, coef = as.vector(fit_beta)))
  }
  
}


# ------------------------------------------------------------
# Calibrate the continuous X model using low-dimensional information
# extracted from the fitted Y | Z regression.
#
# The calibration removes residual structure in X that can be explained
# by the outcome-informed summary g(Z).
# ------------------------------------------------------------

cal_fit_gaussian <- function(gZ_train, y_train, gZ_test, y_test, model = 'linear',
                             RF.num.trees = 30){
  if (model == 'linear'){
    model.fit <- lm(y_train ~ gZ_train)
    res_train <- y_train - as.vector(cbind(1, gZ_train) %*% as.vector(model.fit$coefficient))
    res_test <- y_test - as.vector(cbind(1, gZ_test) %*% as.vector(model.fit$coefficient))
    return(list(res_train = res_train, res_test = res_test, coef = as.vector(model.fit$coefficient)))
  }
  
  if (model == 'RF'){
    gZ_test <- as.data.frame(gZ_test)
    gZ_train <- as.data.frame(gZ_train)
    colnames(gZ_test) <- paste0('V', 1:ncol(gZ_test))
    colnames(gZ_train) <- paste0('V', 1:ncol(gZ_train))
    
    model.fit <- randomForest(x = gZ_train, y = y_train, 
                              ntree = RF.num.trees, mtry = ncol(gZ_train))
    res_train <- y_train - predict(model.fit, gZ_train)
    res_test <- y_test - predict(model.fit, gZ_test)
    
    return(list(res_train = res_train, res_test = res_test, model = model.fit))
  }
}


# ------------------------------------------------------------
# Calibrate fitted probabilities for binary X using the same
# outcome-informed low-dimensional summary.
# ------------------------------------------------------------

cal_fit_binary <- function(gZ_train, y_train, pred_train,
                           gZ_test, y_test, pred_test, model = 'linear',
                           RF.num.trees = 30){
  if (model == 'linear'){
    model.fit <- glm(y_train ~ gZ_train, offset = logit(pred_train), family = binomial())
    pred_test_cal <- glogit(logit(pred_test) +
                              as.vector(cbind(1, gZ_test) %*% as.vector(model.fit$coefficient))) 
    
    return(list(pred_test = pred_test_cal, coef = as.vector(model.fit$coefficient)))
  }
}




# ============================================================
# Test statistics and p-values for the dCRT procedures
# ============================================================

# ------------------------------------------------------------
# dCRT for continuous X and continuous Y.
#
# The standardized X residual is combined with the Y residual to form
# either the basic residual-product statistic or an interaction-augmented
# statistic using the selected covariates Z_sub.
# ------------------------------------------------------------
dCRT_res <- function(res_x, res_y, Z_sub, k = 0, d.interaction = F){
  # Standardize the fitted X residuals.
  res_x <- res_x / sd(res_x) 
  n <- length(res_x)
  if (d.interaction == F | k == 0){
    imp_obe <- abs(mean(res_x * res_y))
    emp_var <- mean(res_y^2) 
    pvl <- 2 * pnorm(- sqrt(n) * imp_obe / sqrt(emp_var))
  }else{
    # Interaction-augmented dCRT statistic using Z_sub * res_y terms.
    weight_inter <- 1 / sqrt(k)
    W_inter <- res_y
    
    for (l in 1:k){
      W_inter <- cbind(W_inter, Z_sub[,l] * res_y)
    }
    
    XTX <- t(cbind(1, Z_sub)) %*% cbind(1, Z_sub)
    XTX_inv <- solve(XTX)
    
    Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% res_x
    imp_obe <- sum(Z_dI^2) 
    
    WTW <- t(W_inter) %*% W_inter
    svd_WTW <- svd(WTW)
    root_WTW <- svd_WTW$u %*% diag(sqrt(svd_WTW$d)) %*% t(svd_WTW$v)
    lambda_W <- svd(root_WTW %*% XTX_inv %*% 
                      diag(c(1, rep(weight_inter^2, k))) %*% XTX_inv %*% root_WTW)$d
    pvl <- tryCatch({
      pvl <- imhof(imp_obe, lambda_W, epsabs = 1e-8)
      pvl <- as.vector(abs(pvl$Qq))
    }, error = function(e){return(0)})
  }
  
  return(pvl)
}


# ------------------------------------------------------------
# dCRT for binary X and continuous Y.
#
# Conditional resamples of X are drawn from the fitted Bernoulli
# probabilities, and the p-value is computed by Monte Carlo comparison
# with the observed statistic.
# ------------------------------------------------------------

dCRT_res_binary <- function(A_label, pred_x_lab_cal, 
                            res_y_lab, Z_sub, k = 0, d.interaction = F, M = 500){
  n <- length(A_label)
  
  res_x_lab_cal_upsample <- (A_label - pred_x_lab_cal) / sqrt(pred_x_lab_cal * (1 - pred_x_lab_cal))
  Z_sub_upsample <- Z_sub
  res_y_lab_upsample <- res_y_lab
  
  if (d.interaction == F | k == 0){
    imp_obe <- abs(mean(res_x_lab_cal_upsample * res_y_lab_upsample))
    
    imp_sample <- unlist(lapply(c(1:M), function(j){
      A_sample <- rbinom(n, 1, pred_x_lab_cal)
      res_x_lab_cal_upsample <- (A_sample - pred_x_lab_cal) / sqrt(pred_x_lab_cal * (1 - pred_x_lab_cal))
      abs(mean(res_x_lab_cal_upsample * res_y_lab_upsample))
    })) 
    pvl <- (1 + sum(imp_sample >= imp_obe)) / (1 + M)

  }else{
    weight_inter <- 1 / sqrt(k)
    W_inter <- res_y_lab_upsample
    for (l in 1:k){
      W_inter <- cbind(W_inter, Z_sub_upsample[,l] * res_y_lab_upsample)
    }
    XTX <- t(cbind(1, Z_sub_upsample)) %*% cbind(1, Z_sub_upsample)
    XTX_inv <- solve(XTX)
    
    Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% res_x_lab_cal_upsample
    imp_obe <- sum(Z_dI^2) 
    
    imp_sample <- unlist(lapply(c(1:M), function(j){
      A_sample <- rbinom(n, 1, pred_x_lab_cal)
      res_x_lab_cal_upsample <- (A_sample - pred_x_lab_cal) / sqrt(pred_x_lab_cal * (1 - pred_x_lab_cal))
      
      Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% res_x_lab_cal_upsample
      sum(Z_dI^2) 
    })) 
    pvl <- (1 + sum(imp_sample >= imp_obe)) / (1 + M)
  }
  return(pvl)
  
}

# ------------------------------------------------------------
# dCRT for binary X and binary Y.
#
# Logistic-regression coefficient statistics are recomputed after
# Monte Carlo resampling of X from its fitted conditional distribution.
# ------------------------------------------------------------

dCRT_binary_Y <- function(A_label, pred_x_lab, y_label, offset_label,
                          Z_sub, k = 0, d.interaction = F, M = 500){
  n <- length(A_label)
  if (d.interaction == F | k == 0){
    data <- A_label
    glm_obe <- glm(y_label ~ data, offset = offset_label, family = binomial())
    imp_obe <- abs(coef(glm_obe)[2])
    
    imp_sample <- unlist(lapply(c(1:M), function(j){
      A_sample <- rbinom(n, 1, pred_x_lab)
      data_sample <- A_sample
      glm_sample <- glm(y_label ~ data_sample, offset = offset_label, family = binomial())
      abs(coef(glm_sample)[2])
    })) 
    pvl <- (1 + sum(imp_sample >= imp_obe)) / (1 + M)
    
  }else{
    data <- cbind(A_label, A_label * Z_sub, Z_sub)
    glm_obe <- glm(y_label ~ data, offset = offset_label, family = binomial())
    imp_obe <- (coef(glm_obe)[2])^2 + 1 / k * sum((coef(glm_obe)[3:(k+2)])^2)
  
    imp_sample <- unlist(lapply(c(1:M), function(j){
      A_sample <- rbinom(n, 1, pred_x_lab)
      data_sample <- cbind(A_sample, A_sample * Z_sub, Z_sub)
      glm_sample <- glm(y_label ~ data_sample, 
                        offset = offset_label, family = binomial())
      (coef(glm_sample)[2])^2 + 1 / k * sum((coef(glm_sample)[3:(k + 2)])^2)
    })) 
    pvl <- (1 + sum(imp_sample >= imp_obe)) / (1 + M)
  }
  return(pvl)
}


# ------------------------------------------------------------
# Nonparametric dCRT for continuous X and Y.
#
# The default random-forest version uses the importance of the X residual
# as the observed statistic and compares it with statistics obtained after
# Gaussian resampling of the standardized X residual.
# ------------------------------------------------------------

dCRT_res_nonp <- function(res_x, res_y, Z_sub, 
                          method = 'RF', RF.num.trees = 30, M = 500){
  
  res_x <- res_x / sd(res_x)
  n <- length(res_x)
  
  if (method == 'RF'){
    data <- as.data.frame(cbind(res_x, Z_sub))
    model.fit <- randomForest(x = data, y = res_y, ntree = RF.num.trees, mtry = ncol(data))
    imp_obe <- model.fit$importance[1]
    
    res_x_sample <- mvrnorm(M, rep(0, n), diag(rep(1, n)))
    imp_sample <- unlist(lapply(c(1:M), function(j){
      res_x_j <- res_x_sample[j,]
      data_sample <- as.data.frame(cbind(res_x_j, Z_sub))
      model.fit <- randomForest(x = data_sample, y = res_y, ntree = RF.num.trees, mtry = ncol(data))
      model.fit$importance[1]
      })) 
  }
  
  # Optional kernel-regression importance statistic.
  if (method == 'kernel'){
    kern.fit <- ksmooth(x = res_x, y = res_y,
                        kernel = "normal", bandwidth = n^(-1/5), x.points = res_x)
    
    imp_obe <- - sd(res_y[sort(res_x, decreasing = F, index.return = T)$ix] - kern.fit$y)
    res_x_sample <- mvrnorm(M, rep(0, n), diag(rep(1, n)))
    imp_sample <- unlist(lapply(c(1:M), function(j){
      res_x_j <- res_x_sample[j,]
      kern.fit.j <- ksmooth(x = res_x_j, y = res_y,
                            kernel = "normal", bandwidth = n^(-1/5), x.points = res_x_j)
      - sd(res_y[sort(res_x_j, decreasing = F, index.return = T)$ix] - kern.fit.j$y)
    }))
  }
  
  if (method == 'rank'){
    
  }
  
  pvl <- (1 + sum(imp_sample >= imp_obe)) / (1 + M)
  return(pvl)
}



# ============================================================
# Auxiliary functions
# ============================================================

# Combine selected covariates with the component of g(Z) orthogonal
# to their linear span.
orthogonalize <- function(gZ, Z_sub){
  model.fit <- lm(gZ ~ Z_sub)
  Z_convert <- cbind(Z_sub, model.fit$residual)
  return(Z_convert)
}

# Construct an autoregressive covariance matrix with
# Cov(Z_j, Z_k) = ar^{|j-k|}.
AR_cov <- function(p, ar){
  ar_series <- ar^(c(1:p) - 1)
  cov_mat <- ar_series
  for (i in 1:(p - 1)){
    cov_mat <- cbind(cov_mat, ar_series[c((p - i + 1):p ,1:(p - i))])
  }
  for (i in 1:(p - 1)){
    for (j in (i + 1):p){
      cov_mat[i, j] <- cov_mat[j, i]
    }
  }
  rownames(cov_mat) <- colnames(cov_mat)
  return(cov_mat)
}

# Logistic and logit transformations.
glogit <- function(x){
  return(1 / (1 + exp(-x)))
}

logit <- function(x){
  return(log(x / (1 - x)))
}



