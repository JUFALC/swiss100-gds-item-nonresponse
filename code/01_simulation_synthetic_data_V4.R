#### V4: Holm brackets removed, red zero line, optional RF-MI drop, results saved (.rds)
#### MI rf added

## -----------------------------------------------------------
## Synthetic GDS-15 simulation study:
## MCAR, MAR(dep + MMSE), MbD (L1 short-form vs L2/L3 long-form)
## Replacement: Complete-case, Rule-8, Sample Mean Sub, Ipsative,
##              MI (logistic FCS), MI_rf (RF-FCS),
##              MI4, MI4_rf for GDS-4
## -----------------------------------------------------------

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(mice)
  library(randomForest)  # for method = "rf" in mice
  library(ggplot2)
  library(viridisLite)
  library(ggstatsplot)
})

set.seed(2025)

## --------------------------
## Synthetic MMSE (0–21) + MbD by level
## --------------------------
## Levels by MMSE cutoffs:
##   L0: 0–5, L1: 6–9, L2: 10–14, L3: 15–21
## --------------------------

.mmse_bounds <- c(-Inf, 5.5, 9.5, 14.5, 21.5)  # continuity-corrected boundaries

fit_mmse_normal_to_levels <- function(target_props = c(0.3863, 0.1372, 0.1047, 0.3718),
                                      start_mu = 12, start_sd = 5) {
  stopifnot(abs(sum(target_props) - 1) < 1e-6, length(target_props) == 4)
  loss <- function(par) {
    mu <- par[1]
    sd <- exp(par[2])
    pk <- diff(pnorm(.mmse_bounds, mean = mu, sd = sd))  # length 4
    sum((pk - target_props)^2)
  }
  o <- optim(c(start_mu, log(start_sd)), loss, method = "Nelder-Mead",
             control = list(reltol = 1e-12, maxit = 5000))
  list(mu = o$par[1], sd = exp(o$par[2]), conv = o$convergence, value = o$value)
}

gen_mmse_from_theta <- function(theta, mu, sd, rho = -0.35) {
  z_th <- scale(theta)[, 1]
  z_ep <- rnorm(length(theta))
  z    <- rho * z_th + sqrt(pmax(0, 1 - rho^2)) * z_ep
  mmse_cont <- mu + sd * z
  mmse_int  <- pmin(21, pmax(0, round(mmse_cont)))
  mmse_int
}

gen_synthetic_optionA <- function(N = 277,
                                  K = 15,
                                  common4_idx = c(3,4,8,14),
                                  theta_sd = 1,
                                  slope_mean = 1.2, slope_sd = 0.3,
                                  target_item_ps = seq(0.10, 0.35, length.out = K),
                                  target_level_props = c(0.3863, 0.1372, 0.1047, 0.3718),
                                  mmse_rho_theta = -0.35,
                                  seed = 123) {
  set.seed(seed)
  
  fit <- fit_mmse_normal_to_levels(target_level_props, start_mu = 12, start_sd = 5)
  if (fit$conv != 0) warning("MMSE fitting did not fully converge; result used anyway.")
  mu_mmse <- fit$mu
  sd_mmse <- fit$sd
  
  theta <- rnorm(N, 0, theta_sd)
  
  beta  <- pmax(0.6, rnorm(K, slope_mean, slope_sd))
  alpha <- qlogis(target_item_ps)
  invlogit <- function(x) 1/(1+exp(-x))
  eta  <- outer(theta, beta) + matrix(rep(alpha, each = N), nrow = N)
  pmat <- invlogit(eta)
  X    <- matrix(rbinom(N*K, 1, c(pmat)), nrow = N, ncol = K)
  colnames(X) <- paste0("gds", 1:K)
  
  mmse <- gen_mmse_from_theta(theta, mu_mmse, sd_mmse, rho = mmse_rho_theta)
  
  level <- cut(
    mmse,
    breaks = c(-Inf, 5, 9, 14, 21),
    labels = c("L0","L1","L2","L3"),
    right  = TRUE
  )
  
  X_masked <- X
  keep4 <- common4_idx
  for (i in seq_len(N)) {
    if (level[i] == "L0") {
      X_masked[i, ] <- NA
    } else if (level[i] == "L1") {
      X_masked[i, -keep4] <- NA
    }
  }
  
  lev_share <- prop.table(table(level))
  names(lev_share) <- paste0("obs_", names(lev_share))
  target_props <- setNames(target_level_props, c("L0","L1","L2","L3"))
  
  list(
    X_truth  = X,
    X_masked = X_masked,
    theta    = theta,
    mmse     = mmse,
    level    = factor(level, levels = c("L0","L1","L2","L3")),
    cuts     = c("L0:0–5","L1:6–9","L2:10–14","L3:15–21"),
    mmse_fit = list(mu = mu_mmse, sd = sd_mmse, target = target_props, observed = lev_share)
  )
}

## ---- Generate synthetic data ----
synA <- gen_synthetic_optionA(N = 500, seed = 2025)

keep <- synA$mmse >= 6
X_truth_keep <- synA$X_truth[keep, , drop = FALSE]
X_obs        <- synA$X_masked[keep, , drop = FALSE]
truth_total  <- rowSums(X_truth_keep)
mmse         <- synA$mmse[keep]
level        <- droplevels(synA$level[keep])

idx4 <- c(3, 4, 8, 14)
X_truth4_keep <- X_truth_keep[, idx4, drop = FALSE]
X_obs4        <- X_obs[, idx4, drop = FALSE]

truth_total4  <- rowSums(X_truth4_keep)
truth_dep15   <- as.integer(truth_total >= 5)

.add_iter_row <- function(lst, iter, miss_lab, method_key, em, scenario = NA_character_) {
  tb <- tibble::tibble(
    iter       = iter,
    miss       = miss_lab,
    method     = method_key,
    MuHat      = as.numeric(em["MuHat"]),
    Bias       = as.numeric(em["Bias"]),
    MuTrue     = as.numeric(em["MuTrue"]),
    UsableRate = as.numeric(em["UsableRate"]),
    scenario   = scenario
  )
  lst[[length(lst) + 1L]] <- tb
  lst
}

## --------------------------
## 1) Scorers 
## --------------------------
score_complete_case <- function(X) {
  ok <- rowSums(!is.na(X)) == ncol(X)
  out <- rowSums(X, na.rm = TRUE)
  out[!ok] <- NA_real_
  out
}

score_rule8 <- function(X) {
  nonmiss <- rowSums(!is.na(X))
  s <- rowSums(X, na.rm = TRUE)
  s[nonmiss < 8] <- NA_real_
  s
}

score_mean_sub <- function(X) {
  cm <- colMeans(X, na.rm = TRUE)
  X2 <- X
  for (j in seq_len(ncol(X2))) {
    miss <- is.na(X2[, j])
    if (any(miss)) X2[miss, j] <- cm[j]
  }
  round(rowSums(X2), 1)
}

score_ipsative_x15 <- function(X) {
  round(rowMeans(X, na.rm = TRUE) * ncol(X), 1)
}

as_binary_factor_df <- function(X) {
  df <- as.data.frame(X)
  for (j in seq_len(ncol(df))) {
    xj <- df[[j]]
    if (!all(is.na(xj) | xj %in% c(0,1))) {
      stop(sprintf("Column %s has values outside {0,1,NA}", colnames(df)[j]))
    }
    df[[j]] <- factor(xj, levels = c(0,1))
  }
  df
}

is_constant_binary <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 1
}

## --- MI: default (logreg) ---
score_mi <- function(X, AUX = NULL, m = 20, seed = NULL, maxit = 10) {
  if (!is.null(seed)) set.seed(seed)
  
  items <- as_binary_factor_df(X)
  cn <- colnames(items)
  
  const <- vapply(items, is_constant_binary, logical(1))
  meth_items <- ifelse(const, "", "logreg")
  
  if (!is.null(AUX)) {
    dat <- cbind(items, AUX = as.numeric(AUX))
    meth <- c(setNames(meth_items, cn), AUX = "")
  } else {
    dat <- items
    meth <- setNames(meth_items, cn)
  }
  
  pred <- make.predictorMatrix(dat)
  diag(pred) <- 0
  if (!is.null(AUX)) {
    pred["AUX", ] <- 0
    pred[, "AUX"] <- 1
  }
  
  imp <- mice(dat, m = m, maxit = maxit, method = meth,
              predictorMatrix = pred, printFlag = FALSE)
  
  comp_sums <- sapply(1:m, function(k) {
    comp <- complete(imp, k)[, cn, drop = FALSE]
    comp_num <- as.data.frame(lapply(comp, function(col) as.numeric(as.character(col))))
    rowSums(comp_num)
  })
  round(rowMeans(comp_sums), 1)
}

## --- MI_rf: random forest FCS ---
score_mi_rf <- function(X, AUX = NULL, m = 20, seed = NULL, maxit = 10) {
  if (!is.null(seed)) set.seed(seed)
  
  items <- as_binary_factor_df(X)
  cn <- colnames(items)
  
  const <- vapply(items, is_constant_binary, logical(1))
  meth_items <- ifelse(const, "", "rf")
  
  if (!is.null(AUX)) {
    dat <- cbind(items, AUX = as.numeric(AUX))
    meth <- c(setNames(meth_items, cn), AUX = "")
  } else {
    dat <- items
    meth <- setNames(meth_items, cn)
  }
  
  pred <- make.predictorMatrix(dat)
  diag(pred) <- 0
  if (!is.null(AUX)) {
    pred["AUX", ] <- 0
    pred[, "AUX"] <- 1
  }
  
  imp <- mice(dat, m = m, maxit = maxit, method = meth,
              predictorMatrix = pred, printFlag = FALSE)
  
  comp_sums <- sapply(1:m, function(k) {
    comp <- complete(imp, k)[, cn, drop = FALSE]
    comp_num <- as.data.frame(lapply(comp, function(col) as.numeric(as.character(col))))
    rowSums(comp_num)
  })
  round(rowMeans(comp_sums), 1)
}

apply_all_methods <- function(X_obs, use_mi = FALSE, aux_vec = NULL, seed_mi = NULL) {
  out <- cbind(
    CompleteCase = score_complete_case(X_obs),
    Rule8        = score_rule8(X_obs),
    MeanSub      = score_mean_sub(X_obs),
    IpsativeX15  = score_ipsative_x15(X_obs)
  )
  if (use_mi) {
    mi_score    <- score_mi(X_obs,    AUX = aux_vec, m = 20, seed = seed_mi)
    mi_rf_score <- score_mi_rf(X_obs, AUX = aux_vec, m = 20, seed = seed_mi)
    out <- cbind(out, MI = mi_score, MI_rf = mi_rf_score)
  }
  out
}

## --- GDS-4 scorers ---
score_cc4 <- function(X4) {
  ok <- rowSums(!is.na(X4)) == ncol(X4)
  s  <- rowSums(X4, na.rm = TRUE)
  s[!ok] <- NA_real_
  s
}

score_ipsative_x4 <- function(X4) {
  round(rowMeans(X4, na.rm = TRUE) * ncol(X4), 1)
}

score_mi4 <- function(X4, AUX = NULL, m = 20, seed = NULL, maxit = 10) {
  score_mi(X4, AUX = AUX, m = m, seed = seed, maxit = maxit)
}

## --- MI4_rf for GDS-4 ---
score_mi4_rf <- function(X4, AUX = NULL, m = 20, seed = NULL, maxit = 10) {
  score_mi_rf(X4, AUX = AUX, m = m, seed = seed, maxit = maxit)
}

apply_methods_4 <- function(X4, use_mi = TRUE, aux_vec = NULL, seed_mi = NULL) {
  out <- cbind(
    CC4        = score_cc4(X4),
    IpsativeX4 = score_ipsative_x4(X4)
  )
  if (use_mi) {
    mi4    <- score_mi4(X4,    AUX = aux_vec, m = 20, seed = seed_mi)
    mi4_rf <- score_mi4_rf(X4, AUX = aux_vec, m = 20, seed = seed_mi)
    out <- cbind(out, MI4 = mi4, MI4_rf = mi4_rf)
  }
  out
}

## --------------------------
## 2) Metrics
## --------------------------
eval_metrics_mean <- function(est, truth) {
  ok <- !is.na(est)
  usable <- mean(ok)
  mu_true <- mean(truth)
  if (!any(ok)) {
    return(c(Bias = NA, RMSE = NA, MAE = NA, UsableRate = 0,
             MuHat = NA, MuTrue = mu_true))
  }
  mu_hat <- mean(est[ok])
  err <- est[ok] - mu_true
  c(Bias = mu_hat - mu_true,
    RMSE = sqrt(mean(err^2)),
    MAE  = mean(abs(err)),
    UsableRate = usable,
    MuHat = mu_hat,
    MuTrue = mu_true)
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

eval_metrics_mean4 <- function(est4, truth4) {
  ok <- !is.na(est4)
  mu_true <- safe_mean(truth4)
  if (!any(ok)) {
    return(c(Bias = NA, RMSE = NA, MAE = NA, UsableRate = 0,
             MuHat = NA, MuTrue = mu_true))
  }
  mu_hat <- mean(est4[ok])
  err    <- est4[ok] - mu_true
  c(Bias = mu_hat - mu_true,
    RMSE = sqrt(mean(err^2)),
    MAE  = mean(abs(err)),
    UsableRate = mean(ok),
    MuHat = mu_hat,
    MuTrue = mu_true)
}

auc_from_score <- function(y_true, score) {
  ok <- !is.na(score) & !is.na(y_true)
  y  <- y_true[ok]
  s  <- score[ok]
  if (length(unique(y)) != 2) return(NA_real_)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r  <- rank(s, ties.method = "average")
  auc <- (sum(r[y == 1]) - n1*(n1 + 1)/2) / (n1 * n0)
  as.numeric(auc)
}

eval_metrics_classif4 <- function(est4, truth_dep15, cutoff = 2) {
  ok <- !is.na(est4)
  if (!any(ok)) return(c(Sens = NA, Spec = NA, AUROC = NA, CalSlope = NA))
  y  <- truth_dep15[ok]
  s  <- est4[ok]
  yhat <- as.integer(s >= cutoff)
  
  tp <- sum(yhat == 1 & y == 1)
  fn <- sum(yhat == 0 & y == 1)
  tn <- sum(yhat == 0 & y == 0)
  fp <- sum(yhat == 1 & y == 0)
  sens <- if ((tp + fn) == 0) NA_real_ else tp/(tp+fn)
  spec <- if ((tn + fp) == 0) NA_real_ else tn/(tn+fp)
  
  auc <- auc_from_score(y, s)
  
  fit <- tryCatch(glm(y ~ s, family = binomial()), error = function(e) NULL)
  slope <- if (is.null(fit)) NA_real_ else unname(coef(fit)[2])
  
  c(Sens = sens, Spec = spec, AUROC = auc, CalSlope = slope)
}

## --------------------------
## 3) Missingness mechanisms
## --------------------------
logit <- function(x) 1/(1+exp(-x))

make_MCAR <- function(X, p_miss) {
  mask <- matrix(rbinom(length(X), 1, p_miss) == 1, nrow = nrow(X))
  Xm <- X
  Xm[mask] <- NA
  Xm
}

calibrate_intercept <- function(s, mmse, beta_s, beta_mmse, target) {
  z_s  <- scale(s)[, 1]
  z_mm <- scale(mmse)[, 1]
  f <- function(a0) mean(logit(a0 + beta_s * z_s + beta_mmse * z_mm))
  uniroot(function(a0) f(a0) - target, c(-8, 8))$root
}

make_MAR <- function(X, total, mmse,
                     target_miss = 0.20,
                     beta_s = 1.0, beta_mmse = -0.8) {
  a0 <- calibrate_intercept(total, mmse, beta_s, beta_mmse, target_miss)
  Xout <- X
  p_i <- logit(a0 + beta_s * scale(total)[, 1] + beta_mmse * scale(mmse)[, 1])
  for (i in seq_len(nrow(X))) {
    drop_i <- rbinom(ncol(X), 1, p_i[i]) == 1
    if (any(drop_i)) Xout[i, drop_i] <- NA
  }
  Xout
}

## --------------------------
## 4) Simulation engines
## --------------------------

## Part A: L2/L3 only (GDS-15)
run_sims_A_L23 <- function(p_mech = c("MCAR","MAR","MbD"),
                           ns = 200,
                           miss_rates = c(0.10, 0.20, 0.30),
                           seed = 1234,
                           beta_s = 1.0, beta_mmse = -0.8,
                           use_mi = TRUE) {
  
  p_mech <- match.arg(p_mech)
  set.seed(seed)
  
  METHODS <- c("CompleteCase","Rule8","MeanSub","IpsativeX15",
               if (use_mi) c("MI","MI_rf"))
  METRICS <- c("Bias","UsableRate","MuHat","MuTrue")
  
  res <- array(NA_real_,
               dim = c(ns, length(METRICS), length(METHODS), length(miss_rates)),
               dimnames = list(iter  = 1:ns,
                               metric = METRICS,
                               method = METHODS,
                               miss   = paste0("p", miss_rates)))
  
  idx    <- level %in% c("L2","L3")
  X0     <- X_obs[idx, , drop = FALSE]
  truth  <- truth_total[idx]
  aux_mm <- mmse[idx]
  
  iter_rows <- list()
  
  for (m in seq_along(miss_rates)) {
    pm <- miss_rates[m]
    miss_lab <- paste0("p", pm)
    
    for (it in 1:ns) {
      Xmiss <- switch(
        p_mech,
        MbD  = X0,
        MCAR = make_MCAR(X0, pm),
        MAR  = make_MAR(X0, total = truth, mmse = aux_mm,
                        target_miss = pm, beta_s = beta_s, beta_mmse = beta_mmse)
      )
      
      est <- apply_all_methods(Xmiss, use_mi = use_mi, aux_vec = aux_mm,
                               seed_mi = sample.int(1e6, 1))
      est <- est[, METHODS, drop = FALSE]
      
      for (k in seq_along(METHODS)) {
        em <- eval_metrics_mean(est[, k], truth)
        res[it, , k, m] <- c(
          Bias       = em["Bias"],
          UsableRate = em["UsableRate"],
          MuHat      = em["MuHat"],
          MuTrue     = em["MuTrue"]
        )
        iter_rows <- .add_iter_row(iter_rows, iter = it, miss_lab = miss_lab,
                                   method_key = METHODS[k], em = em,
                                   scenario = p_mech)
      }
    }
  }
  
  attr(res, "iter_df") <- dplyr::bind_rows(iter_rows)
  res
}

## Part B: L1–L3 combined (GDS-15)
run_sims_B_all <- function(p_mech = c("MCAR","MAR","MbD"),
                           ns = 200,
                           miss_rates = c(0.10, 0.20, 0.30),
                           seed = 1234,
                           beta_s = 1.0, beta_mmse = -0.8,
                           use_mi = TRUE) {
  
  p_mech <- match.arg(p_mech)
  set.seed(seed)
  
  METHODS <- c("MeanSub","IpsativeX15", if (use_mi) c("MI","MI_rf"))
  METRICS <- c("Bias","UsableRate","MuHat","MuTrue")
  
  res <- array(NA_real_,
               dim = c(ns, length(METRICS), length(METHODS), length(miss_rates)),
               dimnames = list(iter  = 1:ns,
                               metric = METRICS,
                               method = METHODS,
                               miss   = paste0("p", miss_rates)))
  
  X0     <- X_obs
  truth  <- truth_total
  aux_mm <- mmse
  
  iter_rows <- list()
  
  for (m in seq_along(miss_rates)) {
    pm <- miss_rates[m]
    miss_lab <- paste0("p", pm)
    
    for (it in 1:ns) {
      Xmiss <- switch(
        p_mech,
        MbD  = X0,
        MCAR = make_MCAR(X0, pm),
        MAR  = make_MAR(X0, total = truth, mmse = aux_mm,
                        target_miss = pm, beta_s = beta_s, beta_mmse = beta_mmse)
      )
      
      est_all <- apply_all_methods(Xmiss, use_mi = use_mi, aux_vec = aux_mm,
                                   seed_mi = sample.int(1e6, 1))
      est <- est_all[, METHODS, drop = FALSE]
      
      for (k in seq_along(METHODS)) {
        em <- eval_metrics_mean(est[, k], truth)
        res[it, , k, m] <- c(
          Bias       = em["Bias"],
          UsableRate = em["UsableRate"],
          MuHat      = em["MuHat"],
          MuTrue     = em["MuTrue"]
        )
        iter_rows <- .add_iter_row(iter_rows, iter = it, miss_lab = miss_lab,
                                   method_key = METHODS[k], em = em,
                                   scenario = p_mech)
      }
    }
  }
  
  attr(res, "iter_df") <- dplyr::bind_rows(iter_rows)
  res
}

## Part B4: L1–L3 combined, GDS-4 only
run_sims_B4_only <- function(p_mech = c("MCAR","MAR","MbD"),
                             ns = 200,
                             miss_rates = c(0.10, 0.20, 0.30),
                             seed = 1234,
                             beta_s = 1.0, beta_mmse = -0.8,
                             use_mi = TRUE) {
  
  p_mech <- match.arg(p_mech)
  set.seed(seed)
  
  METHODS <- c("CC4","IpsativeX4", if (use_mi) c("MI4","MI4_rf"))
  METRICS_MEAN <- c("Bias","RMSE","MAE","UsableRate","MuHat","MuTrue")
  METRICS_CLF  <- c("Sens","Spec","AUROC","CalSlope")
  
  res_mean <- array(NA_real_,
                    dim = c(ns, length(METRICS_MEAN), length(METHODS), length(miss_rates)),
                    dimnames = list(iter = 1:ns, metric = METRICS_MEAN,
                                    method = METHODS, miss = paste0("p", miss_rates)))
  
  res_clf  <- array(NA_real_,
                    dim = c(ns, length(METRICS_CLF), length(METHODS), length(miss_rates)),
                    dimnames = list(iter = 1:ns, metric = METRICS_CLF,
                                    method = METHODS, miss = paste0("p", miss_rates)))
  
  X0     <- X_obs4
  truth4 <- truth_total4
  y15    <- truth_dep15
  aux_mm <- mmse
  
  iter_rows <- list()
  
  for (m in seq_along(miss_rates)) {
    pm <- miss_rates[m]
    miss_lab <- paste0("p", pm)
    
    for (it in 1:ns) {
      Xmiss4 <- switch(
        p_mech,
        MbD  = X0,
        MCAR = make_MCAR(X0, pm),
        MAR  = make_MAR(X0, total = truth_total, mmse = aux_mm,
                        target_miss = pm, beta_s = beta_s, beta_mmse = beta_mmse)
      )
      
      est4 <- apply_methods_4(Xmiss4, use_mi = use_mi, aux_vec = aux_mm,
                              seed_mi = sample.int(1e6, 1))
      est4 <- est4[, METHODS, drop = FALSE]
      
      for (k in seq_along(METHODS)) {
        em <- eval_metrics_mean4(est4[, k], truth4)
        res_mean[it, , k, m] <- em
        
        ec <- eval_metrics_classif4(est4[, k], y15, cutoff = 2)
        res_clf[it, , k, m] <- ec
        
        iter_rows <- .add_iter_row(iter_rows, iter = it, miss_lab = miss_lab,
                                   method_key = METHODS[k], em = em,
                                   scenario = p_mech)
      }
    }
  }
  
  attr(res_mean, "iter_df") <- dplyr::bind_rows(iter_rows)
  list(mean = res_mean, clf = res_clf)
}

## --------------------------
## 5) Run scenarios
## --------------------------
resA_MCAR <- run_sims_A_L23("MCAR", ns = 100, miss_rates = 0.20,
                            use_mi = TRUE, seed = 2025)
resA_MAR  <- run_sims_A_L23("MAR",  ns = 100, miss_rates = 0.20,
                            use_mi = TRUE, seed = 2025,
                            beta_s = 1.1, beta_mmse = -0.9)

resB_MbD  <- run_sims_B_all("MbD",  ns = 100, miss_rates = 0.00,
                            use_mi = TRUE, seed = 2025)
resB_MCAR <- run_sims_B_all("MCAR", ns = 100, miss_rates = 0.20,
                            use_mi = TRUE, seed = 2025)
resB_MAR  <- run_sims_B_all("MAR",  ns = 100, miss_rates = 0.20,
                            use_mi = TRUE, seed = 2025,
                            beta_s = 1.1, beta_mmse = -0.9)

resB4_MbD  <- run_sims_B4_only("MbD",  ns = 100, miss_rates = 0.00,
                               use_mi = TRUE, seed = 2025)
resB4_MCAR <- run_sims_B4_only("MCAR", ns = 100, miss_rates = 0.20,
                               use_mi = TRUE, seed = 2025)
resB4_MAR  <- run_sims_B4_only("MAR",  ns = 100, miss_rates = 0.20,
                               use_mi = TRUE, seed = 2025,
                               beta_s = 1.1, beta_mmse = -0.9)

## V4: save raw results -> re-plot later with readRDS() instead of re-simulating
saveRDS(list(resA_MCAR=resA_MCAR, resA_MAR=resA_MAR,
             resB_MbD=resB_MbD, resB_MCAR=resB_MCAR, resB_MAR=resB_MAR,
             resB4_MbD=resB4_MbD, resB4_MCAR=resB4_MCAR, resB4_MAR=resB4_MAR),
        "sim_results_V4.rds")

## --------------------------
## 6) Summaries
## --------------------------
summarize_results <- function(res) {
  stopifnot(
    is.array(res),
    identical(names(dimnames(res)), c("iter","metric","method","miss"))
  )
  bias   <- res[, "Bias",       , , drop = FALSE]
  usable <- res[, "UsableRate", , , drop = FALSE]
  
  Bias_mean <- apply(bias,        c(3, 4), mean, na.rm = TRUE)
  RMSE      <- sqrt(apply(bias^2, c(3, 4), mean, na.rm = TRUE))
  MAE       <- apply(abs(bias),   c(3, 4), mean, na.rm = TRUE)
  Use_mean  <- apply(usable,      c(3, 4), mean, na.rm = TRUE)
  
  metrics <- c("Bias","RMSE","MAE","UsableRate")
  out <- array(NA_real_,
               dim = c(length(dimnames(res)$method), length(metrics),
                       length(dimnames(res)$miss)),
               dimnames = list(method = dimnames(res)$method,
                               metric = metrics,
                               miss   = dimnames(res)$miss))
  out[, "Bias",       ] <- Bias_mean
  out[, "RMSE",       ] <- RMSE
  out[, "MAE",        ] <- MAE
  out[, "UsableRate", ] <- Use_mean
  out
}

print_summary <- function(res, digits = 3, title = deparse(substitute(res))) {
  cat("\n=== ", title, " ===\n", sep = "")
  out <- summarize_results(res)
  print(round(out, digits))
  invisible(out)
}

print_summary_if_exists <- function(obj_name, title, digits = 3) {
  if (exists(obj_name, inherits = TRUE)) {
    print_summary(get(obj_name), digits = digits, title = title)
  } else {
    cat("\n=== ", title, " ===\n", sep = "")
    cat("[skipped: object '", obj_name, "' not found]\n", sep = "")
  }
}

summarize_results_mean <- function(res_mean) summarize_results(res_mean)

summarize_results_clf <- function(res_clf) {
  stopifnot(
    is.array(res_clf),
    identical(names(dimnames(res_clf)), c("iter","metric","method","miss"))
  )
  out <- apply(res_clf, c(2,3,4), function(x) mean(x, na.rm = TRUE))
  aperm(out, c(2,1,3))
}

print_summary4 <- function(res_list, title_prefix = "", cutoff_gds4 = 2) {
  if (is.null(res_list)) { cat("[no results]\n"); return(invisible(NULL)) }
  
  cat("\n=== ", title_prefix, " — GDS-4 (primary estimand) ===\n", sep = "")
  cat("Primary estimand: depression classification relative to the true GDS-15 status (≥5), ",
      "treating GDS-4 ≥", cutoff_gds4, " as the index test. ",
      "We report sensitivity, specificity, AUROC, and calibration slope.\n\n",
      sep = "")
  
  cat("— GDS-4 mean performance (secondary estimand) —\n",
      "Bias, RMSE, MAE, and usable case proportion are reported against the true GDS-4 mean.\n\n", sep = "")
  
  cat(">>> GDS-4 Mean Metrics (vs true GDS-4) <<<\n")
  print(round(summarize_results_mean(res_list$mean), 3))
  
  cat("\n>>> GDS-4 Classification Metrics (vs true GDS-15 ≥5) <<<\n")
  print(round(summarize_results_clf(res_list$clf), 3))
  
  invisible(NULL)
}

## Print summaries
print_summary_if_exists("resA_MCAR", "A: L2/L3 — MCAR (p=0.20)")
print_summary_if_exists("resA_MAR",  "A: L2/L3 — MAR (p≈0.20; beta_s=1.1, beta_mmse=-0.9)")

print_summary_if_exists("resB_MbD",  "B: L1–L3 — MbD only")
print_summary_if_exists("resB_MCAR", "B: L1–L3 — MbD + MCAR (p=0.20)")
print_summary_if_exists("resB_MAR",  "B: L1–L3 — MbD + MAR (p≈0.20; beta_s=1.1, beta_mmse=-0.9)")

print_summary4(resB4_MbD,  "B4: L1–L3 — MbD only")
print_summary4(resB4_MCAR, "B4: L1–L3 — MbD + MCAR (p=0.20)")
print_summary4(resB4_MAR,  "B4: L1–L3 — MbD + MAR (p≈0.20; beta_s=1.1, beta_mmse=-0.9)")

## --------------------------
## 7) Visualization helpers
## --------------------------
to_metric_df <- function(res, metric) {
  tbl <- as.data.frame(as.table(res[, metric, , , drop = FALSE]))
  names(tbl) <- c("iter","metric_name","method","miss","value")
  tbl$metric_name <- NULL
  names(tbl)[names(tbl) == "value"] <- metric
  tbl
}

array_to_df <- function(res, scenario) {
  stopifnot(
    is.array(res),
    identical(names(dimnames(res)), c("iter","metric","method","miss"))
  )
  df <- Reduce(function(x, y) merge(x, y, by = c("iter","method","miss")),
               list(to_metric_df(res, "MuHat"),
                    to_metric_df(res, "Bias"),
                    to_metric_df(res, "MuTrue"),
                    to_metric_df(res, "UsableRate")))
  df$scenario <- scenario
  df
}

arrays_to_long <- function(named_list) {
  named_list <- Filter(Negate(is.null), named_list)
  if (length(named_list) == 0) return(data.frame())
  bind_rows(lapply(names(named_list), function(scn) array_to_df(named_list[[scn]], scn)))
}

nice_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(size = 0.3, linetype = "dotted"),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = "grey85"),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 20, hjust = 1),
    legend.position    = "none",
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(colour = "grey30")
  )

pretty_method <- function(x) {
  mp <- c(
    "CompleteCase" = "Complete Case",
    "Rule8"        = "Rule-8",
    "MeanSub"      = "Sample Mean Sub",
    "IpsativeX15"  = "Ipsative",
    "MI"           = "Multiple Imputation",
    "MI_rf"        = "MI (RF-FCS)",
    # ---- GDS-4 ----
    "CC4"          = "Complete Case",
    "IpsativeX4"   = "Ipsative",
    "MI4"          = "Multiple Imputation",
    "MI4_rf"       = "MI (RF-FCS)"
  )
  as.character(mp[as.character(x)])
}

plot_part <- function(df_long, part_title, one_column = FALSE, scale_name = "GDS-15") {
  if (nrow(df_long) == 0) {
    warning(sprintf("No data to plot for '%s' (df_long is empty).", part_title))
    return(list(mu = NULL, bias = NULL))
  }
  
  scn_levels <- intersect(c("MCAR","MAR","MbD"), unique(df_long$scenario))
  df_long$scenario <- factor(df_long$scenario, levels = scn_levels)
  df_long$miss     <- factor(df_long$miss, levels = sort(unique(df_long$miss)))
  
  df_long$method <- pretty_method(df_long$method)
  
  meth_vals <- unique(df_long$method)
  if (any(grepl("GDS-4", meth_vals))) {
    desired_order <- c("Complete Case","Ipsative","Multiple Imputation","MI (RF-FCS)")
  } else {
    desired_order <- c("Complete Case","Rule-8","Sample Mean Sub",
                       "Ipsative","Multiple Imputation","MI (RF-FCS)")
  }
  present_order  <- desired_order[desired_order %in% meth_vals]
  if (length(present_order) == 0L) present_order <- sort(meth_vals)
  df_long$method <- factor(df_long$method, levels = present_order, ordered = TRUE)
  
  if (one_column) {
    df_long$facet_lab <- interaction(df_long$scenario, df_long$miss, drop = TRUE, sep = " — ")
    df_truth <- df_long |>
      dplyr::group_by(facet_lab) |>
      dplyr::summarise(MuTrue = unique(MuTrue)[1], .groups = "drop")
  } else {
    df_truth <- df_long |>
      dplyr::group_by(scenario, miss) |>
      dplyr::summarise(MuTrue = unique(MuTrue)[1], .groups = "drop")
  }
  
  meth_levels <- levels(df_long$method)
  method_cols <- setNames(
    viridisLite::viridis(length(meth_levels), option = "D", begin = 0.1, end = 0.9),
    meth_levels
  )
  
  p_mu <- ggplot(df_long |> dplyr::filter(!is.na(MuHat)),
                 aes(x = method, y = MuHat, fill = method)) +
    geom_hline(data = df_truth,
               aes(yintercept = MuTrue),
               linetype = "dashed", linewidth = 0.4) +
    geom_boxplot(width = 0.65, outlier.alpha = 0.35,
                 alpha = 0.8, colour = "grey25") +
    stat_summary(fun = mean, geom = "point", shape = 21,
                 size = 2.5, fill = "white", colour = "black") +
    { if (one_column) facet_grid(facet_lab ~ ., scales = "free_y")
      else facet_grid(scenario ~ miss, scales = "free_y") } +
    scale_fill_manual(values = method_cols) +
    labs(x = NULL, y = sprintf("Estimated mean %s (MuHat)", scale_name),
         title = sprintf("%s: Summaries of %s scores", part_title, scale_name),
         subtitle = "Dashed line = MuTrue") +
    nice_theme
  
  p_bias <- ggplot(df_long |> dplyr::filter(!is.na(Bias)),
                   aes(x = method, y = Bias, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_boxplot(width = 0.65, outlier.alpha = 0.35,
                 alpha = 0.8, colour = "grey25") +
    stat_summary(fun = mean, geom = "point", shape = 21,
                 size = 2.5, fill = "white", colour = "black") +
    { if (one_column) facet_grid(facet_lab ~ ., scales = "free_y")
      else facet_grid(scenario ~ miss, scales = "free_y") } +
    scale_fill_manual(values = method_cols) +
    labs(x = NULL, y = "Bias = MuHat − MuTrue",
         title = sprintf("%s: Bias of %s scores", part_title, scale_name)) +
    nice_theme
  
  list(mu = p_mu, bias = p_bias)
}

## Build long data for panels
df_A <- arrays_to_long(list(
  MCAR = if (exists("resA_MCAR")) resA_MCAR else NULL,
  MAR  = if (exists("resA_MAR"))  resA_MAR  else NULL
))

df_B <- arrays_to_long(list(
  MbD  = if (exists("resB_MbD"))  resB_MbD  else NULL,
  MCAR = if (exists("resB_MCAR")) resB_MCAR else NULL,
  MAR  = if (exists("resB_MAR"))  resB_MAR  else NULL
))

df_C <- arrays_to_long(list(
  MbD  = if (exists("resB4_MbD"))  resB4_MbD$mean  else NULL,
  MCAR = if (exists("resB4_MCAR")) resB4_MCAR$mean else NULL,
  MAR  = if (exists("resB4_MAR"))  resB4_MAR$mean  else NULL
))

## Panel A (L2/L3 only, GDS-15)
plots_A <- plot_part(df_A, part_title = "Panel A (L2/L3 only), GDS-15", one_column = FALSE)
if (!is.null(plots_A$mu))  print(plots_A$mu)
if (!is.null(plots_A$bias)) print(plots_A$bias)

## Panel B (L1–L3 combined, GDS-15)
plots_B <- plot_part(df_B, part_title = "Panel B (L1–L3 combined), GDS-15", one_column = TRUE)
if (!is.null(plots_B$mu))  print(plots_B$mu)
if (!is.null(plots_B$bias)) print(plots_B$bias)

## Panel C (L1–L3 combined, GDS-4)
plots_C <- plot_part(df_C, part_title = "Panel C (L1–L3 combined)", scale_name = "GDS-4")
if (!is.null(plots_C$mu))   print(plots_C$mu)
if (!is.null(plots_C$bias)) print(plots_C$bias)

DROP_RF <- TRUE   # V4: TRUE = figures show the 5 methods of the manuscript; FALSE = add MI (RF-FCS)

######################################################
## 8) ggstatsplot figures (V4: no pairwise tests)
######################################################
prep_for_stats <- function(df_long, one_column = FALSE) {
  if (DROP_RF) df_long <- df_long[!as.character(df_long$method) %in% c("MI_rf","MI4_rf"), , drop = FALSE]
  df_long$scenario <- factor(
    df_long$scenario,
    levels = intersect(c("MCAR","MAR","MbD"), unique(df_long$scenario))
  )
  miss_chr <- as.character(df_long$miss)
  miss_num <- suppressWarnings(as.numeric(sub("^p", "", miss_chr)))
  df_long$miss <- factor(
    df_long$miss,
    levels = unique(miss_chr)[order(miss_num, na.last = TRUE)],
    ordered = TRUE
  )
  
  pretty <- c(
    "CompleteCase" = "Complete Case",
    "Rule8"        = "Rule-8",
    "MeanSub"      = "Sample Mean Sub",
    "IpsativeX15"  = "Ipsative",
    "MI"           = "Multiple Imputation",
    "MI_rf"        = "Multiple Imputation (RF-FCS)",
    # GDS-4
    "CC4"          = "Complete Case (GDS-4)",
    "IpsativeX4"   = "Ipsative (GDS-4)",
    "MI4"          = "Multiple Imputation (GDS-4)",
    "MI4_rf"       = "Multiple Imputation (GDS-4, RF-FCS)"
  )
  m <- as.character(df_long$method)
  m <- ifelse(m %in% names(pretty), pretty[m], m)
  wanted <- c(
    "Complete Case","Rule-8","Sample Mean Sub",
    "Ipsative","Multiple Imputation","Multiple Imputation (RF-FCS)",
    "Complete Case (GDS-4)","Ipsative (GDS-4)",
    "Multiple Imputation (GDS-4)","Multiple Imputation (GDS-4, RF-FCS)"
  )
  df_long$method <- factor(m, levels = wanted[wanted %in% unique(m)], ordered = TRUE)
  
  if (one_column) {
    df_long$facet_lab <- interaction(df_long$scenario, df_long$miss, drop = TRUE, sep = " — ")
  }
  df_long
}

bias_betweenstats <- function(df_long, title, one_column = FALSE, k_digits = 2, ncol_facets = NULL) {
  d <- prep_for_stats(df_long, one_column = one_column)
  d <- d[!is.na(d$Bias), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  
  grouping_var <- if (one_column) "facet_lab" else "scenario"
  n_groups <- length(unique(d[[grouping_var]]))
  ncol_use <- if (!is.null(ncol_facets)) ncol_facets else if (one_column) 1 else min(3, n_groups)
  
  old_opt <- options(scipen = 999)
  on.exit(options(old_opt), add = TRUE)
  
  ggstatsplot::grouped_ggbetweenstats(
    data                 = d,
    x                    = method,
    y                    = Bias,
    grouping.var         = !!rlang::sym(grouping_var),
    plotgrid.args        = list(ncol = ncol_use),
    results.subtitle     = FALSE,
    pairwise.display     = "none",   # V4: no Holm brackets (current ggstatsplot API)
    sample.size.label    = FALSE,    # V4: no "(n = 100)" under axis labels
    effsize.type         = "g",
    mean.plotting        = TRUE,
    mean.ci              = TRUE,
    k                    = k_digits,
    ggtheme              = ggplot2::theme_bw(base_size = 12),
    title.text           = title,
    xlab                 = NULL,
    ylab                 = "Bias",
    ggplot.component = list(
      ggplot2::geom_hline(yintercept = 0, colour = "red"),   # V4: no-bias reference line
      ggplot2::scale_x_discrete(
        labels = function(lbl) gsub("\\(n\\s*=\\s*([0-9]+)\\)", "ns = \\1", lbl)
      ),
      ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
    )
  )
}

muhat_betweenstats <- function(df_long, title, one_column = FALSE, k_digits = 2, ncol_facets = NULL) {
  d <- prep_for_stats(df_long, one_column = one_column)
  d <- d[!is.na(d$MuHat), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  
  grouping_var <- if (one_column) "facet_lab" else "scenario"
  n_groups <- length(unique(d[[grouping_var]]))
  ncol_use <- if (!is.null(ncol_facets)) ncol_facets else if (one_column) 1 else min(3, n_groups)
  
  old_opt <- options(scipen = 999)
  on.exit(options(old_opt), add = TRUE)
  
  ggstatsplot::grouped_ggbetweenstats(
    data                 = d,
    x                    = method,
    y                    = MuHat,
    grouping.var         = !!rlang::sym(grouping_var),
    plotgrid.args        = list(ncol = ncol_use),
    results.subtitle     = FALSE,
    pairwise.display     = "none",   # V4: no Holm brackets (current ggstatsplot API)
    sample.size.label    = FALSE,    # V4: no "(n = 100)" under axis labels
    effsize.type         = "g",
    mean.plotting        = TRUE,
    mean.ci              = TRUE,
    k                    = k_digits,
    ggtheme              = ggplot2::theme_bw(base_size = 12),
    title.text           = title,
    xlab                 = NULL,
    ylab                 = "GDS score",
    ggplot.component = list(
      ggplot2::scale_x_discrete(
        labels = function(lbl) gsub("\\(n\\s*=\\s*([0-9]+)\\)", "ns = \\1", lbl)
      ),
      ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
    )
  )
}

## Part A (L2/L3, GDS-15)
if (exists("df_A") && nrow(df_A)) {
  pA_bias <- bias_betweenstats(
    df_A,
    title      = "A (L2/L3): Bias by method within MCAR/MAR",
    one_column = FALSE
  )
  if (!is.null(pA_bias)) print(pA_bias)
  
  pA_muhat <- muhat_betweenstats(
    df_A,
    title      = "A (L2/L3): MuHat by method within MCAR/MAR",
    one_column = FALSE
  )
  if (!is.null(pA_muhat)) print(pA_muhat)
}

## Part B (L1–L3, GDS-15)
if (exists("df_B") && nrow(df_B)) {
  pB_bias <- bias_betweenstats(
    df_B,
    title      = "B (L1–L3): Bias by method — MbD | MCAR | MAR",
    one_column = FALSE,
    k_digits   = 3
  )
  if (!is.null(pB_bias)) print(pB_bias)
  
  pB_muhat <- muhat_betweenstats(
    df_B,
    title      = "B (L1–L3): MuHat by method — MbD | MCAR | MAR",
    one_column = FALSE,
    k_digits   = 3
  )
  if (!is.null(pB_muhat)) print(pB_muhat)
}

## Part C (L1–L3, GDS-4)
if (exists("df_C") && nrow(df_C)) {
  df_C_sub <- df_C[df_C$scenario != "MbD", , drop = FALSE]
  
  pC_bias <- bias_betweenstats(
    df_C_sub,
    title      = "C (L1–L3, GDS-4): Bias by method — MCAR | MAR",
    one_column = FALSE,
    k_digits   = 3
  )
  if (!is.null(pC_bias)) print(pC_bias)
  
  pC_muhat <- muhat_betweenstats(
    df_C_sub,
    title      = "C (L1–L3, GDS-4 mean): MuHat by method — MCAR | MAR",
    one_column = FALSE,
    k_digits   = 3
  )
  if (!is.null(pC_muhat)) print(pC_muhat)
}
