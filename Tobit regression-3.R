#********************************************************************************
#* Model:  tobit model
#* Author: Abdulbaki Bilgic
#*******************************************************************************
#* Clear environment and graphics:
rm(list = ls()); cat("\f"); graphics.off()
#*******************************************************************************
#* load required libraries:
library(maxLik); library(AER)
#*******************************************************************************
# --- Log-likelihood function --------------------------------------------------
gaussian_tobit_ll <- function(bs, data) {
  np <- length(bs)
  nx <- np - 1
  b  <- bs[1:nx]
  s  <- bs[np]
  
  # censored observations
  x  <- data[data[, 1] == 0, , drop = FALSE]
  f0 <- if(nrow(x) > 0) pnorm((x[,1] - x[,2:np] %*% b)/s, log.p = TRUE) else numeric(0)
  
  # uncensored observations
  x  <- data[data[, 1] > 0, , drop = FALSE]
  f1 <- if(nrow(x) > 0) -log(s) + dnorm((x[,1] - x[,2:np] %*% b)/s, log=TRUE) else numeric(0)
  
  return(sum(c(f0, f1)))
}
# --- Custom summary function --------------------------------------------------
summary.tobitMaxLik <- function(result, data) {
  # Observations
  n_total    <- nrow(data)
  n_left     <- sum(data[,1] == 0)
  n_uncens   <- sum(data[,1] > 0)
  n_right    <- 0
  obs_summary<- c("Total"          = n_total,
                  "Left-censored"  = n_left,
                  "Uncensored"     = n_uncens,
                  "Right-censored" = n_right)
  cat("Observations:\n")
  print(obs_summary)
  
  # Estimates
  est  <- result$estimate
  se   <- sqrt(diag(vcov(result)))
  tval <- est / se
  pval <- 2 * pnorm(-abs(tval))
  
  # Significance codes
  signif_code      <- rep(" ", length(pval))
  signif_code[pval < 0.1]   <- "."
  signif_code[pval < 0.05]  <- "*"
  signif_code[pval < 0.01]  <- "**"
  signif_code[pval < 0.001] <- "***"
  
  pval_formatted <- format.pval(pval, digits = 3, eps = 0.001)
  
  coef_table <- data.frame(
    Estimate     = round(est,5),
    `Std. Error` = round(se,5),
    `t value`    = round(tval,3),
    `Pr(>|t|)`   = pval_formatted,
    Signif       = signif_code,
    row.names    = names(est),
    check.names  = FALSE
  )
  
  cat("\nCoefficients:\n")
  print(coef_table)
  
  # Wald statistic (excluding sigma)
  beta_hat  <- est[1:(length(est)-1)]
  vcov_mat  <- vcov(result)[1:(length(est)-1), 1:(length(est)-1), drop=FALSE]
  wald_stat <- t(beta_hat) %*% solve(vcov_mat) %*% beta_hat
  wald_stat <- as.numeric(wald_stat)
  df        <- length(beta_hat)
  p_val     <- 1 - pchisq(wald_stat, df)
  cat("\nWald-statistic:", round(wald_stat,2), "on", df, "Df, p-value:", 
      format.pval(p_val, digits=3, eps=0.001), "\n")
  
  # ---------------------------------------------------------------------------
  # Marginal Effects (at means) + Delta Method for SEs
  # ---------------------------------------------------------------------------
  X         <- data[, -1]                                                       # X matrix with intercept
  xbar      <- colMeans(X)                                                      # mean of regressors
  mu_hat    <- sum(xbar * beta_hat)
  sigma_hat <- est["sigma"]
  
  # Core components
  z   <- mu_hat / sigma_hat
  Phi <- pnorm(z)
  phi <- dnorm(z)
  
  # Marginal effects
  ME_uncond <- beta_hat * Phi
  ME_cond   <- beta_hat
  ME_prob   <- beta_hat * phi / sigma_hat
  
  # Delta method for standard errors
  vcov_all<- vcov(result)
  get_se  <- function(grad_mat) {
    sqrt(diag(grad_mat %*% vcov_all %*% t(grad_mat)))
  }
  
  k <- length(beta_hat)
  # Gradient matrices
  grad_uncond <- matrix(0, nrow=k, ncol=length(est))
  grad_cond   <- matrix(0, nrow=k, ncol=length(est))
  grad_prob   <- matrix(0, nrow=k, ncol=length(est))
  
  for (j in 1:k) {
    # For unconditional mean effect
    grad_uncond[j, j]   <- Phi
    dPhi_dmu            <- phi / sigma_hat
    dPhi_dsigma         <- -phi * mu_hat / (sigma_hat^2)
    grad_uncond[j, 1:k] <- grad_uncond[j, 1:k] + beta_hat[j] * dPhi_dmu * xbar
    grad_uncond[j, k+1] <- grad_uncond[j, k+1] + beta_hat[j] * dPhi_dsigma
    
    # For conditional mean effect (just beta_j)
    grad_cond[j, j] <- 1
    
    # For probability effect
    grad_prob[j, j]   <- phi / sigma_hat
    dterm_dmu         <- -z * phi / (sigma_hat^2)
    dterm_dsigma      <- -phi / (sigma_hat^2) + z * phi * mu_hat / (sigma_hat^3)
    grad_prob[j, 1:k] <- grad_prob[j, 1:k] + beta_hat[j] * dterm_dmu * xbar
    grad_prob[j, k+1] <- grad_prob[j, k+1] + beta_hat[j] * dterm_dsigma
  }
  
  se_uncond <- get_se(grad_uncond)
  se_cond   <- get_se(grad_cond)
  se_prob   <- get_se(grad_prob)
  
  make_table <- function(est, se) {
    tval             <- est / se
    pval             <- 2 * pnorm(-abs(tval))
    signif_code      <- rep(" ", length(pval))
    signif_code[pval < 0.1]   <- "."
    signif_code[pval < 0.05]  <- "*"
    signif_code[pval < 0.01]  <- "**"
    signif_code[pval < 0.001] <- "***"
    data.frame(
      Estimate     = round(est,5),
      `Std. Error` = round(se,5),
      `t value`    = round(tval,3),
      `Pr(>|t|)`   = format.pval(pval, digits=3, eps=0.001),
      Signif       = signif_code,
      row.names    = names(beta_hat),
      check.names  = FALSE
    )
  }
  
  # --- Print in desired order ---
  cat("\nMarginal Effects on Pr(y>0):\n")
  print(make_table(ME_prob, se_prob))
  
  cat("\nMarginal Effects on E[y|y>0]:\n")
  print(make_table(ME_cond, se_cond))
  
  cat("\nMarginal Effects on E[y]:\n")
  print(make_table(ME_uncond, se_uncond))
}
#*******************************************************************************
# --- Example usage ------------------------------------------------------------
set.seed(123)
n      <- 200
x1     <- rnorm(n)
x2     <- runif(n)
y_star <- 1 + 2*x1 - 1.5*x2 + rnorm(n, sd=0.5)
y      <- ifelse(y_star > 0, y_star, 0)
X      <- cbind(1, x1, x2)
data   <- cbind(y, X)
#*******************************************************************************
# maxLik estimation
logLik_fn <- function(par) gaussian_tobit_ll(par, data)
start     <- c(0,0,0,1)
result    <- maxLik(logLik=logLik_fn, start=start)
names(result$estimate) <- c("Intercept","x1","x2","sigma")

# --- Use custom summary -------------------------------------------------------
summary.tobitMaxLik(result, data)
#*******************************************************************************
#*# --- AER::tobit estimation (left-censoring at 0) ----------------------------
data_df   <- data.frame(y, x1, x2)
tobit_fit <- tobit(y ~ x1 + x2, left = 0, data = data_df)

cat("\n=== AER::tobit Results ===\n")
print(summary(tobit_fit))
#*******************************************************************************
#* Marginal effects for AER package:
library(marginaleffects)
mfx <- avg_slopes(tobit_fit)
cat("\n=== AER:: marginal effects for Tobit Results ===\n")
print(mfx)
#*******************************************************************************