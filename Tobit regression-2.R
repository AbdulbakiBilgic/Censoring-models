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
  obs_summary<- c("Total" = n_total,
                  "Left-censored" = n_left,
                  "Uncensored" = n_uncens,
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
  
  # Format p-values with 3 decimal places
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

# --- Use custom summary ---
summary.tobitMaxLik(result, data)
#*******************************************************************************
#*# --- AER::tobit estimation (left-censoring at 0) ---
data_df   <- data.frame(y, x1, x2)
tobit_fit <- tobit(y ~ x1 + x2, left = 0, data = data_df)

cat("\n=== AER::tobit Results ===\n")
print(summary(tobit_fit))
#*******************************************************************************