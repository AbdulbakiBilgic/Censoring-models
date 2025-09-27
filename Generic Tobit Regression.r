#**************************************************************************************************
# --- Generic Tobit Model Class -------------------------------------------------------------------
tobitModel <- function(formula, data, left = 0, right = Inf) {
  #************************************************************************************************
  #* Generic Tobit Model Class with Analytic/Numeric Gradients for Marginal Effects
  #* Author: Abdulbaki Bilgic
  #************************************************************************************************
  #************************************************************************************************
  #* load required libraries:
  library(maxLik); library(AER); library(censReg); library(numDeriv)
  #************************************************************************************************
  # Model frame:
  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  
  # Censoring indicator
  if (is.finite(left) && is.infinite(right)) {
    # Left censoring:
    censored     <- y <= left
    y_uncensored <- ifelse(censored, left, y)
  } else if (is.finite(right) && is.infinite(left)) {
    # Right censoring:
    censored     <- y >= right
    y_uncensored <- ifelse(censored, right, y)
  } else {
    # Interval censoring (both finite):
    censored     <- y <= left | y >= right
    y_uncensored <- ifelse(y <= left, left, ifelse(y >= right, right, y))
  }
  
  # Prepare data matrix for likelihood function:
  data_matrix    <- cbind(y_uncensored, X)
  
  # Calculate censoring indicators for the output:
  left_censored  <- y <= left
  right_censored <- ifelse(is.finite(right), y >= right, FALSE)
  uncensored     <- !left_censored & !right_censored
  
  # Likelihood function for Gaussian Tobit:
  gaussian_tobit_ll <- function(bs, data) {
    np <- length(bs)
    nx <- np - 1
    b  <- bs[1:nx]
    s  <- exp(bs[np])                                                           # guaranteeing positive sigma
    
    y  <- data[, 1]
    X  <- data[, -1, drop = FALSE]
    
    # Censoring indicators (recalculate inside function):
    left_censored  <- y <= left
    right_censored <- ifelse(is.finite(right), y >= right, FALSE)
    uncensored     <- !left_censored & !right_censored
    
    # Left-censored observations:
    if (any(left_censored)) {
      x_left <- X[left_censored, , drop = FALSE]
      f_left <- pnorm((left - x_left %*% b)/s, log.p = TRUE)
    } else {
      f_left <- numeric(0)
    }
    
    # Right-censored observations:
    if (any(right_censored) && is.finite(right)) {
      x_right <- X[right_censored, , drop = FALSE]
      f_right <- pnorm((x_right %*% b - right)/s, log.p = TRUE)
    } else {
      f_right <- numeric(0)
    }
    
    # Uncensored observations:
    if (any(uncensored)) {
      x_uncens <- X[uncensored, , drop = FALSE]
      y_uncens <- y[uncensored]
      f_uncens <- -log(s) + dnorm((y_uncens - x_uncens %*% b)/s, log = TRUE)
    } else {
      f_uncens <- numeric(0)
    }
    
    return(sum(c(f_left, f_right, f_uncens)))
  }
  
  # Initial values:
  start <- c(rep(0, ncol(X)), log(sd(y_uncensored)))
  
  # Estimation:
  result <- maxLik(logLik = function(par) gaussian_tobit_ll(par, data_matrix), 
                   start = start)
  
  # Name parameters:
  coef_names <- colnames(X)
  names(result$estimate) <- c(coef_names, "log_sigma")
  
  # Return model object:
  structure(list(
    result       = result,
    formula      = formula,
    data         = data_matrix,
    left         = left,
    right        = right,
    y_name       = all.vars(formula)[1],
    x_names      = coef_names,
    n_obs        = nrow(data_matrix),
    n_left       = sum(left_censored),
    n_right      = sum(right_censored),
    n_uncensored = sum(uncensored)
  ), class = "tobitModel")
}
# --- Summary Method for tobitModel --------------------------------------------
summary.tobitModel <- function(object, method = c("analytic", "numeric"), ...) {
  method <- match.arg(method)
  result <- object$result
  data   <- object$data
  left   <- object$left
  right  <- object$right
  
  # Observations summary
  obs_summary <- c(
    "Total"          = object$n_obs,
    "Left-censored"  = object$n_left,
    "Uncensored"     = object$n_uncensored,
    "Right-censored" = object$n_right
  )
  
  cat("Tobit Model Summary\n")
  cat("Formula:", deparse(object$formula), "\n")
  cat("\nObservations:\n")
  print(obs_summary)
  
  # Coefficient table
  est  <- result$estimate
  se   <- sqrt(diag(vcov(result)))
  tval <- est / se
  pval <- 2 * pnorm(-abs(tval))
  
  signif_code      <- rep(" ", length(pval))
  signif_code[pval < 0.1]   <- "."
  signif_code[pval < 0.05]  <- "*"
  signif_code[pval < 0.01]  <- "**"
  signif_code[pval < 0.001] <- "***"
  
  pval_formatted <- format.pval(pval, digits = 3, eps = 0.001)
  
  coef_table <- data.frame(
    Estimate     = round(est, 5),
    `Std. Error` = round(se, 5),
    `t value`    = round(tval, 3),
    `Pr(>|t|)`   = pval_formatted,
    Signif       = signif_code,
    row.names    = names(est),
    check.names  = FALSE
  )
  
  cat("\nCoefficients:\n")
  print(coef_table)
  
  # Wald statistic (excluding sigma):
  k         <- length(object$x_names)
  beta_hat  <- est[1:k]
  vcov_mat  <- vcov(result)[1:k, 1:k, drop = FALSE]
  wald_stat <- t(beta_hat) %*% solve(vcov_mat) %*% beta_hat
  wald_stat <- as.numeric(wald_stat)
  df        <- length(beta_hat)
  p_val     <- 1 - pchisq(wald_stat, df)
  
  cat("\nWald-statistic:", round(wald_stat, 2), "on", df, "Df, p-value:", 
      format.pval(p_val, digits = 3, eps = 0.001), "\n")
  
  # Log-likelihood
  cat("Log-likelihood:", round(logLik(result), 2), "on", object$n_obs, "observations\n")
  
  # Marginal Effects (only for left-censored at 0 case)
  if (left == 0 && is.infinite(right)) {
    cat("\n--- Marginal Effects (at means) ---\n")
    
    X             <- data[, -1, drop = FALSE]
    xbar          <- colMeans(X)
    mu_hat        <- sum(xbar * beta_hat)
    log_sigma_hat <- est["log_sigma"]
    sigma_hat     <- exp(log_sigma_hat)
    
    z         <- mu_hat / sigma_hat
    Phi       <- pnorm(z)
    phi       <- dnorm(z)
    lambda    <- phi / Phi
    
    # Marginal effects:
    ME_uncond <- beta_hat * Phi
    ME_prob   <- beta_hat * (phi / sigma_hat)
    C_factor  <- 1 - lambda * (z + lambda)
    ME_cond   <- beta_hat * C_factor
    
    vcov_all  <- vcov(result)
    k         <- length(beta_hat)
    
    if (method == "analytic") {
      # Analytic gradients:
      grad_uncond  <- matrix(0, nrow = k, ncol = length(est))
      grad_prob    <- matrix(0, nrow = k, ncol = length(est))
      grad_cond    <- matrix(0, nrow = k, ncol = length(est))
      
      dPhi_dmu     <- phi / sigma_hat
      dPhi_dsigma  <- -phi * mu_hat / (sigma_hat^2)
      
      lambda_prime <- -lambda * (z + lambda)
      dC_dz        <- lambda * ((z + lambda) * (z + 2 * lambda) - 1)
      dC_dmu       <- dC_dz * (1 / sigma_hat)
      dC_dsigma    <- dC_dz * (-mu_hat / (sigma_hat^2))
      
      for (j in 1:k) {
        # Unconditional:
        grad_uncond[j, j]     <- Phi
        grad_uncond[j, 1:k]   <- grad_uncond[j, 1:k] + beta_hat[j] * dPhi_dmu * xbar
        grad_uncond[j, k + 1] <- grad_uncond[j, k + 1] + beta_hat[j] * dPhi_dsigma * sigma_hat  # Correction: * sigma_hat
        
        # Probability:
        T0                  <- phi / sigma_hat
        grad_prob[j, j]     <- T0
        dT_dmu              <- -z * phi / (sigma_hat^2)
        dT_dsigma           <- -phi/(sigma_hat^2) + z * phi * mu_hat / (sigma_hat^3)
        grad_prob[j, 1:k]   <- grad_prob[j, 1:k] + beta_hat[j] * dT_dmu * xbar
        grad_prob[j, k + 1] <- grad_prob[j, k + 1] + beta_hat[j] * dT_dsigma * sigma_hat  # Correction: * sigma_hat
        
        # Conditional:
        C <- C_factor
        for (m in 1:k) {
          grad_cond[j, m] <- grad_cond[j, m] + (if (j == m) C else 0) + beta_hat[j] * dC_dmu * xbar[m]
        }
        grad_cond[j, k + 1] <- grad_cond[j, k + 1] + beta_hat[j] * dC_dsigma * sigma_hat  # Correction: * sigma_hat
      }
      
    } else {
      # Numeric gradients:
      grad_fun_uncond <- function(p) {
        b       <- p[1:k]
        sig     <- exp(p[k + 1])
        mu      <- sum(xbar * b)
        ztmp    <- mu / sig
        Phi_tmp <- pnorm(ztmp)
        b * Phi_tmp
      }
      
      grad_fun_cond <- function(p) {
        b          <- p[1:k]
        sig        <- exp(p[k + 1])
        mu         <- sum(xbar * b)
        ztmp       <- mu / sig
        phi_tmp    <- dnorm(ztmp)
        Phi_tmp    <- pnorm(ztmp)
        lambda_tmp <- phi_tmp / Phi_tmp
        Ctmp       <- 1 - lambda_tmp * (ztmp + lambda_tmp)
        b * Ctmp
      }
      
      grad_fun_prob <- function(p) {
        b       <- p[1:k]
        sig     <- exp(p[k + 1])
        mu      <- sum(xbar * b)
        ztmp    <- mu / sig
        phi_tmp <- dnorm(ztmp)
        b * (phi_tmp / sig)
      }
      
      grad_uncond <- jacobian(grad_fun_uncond, est)
      grad_cond   <- jacobian(grad_fun_cond, est)
      grad_prob   <- jacobian(grad_fun_prob, est)
    }
    
    get_se <- function(grad_mat) sqrt(diag(grad_mat %*% vcov_all %*% t(grad_mat)))
    se_uncond <- get_se(grad_uncond)
    se_cond   <- get_se(grad_cond)
    se_prob   <- get_se(grad_prob)
    
    make_table <- function(est, se, title) {
      tval <- est / se
      pval <- 2 * pnorm(-abs(tval))
      signif_code               <- rep(" ", length(pval))
      signif_code[pval < 0.1]   <- "."
      signif_code[pval < 0.05]  <- "*"
      signif_code[pval < 0.01]  <- "**"
      signif_code[pval < 0.001] <- "***"
      
      df <- data.frame(
        Estimate     = round(est, 5),
        `Std. Error` = round(se, 5),
        `t value`    = round(tval, 3),
        `Pr(>|t|)`   = format.pval(pval, digits = 3, eps = 0.001),
        Signif       = signif_code,
        row.names    = object$x_names,
        check.names  = FALSE
      )
      
      cat("\n", title, ":\n", sep = "")
      print(df)
    }
    
    make_table(ME_prob, se_prob, "Marginal Effects on Pr(y>0)")
    make_table(ME_cond, se_cond, "Marginal Effects on E[y|y>0]")
    make_table(ME_uncond, se_uncond, "Marginal Effects on E[y]")
    
    cat("\nMethod used for standard errors:", method, "\n")
  }
  
  invisible(object)
}

# --- Print Method -------------------------------------------------------------
print.tobitModel <- function(x, ...) {
  cat("Tobit Model\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Left-censored:", x$n_left, "\n")
  cat("Uncensored:", x$n_uncensored, "\n")
  cat("Right-censored:", x$n_right, "\n")
  cat("Log-likelihood:", round(logLik(x$result), 2), "\n")
}

# --- Predict Method -----------------------------------------------------------
predict.tobitModel <- function(object, 
                               newdata = NULL, 
                               type    = c("response", "probability", "conditional")) {
  type <- match.arg(type)
  
  if (is.null(newdata)) {
    X  <- object$data[, -1, drop = FALSE]
  } else {
    mf <- model.frame(object$formula, newdata)
    X  <- model.matrix(object$formula, mf)
  }
  
  coefs <- coef(object$result)[1:ncol(X)]
  sigma <- exp(coef(object$result)["log_sigma"])
  
  linear_pred <- X %*% coefs
  z           <- linear_pred / sigma
  
  switch(type,
         response = { # E[y]
           pnorm(z) * (linear_pred + sigma * dnorm(z)/pnorm(z))
         },
         probability = { # P(y > 0)
           pnorm(z)
         },
         conditional = { # E[y|y>0]
           linear_pred + sigma * (dnorm(z)/pnorm(z))
         })
}
#**************************************************************************************************