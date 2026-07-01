if (!require(tidyverse)) {
  install.packages("tidyverse")
}
if (!require(pracma)) {
  install.packages("pracma")
}
if (!require(lubridate)) {
  install.packages("lubridate")
}
if (!require(fixest)) {
  install.packages("fixest")
}

### Helper Function: time_shift
time_shift <- function(y0, i_index, t_index, L) {
  if (is.null(dim(y0))) {
    y0 <- matrix(y0, nrow = length(y0), ncol = 1)  # coerces vectors into column matrices
  }

  y1 <- matrix(NaN, nrow = nrow(y0), ncol = ncol(y0))
  for (i in unique(i_index)) {
    # extract unit data
    i_unit  <- i_index == i
    t_unit  <- t_index[i_unit, drop = FALSE]
    y0_unit <- y0[i_unit, , drop = FALSE]

    # obtain data lag (negative L) or lead (positive L)
    is_t    <- is.element(t_unit + L, t_unit)
    loc_t   <- match(t_unit + L, t_unit, nomatch = 0)
    y1_unit <- matrix(NaN, nrow = nrow(y0_unit), ncol = ncol(y0_unit))
    y1_unit[is_t, ] <- y0_unit[loc_t[is_t], , drop = FALSE]
    y1[i_unit, ]    <- y1_unit
  }

  return(y1)
}

### Helper Function: regress_HDFE
### Helper Function: warn_if_illconditioned
### Emits a STRONG warning when the (residualized) design matrix is severely
### ill-conditioned -- i.e. a regressor is (near-)collinear with the others, so the
### specification is effectively unidentified. At high condition numbers the OLS
### coefficient estimates are numerically unreliable. kappa here is the 2-norm
### condition number (sigma_max/sigma_min of X_resid). This only DETECTS & REPORTS;
### it does not alter the numerics (no-op on well-conditioned specifications).
### Tunable via kappa_max.
warn_if_illconditioned <- function(X_resid, kappa_max = 1e5) {
  X_resid <- as.matrix(X_resid)
  if (ncol(X_resid) == 0 || nrow(X_resid) == 0) return(invisible(NA_real_))
  sv      <- svd(X_resid, nu = 0, nv = 0)$d
  kappa_X <- if (length(sv) == 0 || min(sv) <= 0) Inf else sv[1] / min(sv)
  if (kappa_X > kappa_max) {
    warning(sprintf(paste0(
      "SINGULARITY SAFEGUARD: the design matrix is severely ill-conditioned ",
      "(condition number = %.3g > %.0g). One or more regressors are nearly collinear ",
      "with the rest, so the specification is effectively unidentified: the OLS ",
      "coefficient estimates are numerically unreliable. Check the specification for ",
      "a redundant or near-duplicate control / interaction term."),
      kappa_X, kappa_max), call. = FALSE)
  }
  invisible(kappa_X)
}

regress_HDFE <- function(y, X, FE, tol = 1e-8, max_iter = 10000) {
  if (ncol(FE) == 0) {
    warn_if_illconditioned(X)
    b <- pracma::mldivide(X, y)
    return(list(
      b = unname(as.matrix(b)),
      y_resid = as.matrix(y),
      X_resid = as.matrix(X)
    ))
  } else {
    fe_list <- list()
    for (i in seq_len(ncol(FE))) {
      fe_list[[paste0("fe", i)]] <- FE[, i]
    }

    # demean y and X column by column
    y_resid <- fixest::demean(y, fe_list, iter = max_iter, tol = tol, na.rm = FALSE)
    y_resid[!is.finite(y_resid)] <- 0

    X_resid <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    for (j in seq_len(ncol(X))) {
      X_resid[, j] <- fixest::demean(X[, j], fe_list, iter = max_iter, tol = tol, na.rm = FALSE)
    }
    X_resid[!is.finite(X_resid)] <- 0

    warn_if_illconditioned(X_resid)
    b <- pracma::mldivide(X_resid, y_resid)

    return(list(
      b = unname(as.matrix(b)),
      y_resid = unname(as.matrix(y_resid)),
      X_resid = unname(as.matrix(X_resid))
    ))
  }
}

### Helper Function: pinv_m
### Pseudo-inverse with an EXPLICIT rank tolerance. pracma::pinv defaults to
### tol = .Machine$double.eps^(2/3) (~3.67e-11), applied relatively as s$d > tol*s$d[1].
### On the near-singular Gram matrix X_t'X_t in the Imbens-Kolesar df block
### (kappa(X_t'X_t) ~ 1e11) that loose threshold drops the smallest singular directions,
### which makes the degrees-of-freedom estimate sensitive to the truncation rank. Pinning
### the relative tolerance to max(dim)*eps (~2.2e-15) retains the full numerical rank and
### yields a stable df.
pinv_m <- function(A) pracma::pinv(A, tol = max(dim(A)) * .Machine$double.eps)

### Helper Function: warn_if_illconditioned_ik
### IK small-sample (X_t) ill-conditioning detector. WARNs that a specification is
### ill-conditioned rather than silently rounding it off with a numerical floor. It
### fires on EITHER of two triggers, per (h, s-component):
###   (1) a within-period instrument sum-of-squares s's below ss_tol -- i.e. where the
###       per-period projection denominator is ~0 (the projection is ~0/0);
###   (2) the 2-norm condition number of the IK design X_t exceeding kappa_max.
### DETECT-ONLY: never alters the numerics. Both tolerances are tunable.
warn_if_illconditioned_ik <- function(ss_t, X_t, h, i_s, ss_tol = 1e-12, kappa_max = 1e8) {
  ss_min <- suppressWarnings(min(ss_t))
  if (!is.finite(ss_min) || ss_min < ss_tol) {
    warning(sprintf(paste0(
      "ILL-CONDITIONED IK SPEC (h=%d, s-component %d): a within-period instrument ",
      "sum-of-squares s's = %.3g fell below tolerance %.0g. The per-period projection ",
      "denominator is ~0, so the small-sample (X_t) refinement is unstable and the df/SE ",
      "for this cell are unreliable. Consider rescaling the instrument or dropping the ",
      "degenerate period."),
      h, i_s, ss_min, ss_tol), call. = FALSE)
  }
  if (all(is.finite(X_t))) {
    sv      <- svd(X_t, nu = 0, nv = 0)$d
    kappa_X <- if (length(sv) == 0 || min(sv) <= 0) Inf else sv[1] / min(sv)
  } else {
    kappa_X <- Inf
  }
  if (kappa_X > kappa_max) {
    warning(sprintf(paste0(
      "ILL-CONDITIONED IK SPEC (h=%d, s-component %d): the small-sample design matrix X_t ",
      "is severely ill-conditioned (condition number = %.3g > %.0g). The Imbens-Kolesar ",
      "degrees-of-freedom estimate is numerically unstable; treat the small-sample ",
      "inference for this cell with caution."),
      h, i_s, kappa_X, kappa_max), call. = FALSE)
  }
  invisible(c(ss_min = ss_min, kappa_X = kappa_X))
}

### Main Function
panel_LP <- function(y, X, s = NULL, i_index, t_index,
                     W = NULL, FE = NULL, H = NULL, p_max = NULL, small_sample = FALSE, cumulative = FALSE) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  if (!is.null(s)) {
    s <- as.matrix(s)
  } else {
    s <- matrix(1, nrow = nrow(y), ncol = 1)
  }

  i_index <- as.matrix(i_index)
  t_index <- as.matrix(t_index)

  if (any(nrow(y) != nrow(X), nrow(y) != nrow(t_index), nrow(y) != nrow(i_index))) {
    stop("Row counts of inputs unequal.")
  }

  # normalize time indexes if needed
  if (all(is.Date(t_index))) {
    t_index <- 12 * (year(t_index) - min(year(t_index))) + month(t_index)
  }

  t_min   <- min(t_index)
  t_diff  <- min(diff(sort(unique(t_index), decreasing = FALSE)))
  t_index <- as.matrix((t_index - t_min) / t_diff + 1)

  # drop units with inconsistent time indexes
  keep    <- !(t_index - floor(t_index) > 1e-2)
  y       <- y[keep]
  s       <- s[keep, , drop = FALSE]
  X       <- X[keep, , drop = FALSE]
  i_index <- i_index[keep]
  t_index <- round(t_index[keep])

  # recover dimensions
  n_obs <- length(t_index)
  T_eff <- length(unique(t_index))
  n_s   <- ncol(s) * ncol(X)

  # construct interacted regressor: columns are s[,j_s] * X[,j_x], s-index varies faster
  sX <- matrix(0, n_obs, n_s)
  for (j_x in seq_len(ncol(X))) {
    for (j_s in seq_len(ncol(s))) {
      sX[, (j_x - 1) * ncol(s) + j_s] <- s[, j_s] * X[, j_x]
    }
  }

  # recover optional data
  if (!is.null(W)) {
    if (!is.numeric(W)) {
      stop("Controls W must be numeric.")
    } else {
      W <- as.matrix(W)
      W <- W[keep, , drop = FALSE]
    }
  } else {
    W <- matrix(0, nrow = n_obs, ncol = 0)
  }

  if (!is.null(FE)) {
    FE <- as.matrix(FE)
    FE <- FE[keep, , drop = FALSE]
  } else {
    FE <- matrix(0, nrow = n_obs, ncol = 0)
  }

  # recover optional arguments
  if (is.null(H)) {
    H <- ceiling(0.25 * T_eff)
  }
  if (is.null(p_max)) {
    p_max <- ceiling((T_eff - H)^(1/3))
  }

  # preallocate output
  LP_estimate <- matrix(0, H+1, n_s)
  LP_SE       <- matrix(0, H+1, n_s)
  LP_df       <- matrix(0, H+1, n_s)
  LP_CI90     <- array(0, c(H+1, n_s, 2))
  LP_CI95     <- array(0, c(H+1, n_s, 2))
  LP_CI99     <- array(0, c(H+1, n_s, 2))
  LP_pval     <- matrix(0, H+1, n_s)

  # iterate over horizons
  y_h <- matrix(0, n_obs, 1)
  for (h in 0:H) {
    # lead regressand
    if (cumulative) {
      y_h <- y_h + time_shift(y, i_index, t_index, h)
    } else {
      y_h <- time_shift(y, i_index, t_index, h)
    }

    # construct lagged controls
    p     <- min(h, p_max)
    W_lag <- array(NaN, c(n_obs, p, 1+n_s))
    if (p != 0) {         # construction of lagged controls only takes place if p > 0
      for (j in seq_len(p)) {
        W_lag[, j, 1] <- time_shift(y, i_index, t_index, -j)
        for (i_s in seq_len(n_s)) {
          W_lag[, j, 1 + i_s] <- time_shift(sX[, i_s, drop = FALSE], i_index, t_index, -j)
        }
      }
    }
    W_lag <- array(W_lag, c(length(i_index), (1+n_s)*p))

    # prepare data
    d      <- !apply(is.na(cbind(y_h, sX, W_lag, W, FE)), 1, any)
    y_LP   <- y_h[d, , drop = FALSE]
    X_LP   <- cbind(sX[d, , drop = FALSE], W_lag[d, , drop = FALSE], W[d, , drop = FALSE])
    n_X    <- ncol(X_LP)
    dum_LP <- FE[d, , drop = FALSE]

    # prepare time-series indexes
    t_LP  <- t_index[d]
    t_set <- sort(unique(t_LP))
    T     <- length(t_set)

    # compute LP estimator
    HDFE_output <- regress_HDFE(y_LP, X_LP, dum_LP)
    b_LP <- HDFE_output$b
    y_LP <- HDFE_output$y_resid
    X_LP <- HDFE_output$X_resid
    LP_estimate[h+1, ] <- b_LP[1:n_s]

    # compute score and hessian
    Xv_it <- sweep(X_LP, 1, (y_LP - X_LP %*% b_LP), FUN = "*")
    Xv_t  <- matrix(0, T, n_X)
    for (t in seq_len(T)) {
      t_tmp     <- t_LP == t_set[t]
      Xv_t[t, ] <- colSums(Xv_it[t_tmp, , drop = FALSE])
    }
    XX <- t(X_LP) %*% X_LP

    # compute t-LAHR standard error
    if (small_sample) {
      for (i_s in seq_len(n_s)) {
        # Compute Imbens-Kolesar small-sample refinement
        s_LP <- s[d, i_s, drop = FALSE]
        X_t  <- matrix(0, T, n_X)
        ss_t <- numeric(T)                       # per-period s'~s denominator
        for (t in seq_len(T)) {
          t_tmp    <- t_LP == t_set[t]
          ss_t[t]  <- crossprod(s_LP[t_tmp, 1])
          X_t[t, ] <- (s_LP[t_tmp, 1] %*% X_LP[t_tmp, ]) / ss_t[t]   # ill-conditioning checked below
        }
        warn_if_illconditioned_ik(ss_t, X_t, h, i_s)
        # A truly-degenerate period (s'~s = 0) makes X_t non-finite, on which
        # pracma::pinv -> svd() errors. Set df/SE to NaN for this cell and skip so a
        # single degenerate period does not abort the run.
        if (!all(is.finite(X_t))) {
          LP_SE[h+1, i_s] <- NaN
          LP_df[h+1, i_s] <- NaN
          next
        }
        P0     <- diag(x = 1, T) - X_t %*% pinv_m(X_t)            # explicit rank tolerance
        Xv_var <- t(Xv_t / sqrt(diag(P0))) %*% (Xv_t / sqrt(diag(P0)))
        b_var  <- pracma::pinv(XX) %*% Xv_var %*% pracma::pinv(XX)
        G0     <- matrix(0, T, T)
        XX0    <- pinv_m(t(X_t) %*% X_t)                          # explicit rank tolerance
        for (t in seq_len(T)) {
          G0[, t] <- P0[, t, drop = FALSE] %*% X_t[t, , drop = FALSE] %*% XX0[, 1, drop = FALSE] / sqrt(P0[t, t])
        }
        lam0 <- eigen(t(G0) %*% G0)$values

        # store standard error and degrees of freedom
        LP_SE[h+1, i_s] <- sqrt(pmax(0, diag(b_var)[i_s]))
        LP_df[h+1, i_s] <- sum(lam0)^2 / sum(lam0^2)
      }
    } else {
      # compute time-clustered sandwich formula
      Xv_var <- t(Xv_t) %*% Xv_t
      b_var  <- pracma::pinv(XX) %*% Xv_var %*% pracma::pinv(XX)
      # store standard error and degrees of freedom
      # drop = FALSE: keep the n_s x n_s block a matrix so diag() takes its diagonal
      # even when n_s == 1 (a bare scalar makes R's diag() build an identity matrix).
      LP_SE[h+1, ] <- sqrt(pmax(0, diag(b_var[1:n_s, 1:n_s, drop = FALSE])))
      LP_df[h+1, ] <- Inf
    }

    # compute confidence intervals (LP_df[h+1, ] is a length-n_s vector;
    # qt/pt vectorize over df so each s-component uses its own degrees of freedom)
    cv_tmp            <- qt(1 - (1-0.90)/2, LP_df[h+1, ])
    LP_CI90[h+1, , 1] <- LP_estimate[h+1, ] - cv_tmp * LP_SE[h+1, ]
    LP_CI90[h+1, , 2] <- LP_estimate[h+1, ] + cv_tmp * LP_SE[h+1, ]
    cv_tmp            <- qt(1 - (1-0.95)/2, LP_df[h+1, ])
    LP_CI95[h+1, , 1] <- LP_estimate[h+1, ] - cv_tmp * LP_SE[h+1, ]
    LP_CI95[h+1, , 2] <- LP_estimate[h+1, ] + cv_tmp * LP_SE[h+1, ]
    cv_tmp            <- qt(1 - (1-0.99)/2, LP_df[h+1, ])
    LP_CI99[h+1, , 1] <- LP_estimate[h+1, ] - cv_tmp * LP_SE[h+1, ]
    LP_CI99[h+1, , 2] <- LP_estimate[h+1, ] + cv_tmp * LP_SE[h+1, ]

    # compute p-values of two-sided significance test
    LP_pval[h+1, ] <- 2*(1 - pt(abs(LP_estimate[h+1, ] / LP_SE[h+1, ]), LP_df[h+1, ]))
  }

  # store output
  return(list(estimate = LP_estimate,
              SE = LP_SE,
              df = LP_df,
              CI90 = LP_CI90,
              CI95 = LP_CI95,
              CI99 = LP_CI99,
              pval = LP_pval))
}
