suppressPackageStartupMessages(library(eClosure))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rSEA))

# ---- Helper functions for closed e-BH ----
harmonic_number <- function(m) sum(1 / seq_len(m))

calibrator_BY <- function(p, alpha) {
  # Calibrator of Equation 6
  k <- length(p)
  h_k <- harmonic_number(k)
  
  indicator <- (h_k * p <= alpha)
  ceil <- ceiling((k * h_k * p) / alpha)
  
  ifelse(indicator, k / (alpha * ceil), 0)
}

# ---- Decision boundary for double filtering ----
calc_boundary <- function(df, alpha, method = c("closedeBH", "closedBY", "closedSu", "medianTDP", "TDPbound")) {
  start_time <- Sys.time()
  
  n   <- nrow(df)
  ids <- rownames(df)
  
  # Orders the set
  ord_lf <- ids[order(-abs(df$log2FoldChange))]
  ord_p  <- ids[order(df$pvalue)]
  
  # Calculates e-values for the e-BH procedure
  if (method == "closedeBH") {
    df$evalue <- calibrator_BY(df$pvalue, 0.05) # Independent of global alpha
  }
  
  # Tests a candidate discovery set for a given method
  test_set <- switch(method,
                     closedSu  = function(set) closedSu(df$pvalue,  set = set, alpha = alpha),
                     closedBY  = function(set) closedBY(df$pvalue,  set = set, alpha = alpha),
                     closedeBH = function(set) closedeBH(df$evalue, set = set, alpha = alpha),
                     medianTDP = function(set) {
                       if (length(df$gene_id[set]) == 0) return(TRUE)
                       setTDP(df$pvalue, featureIDs = df$gene_id, set = df$gene_id[set], alpha = 0.5)$TDP.bound >= 0.95
                     },
                     TDPbound = function(set) {
                       if (length(df$gene_id[set]) == 0) return(TRUE)
                       setTDP(df$pvalue, featureIDs = df$gene_id, set = df$gene_id[set], alpha = alpha)$TDP.bound >= 0.95
                     })
  
  # Stores steps taken
  k <- 0
  
  # Stores current direction (not strictly necessary)
  dir <- "start"
  
  # Boundary coordinates
  cx <- integer(2 * n)
  cy <- integer(2 * n)
  
  # Stores visited coordinates
  visited <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  
  # Helper functions for above matrix
  is_visited <- function(p, lf) {
    as.logical(visited[p, lf] != 0)
  }
  
  mark_visited <- function(p, lf) {
    visited[p, lf] <<- 1
    return(TRUE)
  }
  
  # Creates stacks for backtracking
  stack_p  <- integer(0)
  stack_lf <- integer(0)
  
  # Adds starting position to above matrix
  lf <- n
  p  <- 1
  
  # Pushes the starting position onto the stack
  stack_p  <- c(stack_p, p)
  stack_lf <- c(stack_lf, lf)
  
  while(TRUE) {
    # If you want to track progress
    # if (k %% 500 == 0) {
      # cat(sprintf(
        # "\n—--- Progress ---—\nMethod: %s\nTime elapsed: %.0f s\nSteps taken: %d\nLargest |lf|: %g\nLargest p: %g\n",
        # method, as.numeric(Sys.time() - start_time, units = "secs"), k, lf, p
      # ))
    # }
    
    # Ends algorithm once boundary is reached
    if (lf < 1 || p < 1 || lf > n || p > n) {
      break
    }
    
    k <- k + 1
    cx[k] <- lf
    cy[k] <- p
    
    try_step <- function(new_lf, new_p, new_dir) {
      if (new_p >= 1 && new_lf >= 1 && new_p <= n && new_lf <= n) { 
        if (is_visited(new_p, new_lf)) {
          return(FALSE)
        }
      }
      
      set   <- rep(FALSE, n)
      inter <- intersect(ord_lf[1:new_lf], ord_p[1:new_p])
      
      if (length(inter)) {
        set[as.integer(inter)] <- TRUE
      }
      
      if (test_set(set)) {
        if (new_p >= 1 && new_lf >= 1 && new_p <= n && new_lf <= n) { 
          mark_visited(new_p, new_lf)
        }
        
        dir <<- new_dir
        lf  <<- new_lf
        p   <<- new_p
        
        stack_p  <<- c(stack_p, p)
        stack_lf <<- c(stack_lf, lf)
        
        return(TRUE)
      }
      
      return(FALSE)
    }
    
    # Tries to take steps in this specific order
    step_taken <- FALSE
    
    if (!step_taken && !(dir == "up"))    step_taken <- try_step(lf,     p + 1, "down")
    if (!step_taken && !(dir == "right")) step_taken <- try_step(lf - 1, p,     "left")
    if (!step_taken && !(dir == "down"))  step_taken <- try_step(lf,     p - 1, "up")
    # Disabled; see footnote underneath Algorithm 1
    # if (!step_taken && !(dir == "left"))  step_taken <- try_step(lf + 1, p,     "right")
    
    # Backtracks if no steps are valid
    if (!step_taken) {
      mark_visited(p, lf)
      
      if (length(stack_p) > 1) {
        # Goes back one step
        stack_p  <- stack_p[-length(stack_p)]
        stack_lf <- stack_lf[-length(stack_lf)]
        
        # Moves to top of the previous stack
        p  <- stack_p[length(stack_p)]
        lf <- stack_lf[length(stack_lf)]
        
        next
      } else {
        break
      }
    }
  }
  
  cx <- stack_lf
  cy <- stack_p
  
  return(list(method = method, steps = k, cx = cx, cy = cy))
}

# ---- Approximation of alpha_min ----
test_result <- function(vals, set, method = c("closedeBH", "closedBY", "closedSu"), alpha) {
  # Corresponds to IsValidDiscoverySet in Algorithm 2
  if (method == "closedeBH") {
    return(closedeBH(vals, set, alpha))
  } else if (method == "closedBY") {
    return(closedBY(vals,  set, alpha))
  } else if (method == "closedSu") {
    return(closedSu(vals,  set, alpha))
  }
}

find_min_alpha <- function(vals, set, method = c("closedeBH", "closedBY", "closedSu")) {
  # Implementation of Algorithm 2
  lb  <- 0
  ub  <- 1
  
  # Checks validity of pathways at bounds
  always_valid   <- test_result(vals, set, method, alpha = lb)
  always_invalid <- test_result(vals, set, method, alpha = ub)
  
  if (always_valid) {
    return(0)
  }  else if (!always_invalid) {
    return(1)
  }
  
  tol      <- 1e-4
  iter     <- 0
  max_iter <- 50
  
  while ((ub - lb) >= tol && iter <= max_iter) {
    mid   <- (lb + ub) / 2
    valid <- test_result(vals, set, method, alpha = mid)
    
    if (valid) {
      ub <- mid
    } else {
      lb <- mid
    }
    
    iter <- iter + 1
  }
  
  return((lb + ub) / 2)
}

# ---- Calculating power of TDP methods ----
TDP_power <- function(df, alpha) {
  # Implementation of Algorithm 3
  ord_by_p <- df$gene_id[order(df$pvalue)]
  
  find_max_k <- function(lo, hi) {
    while (lo < hi) {
      # Biasses upward
      mid <- ceiling((lo + hi + 1) / 2)
      bound <- setTDP(df$pvalue, featureIDs = df$gene_id,
                      set = ord_by_p[1:mid], alpha = alpha)$TDP.bound
      if (bound >= 0.95) {
        lo <- mid
      } else {
        # Can't exceed mid - 1 by definition
        hi <- mid - 1
      }
    }
    lo
  }
  
  find_max_k(0, nrow(df))
}

# ---- Simulation scheme of equicorrelated data ----
simulate_data <- function(n, pi1, rho, signal, df) {
  # Sets non-null indices
  m <- floor(n * pi1)
  non_null <- 1:m
  
  # Generates equicorrelated variables
  b <- sqrt(rho)
  a <- sqrt(1 - rho)
  
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Z1 <- rnorm(1)
  Z2 <- rnorm(1)
  
  # Adds signal to non-null indices
  mu <- numeric(n)
  mu[non_null] <- signal
  
  G1 <- a * X1 + b * Z1
  G2 <- a * X2 + b * Z2 + mu
  
  # Calculates p-values and fold changes
  z <- (G2 - G1) / sqrt(2)
  pvals <- 2 * pnorm(-abs(z))
  
  sigma <- sqrt(df / rchisq(n, df = df))
  fold_change <- z * sigma
  
  data.frame(
    gene_id        = seq_len(n),
    log2FoldChange = fold_change,
    pvalue         = pvals,
    non_null       = (1:n %in% non_null)
  )
}

# ---- Simulation scheme of pathways ----
simulate_pathways <- function(df, m, n, TDP) {
  non_nulls <- which(df[["non_null"]])
  nulls     <- which(!df[["non_null"]])
  
  k   <- floor(n * TDP)
  out <- vector("list", m)
  
  for (i in seq_len(m)) {
    samp_non <- sample(non_nulls, size = k, replace = FALSE)
    samp_nil <- sample(nulls,     size = n - k, replace = FALSE)
    out[[i]] <- c(samp_non, samp_nil)
  }
  out
}
