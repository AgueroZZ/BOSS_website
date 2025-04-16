# This function computes the pairwise squared distances for matrices
compute_sq_dist <- function(X, Y) {
  if(is.vector(X)) X <- matrix(X, ncol = 1)
  if(is.vector(Y)) Y <- matrix(Y, ncol = 1)
  XX <- rowSums(X^2)
  YY <- rowSums(Y^2)
  XY <- tcrossprod(X, Y)
  sq_dist <- matrix(XX, ncol=nrow(Y), nrow=nrow(X)) +
    t(matrix(YY, ncol=nrow(X), nrow=nrow(Y))) -
    2 * XY
  return(sq_dist)
}

# Square Exponential covariance function
square_exp_cov_generator_nd <- function(length_scale = 1, signal_var = 1) {
  square_exp_cov <- function(x, x_prime) {

    # Check if x and x_prime are the same, if so only compute half the matrix
    same_input <- identical(x, x_prime)

    # Compute the pairwise squared distances for multivariate data
    sq_dist <- compute_sq_dist(x, x_prime)

    # If inputs are the same, make the matrix symmetric
    if (same_input) {
      upper_tri <- upper.tri(sq_dist)
      sq_dist[upper_tri] <- t(sq_dist)[upper_tri]
    }

    # Compute the covariance
    cov_matrix <- signal_var * exp(-sq_dist / (2 * length_scale^2))

    return(t(cov_matrix))
  }

  return(square_exp_cov)
}



# Function to compute conditional mean and variance of a GP
predict_gp <- function(data, x_pred, noise_var = 1e-6, choice_cov) {
  # Extract x and y from the data
  x_obs <- data$x
  y_obs <- data$y
  if(is.vector(x_obs)) {
    N <- length(x_obs)
    D <- 1
  } else {
    N <- nrow(x_obs)
    D <- ncol(x_obs)
  }

  # Compute covariance matrices
  K_obs_obs <- choice_cov(x_obs, x_obs)
  K_obs_pred <- choice_cov(x_obs, x_pred)
  K_pred_pred <- choice_cov(x_pred, x_pred)

  # Add noise to the diagonal and compute Cholesky decomposition
  K_obs_obs <- K_obs_obs + noise_var * diag(N)
  L <- chol(K_obs_obs)  # Upper triangular factor

  # Solve for conditional mean
  Ly <- forwardsolve(t(L), y_obs)
  cond_mean <- K_obs_pred %*% backsolve(L, Ly)

  # Solve for conditional variance
  LK <- forwardsolve(t(L), t(K_obs_pred))
  cond_var <- K_pred_pred - crossprod(LK)

  # Generate simulation
  sim <- MASS::mvrnorm(1, as.vector(cond_mean), cond_var)

  return(list(x = x_pred, mean = as.vector(cond_mean), var = cond_var, sim = sim))
}

# Compute likelihood using Square Exponential:
compute_like <- function(length_scale, x, y, signal_var, noise_var){
  square_exp_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
  C <- square_exp_cov(x = x, x_prime = x)
  A <- C + noise_var * diag(nrow(x))
  L <- chol(A)  # A = L %*% t(L)

  # Compute the log-determinant of A and thus of Q = A^{-1}
  log_det_A <- sum(log(diag(L)))
  log_det_Q <- -log_det_A

  # Compute Q using the Cholesky inverse (more stable than solve())
  Q <- chol2inv(L)

  # Compute the likelihood in a numerically stable way.
  # Here, tcrossprod(t(y), Q) %*% y is equivalent to t(y) %*% Q %*% y.
  like <- as.numeric((t(y) %*% Q %*% y) / 2 - log_det_Q)

  if(is.nan(like) | is.na(like)){
    return(1e20)
  }
  else if(like == -Inf){
    return(-1e20)
  }
  else if(like == Inf){
    return(1e20)
  }
  else{
    return(like)
  }
}


UCB <- function(x, data, cov, nv, D, d){
  fnew <- predict_gp(data, x, choice_cov = cov, noise_var = nv)

  # Compute the UCB acquisition function
  beta <- 2*log((D^2)*(pi^2)/(6*d))
  return(as.numeric(-fnew$mean - sqrt(beta) * sqrt(fnew$var)))
}


obtain_aghq <- function(f, k = 100, startingvalue = NULL, optresult = NULL){
  if(!is.null(optresult)){
    return(aghq::aghq(ff = ff, k = k, startingvalue = startingvalue, optresults = optresult))
  }
  else{
    ff <- list(fn = f, gr = function(x) numDeriv::grad(f, x), he = function(x) numDeriv::hessian(f, x))
    return(aghq::aghq(ff = ff, k = k, startingvalue = startingvalue))
  }
}


#### Making BO adaptive:
BOSS <- function(func, update_step = 5, max_iter = 100, D = 1,
                 lower = rep(0, D), upper = rep(1, D),
                 noise_var = 1e-6,
                 AGHQ_k = 3, AGHQ_iter_check = 10, AGHQ_eps = 0.1, buffer = 1e-3,
                 initial_design = 5, delta = 0.01, optim.n = 5,
                 opt.lengthscale.grid = NULL,  # Grid-based option for lengthscale optimization.
                 opt.grid = NULL,              # Grid-based option for AF optimization.
                 verbose = 3) {              # Verbosity level: 0, 1, 2, or 3

  # Initialize a helper for verbose printing
  vprint <- function(level, msg) {
    if (verbose >= level) print(msg)
  }

  # If verbose == 1, set up a progress bar.
  if (verbose == 1) {
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  }

  # Check if dimensions of lower and upper bounds match.
  if(length(lower) != D || length(upper) != D) {
    stop("lower and upper must have the same length as the function's input dimension")
  }

  # Initialize matrices/vectors to store evaluations.
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()

  if (verbose == 3) {
    print('Initial fixed evaluation phase...')
  }
  # Initial design: uniformly fixed points.
  initial <- matrix(rep(seq(from = lower, to = upper, length.out = initial_design), D),
                    nrow = initial_design, ncol = D, byrow = FALSE)
  for (i in 1:nrow(initial)) {
    xmat_trans <- rbind(xmat_trans, (initial[i,] - lower)/(upper - lower))
    xmat <- rbind(xmat, initial[i,])
    yvec <- c(yvec, func(initial[i,]))
  }
  # Center the function values.
  rel <- mean(yvec)
  yvec <- yvec - rel
  num_initial <- initial_design
  signal_var <- var(yvec)

  # Initial optimization for lengthscale.
  if (!is.null(opt.lengthscale.grid)) {
    length_scale_vec <- seq(0.01, 0.99, length.out = opt.lengthscale.grid)
    like_vec <- sapply(length_scale_vec, function(l)
      compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   signal_var = signal_var, noise_var = noise_var))
    max_idx <- which.max(like_vec)
    length_scale <- length_scale_vec[max_idx]
    lik <- like_vec[max_idx]
  } else {
    opt <- optim(runif(1, 0.01, 0.99), function(l)
      compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   signal_var = signal_var, noise_var = noise_var),
      control = list(maxit = 100), lower = 0.01, upper = 0.9, method = 'L-BFGS-B')
    length_scale <- opt$par
    lik <- opt$value
  }
  vprint(3, paste("The new length.scale:", length_scale))
  vprint(3, paste("The new signal_var:", signal_var))

  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)

  # Initialize AGHQ-related quantities.
  AGHQ_diff <- Inf
  AGHQ_f_moments_old <- list(first = rep(Inf, D), second = matrix(Inf, nrow = D, ncol = D))
  AGHQ_range_old <- matrix(Inf, nrow = 2, ncol = D)

  i <- 1
  while(i <= max_iter && AGHQ_eps < AGHQ_diff){
    if(verbose == 1) setTxtProgressBar(pb, i)
    if(verbose == 3) {
      print(paste("Iteration:", i))
    } else if(verbose == 2) {
      cat(paste("Iteration:", i, "\n"))
    }

    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

    if(i %% update_step == 0){
      if(verbose == 3) print("Time to update the parameters!")
      signal_var <- var(newdata$y)
      if(!is.null(opt.lengthscale.grid)) {
        length_scale_vec <- seq(0.01, 0.99, length.out = opt.lengthscale.grid)
        like_vec <- sapply(length_scale_vec, function(l)
          -compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                       signal_var = signal_var, noise_var = noise_var))
        max_idx <- which.max(like_vec)
        length_scale <- length_scale_vec[max_idx]
        lik <- like_vec[max_idx]
      } else {
        opt <- optim(runif(1, 0.01, 0.99), function(l)
          compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                       signal_var = signal_var, noise_var = noise_var),
          control = list(maxit = 100), lower = 0.01, upper = 0.9, method = 'L-BFGS-B')
        length_scale <- opt$par
        lik <- opt$value
      }
      if(verbose == 3) {
        print(paste("The new length.scale:", length_scale))
        print(paste("The new signal_var:", signal_var))
      }
      choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
    }

    # Optimize the acquisition function (UCB).
    if(verbose == 3) print("Maximize Acquisition Function")
    if(is.null(opt.grid)){
      # Multi-start local optimization.
      initialize_UCB <- matrix(runif(D*optim.n, rep(buffer, D), rep((1 - buffer), D)),
                               nrow = optim.n, ncol = D, byrow = TRUE)
      optimizer <- optimx::multistart(initialize_UCB, function(x)
        UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
            cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta),
        control = list(maxit = 100),
        lower = rep(buffer, D), upper = rep((1 - buffer), D),
        method = 'L-BFGS-B')
      next_point <- unlist(unname(unique(optimizer[which.min(optimizer$value), 1:D])))
    } else {
      # Grid-based search for the acquisition function.
      grid_list <- replicate(D, seq(buffer, 1 - buffer, length.out = opt.grid), simplify = FALSE)
      grid <- as.matrix(expand.grid(grid_list))
      af_values <- apply(grid, 1, function(x)
        -UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
            cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta))
      best_idx <- which.max(af_values)
      next_point <- grid[best_idx, ]
      # if next_point is already covered, add a small perturbation
      if(any(apply(xmat_trans, 1, function(row) all(abs(row - next_point) == 0)))) {
        next_point <- next_point + runif(D, (-buffer/2), (buffer/2))
      }
    }
    # Update design matrices and evaluate the function.
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point_original <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point_original)
    y_new <- func(next_point_original)
    yvec <- c(yvec, y_new)
    rel <- mean(yvec)
    yvec <- yvec - rel

    if(verbose == 3){
      print(paste("Next point:", next_point_original))
      print(paste("Function value:", y_new))
    } else if(verbose == 2) {
      print(paste("Iteration:", i, "Next point:", next_point_original, "Function value:", y_new))
    }

    # Check convergence with AGHQ every AGHQ_iter_check iterations.
    if(i %% AGHQ_iter_check == 0){
      if(verbose == 3) print("Time to check AGHQ difference!")
      surrogate <- function(xvalue, data_to_smooth) {
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D),
                   choice_cov = choice_cov, noise_var = noise_var)$mean
      }
      lf_design <- list(x = xmat_trans, y = yvec)
      lg_design <- list(x = apply(lf_design$x, 2, qnorm),
                        y = lf_design$y + apply(dnorm(apply(lf_design$x, 2, qnorm), log = TRUE), 1, sum))
      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lg_design))
      AGHQ_f_new <- obtain_aghq(f = fn_new, k = AGHQ_k, startingvalue = rep(0, D))$normalized_posterior$nodesandweights
      AGHQ_f_new$prob <- (AGHQ_f_new$weights * exp(AGHQ_f_new$logpost_normalized))
      AGHQ_f_moments_new <- list()
      AGHQ_f_moments_new$first <- colSums(AGHQ_f_new[,1:D, drop = FALSE] * AGHQ_f_new$prob)
      AGHQ_f_moments_new$second <- t(as.matrix(AGHQ_f_new[,1:D, drop = FALSE])) %*%
        as.matrix(AGHQ_f_new[,1:D, drop = FALSE] * AGHQ_f_new$prob) -
        as.matrix(AGHQ_f_moments_new$first, nrow = D) %*% t(as.matrix(AGHQ_f_moments_new$first, nrow = D))
      difference_in_moments <- list()
      difference_in_moments$first <- sqrt(sum((AGHQ_f_moments_new$first - AGHQ_f_moments_old$first)^2)) /
        sqrt(sum(AGHQ_f_moments_new$first^2))
      difference_in_moments$second <- norm(as.matrix(AGHQ_f_moments_new$second - AGHQ_f_moments_old$second),
                                           type = "F")/norm(as.matrix(AGHQ_f_moments_new$second), type = "F")
      if(AGHQ_k == 1) {
        difference_in_moments$second <- 0
      }
      AGHQ_range_new <- apply(AGHQ_f_new[,1:D, drop = FALSE], 2, range)
      AGHQ_range_diff <- max(abs(AGHQ_range_new - AGHQ_range_old)/abs(AGHQ_range_new))
      AGHQ_diff_new <- max(difference_in_moments$first, difference_in_moments$second, AGHQ_range_diff)
      if(verbose == 3) {
        print(paste("AGHQ rel-difference in first moment:", difference_in_moments$first))
        print(paste("AGHQ rel-difference in second moment:", difference_in_moments$second))
        print(paste("AGHQ rel-difference in range:", AGHQ_range_diff))
      }
      AGHQ_f_moments_old <- AGHQ_f_moments_new
      AGHQ_diff <- AGHQ_diff_new
      AGHQ_range_old <- AGHQ_range_new
    }
    i <- i + 1
  }

  if (AGHQ_diff < AGHQ_eps && AGHQ_eps > 0) {
    if(verbose >= 2) print("Posterior surrogate converged based on AGHQ criteria!")
  } else {
    if(verbose >= 2) {
      print(paste0("Maximum iterations reached! AGHQ difference: ", AGHQ_diff))
    }
  }

  if(verbose == 1) close(pb)

  return(list(result = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var))
}

