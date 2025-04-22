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
  # Check if fnew$var is a matrix, if so take out its diagonal
  if(is.matrix(fnew$var)){
    fnew$var <- diag(fnew$var)
  }
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

#### Compute the KL distance:
Compute_KL <- function(x, qx, px){
  to_kept <- which(px > 0)
  x <- x[to_kept]
  qx <- qx[to_kept]
  px <- px[to_kept]
  dx <- diff(x)
  left <- c(0,dx)
  right <- c(dx,0)
  0.5 * sum(left * log(px/qx) * px) + 0.5 * sum(right * log(px/qx) * px)
}

#### Compute the KS distance:
Compute_KS <- function(x, qx, px){
  dx <- c(diff(x),0)
  max(abs(cumsum(qx * dx) - cumsum(px * dx)))
}



#### Making BO adaptive:
#' Single entry point for BOSS
#'
#' @param func The function to optimize.
#' @param criterion One of "modal" (calls BOSS_modal) or "aghq" (calls BOSS_aghq) or "KS" (calls BOSS_KS) or "KL" (calls BOSS_KL).
#' @param ... Further arguments passed on to BOSS_modal or BOSS_aghq.
#' @export
BOSS <- function(func,
                 criterion = c("aghq", "modal", "KS", "KL"),
                 update_step          = 5,
                 max_iter             = 100,
                 D                    = 1,
                 lower                = rep(0, D),
                 upper                = rep(1, D),
                 noise_var            = 1e-6,
                 initial_design       = 5,
                 delta                = 0.01,
                 optim.n              = 5,
                 optim.max.iter       = 1000,
                 opt.lengthscale.grid = NULL,
                 opt.grid             = NULL,
                 # KL-specific
                 KL_iter_check        = 10,
                 KL_check_warmup      = 20,
                 KL_eps               = 0.1,
                 # KS-specific
                 KS_iter_check        = 10,
                 KS_check_warmup      = 20,
                 KS_eps               = 0.1,
                 # modal-specific
                 modal_iter_check     = 10,
                 modal_check_warmup   = 20,
                 modal_k.nn           = 5,
                 modal_eps            = 0.1,
                 # AGHQ-specific
                 AGHQ_k               = 3,
                 AGHQ_iter_check      = 10,
                 AGHQ_check_warmup    = 20,
                 AGHQ_eps             = 0.1,
                 buffer               = 1e-4,
                 verbose              = 3) {

  criterion <- match.arg(criterion)

  # collect the arguments common to both
  args <- list(
    func                 = func,
    update_step          = update_step,
    max_iter             = max_iter,
    D                    = D,
    lower                = lower,
    upper                = upper,
    noise_var            = noise_var,
    initial_design       = initial_design,
    delta                = delta,
    optim.n              = optim.n,
    optim.max.iter       = optim.max.iter,
    opt.lengthscale.grid = opt.lengthscale.grid,
    opt.grid             = opt.grid,
    verbose              = verbose
  )

  if (criterion == "modal") {
    # add the modal‐only args
    args <- c(args, list(
      modal_iter_check   = modal_iter_check,
      modal_check_warmup = modal_check_warmup,
      modal_k.nn         = modal_k.nn,
      modal_eps          = modal_eps
    ))
    return(do.call(BOSS_modal, args))

  }

  else if (criterion == "KL"){
    # add the KL‐only args
    args <- c(args, list(
      KL_iter_check   = KL_iter_check,
      KL_check_warmup = KL_check_warmup,
      KL_eps          = KL_eps
    ))
    return(do.call(BOSS_KL, args))
  }

  else if (criterion == "KS"){
    # add the KS‐only args
    args <- c(args, list(
      KS_iter_check   = KS_iter_check,
      KS_check_warmup = KS_check_warmup,
      KS_eps          = KS_eps
    ))
    return(do.call(BOSS_KS, args))
  }

  else {
    # add the AGHQ‐only args
    args <- c(args, list(
      AGHQ_k            = AGHQ_k,
      AGHQ_iter_check   = AGHQ_iter_check,
      AGHQ_check_warmup = AGHQ_check_warmup,
      AGHQ_eps          = AGHQ_eps,
      buffer            = buffer
    ))
    return(do.call(BOSS_aghq, args))
  }
}


BOSS_modal <- function(func, update_step = 5, max_iter = 100, D = 1,
                       lower = rep(0, D), upper = rep(1, D),
                       noise_var = 1e-6,
                       modal_iter_check = 10,  modal_check_warmup = 20,
                       modal_k.nn = 5, modal_eps = 0.1, # The number of nearest neighbor K should be smaller than modal_iter_check, otherwise there may not be enough changes in the mode
                       initial_design = 5, delta = 0.01,
                       optim.n = 5, optim.max.iter = 1000,
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
    print('Initial evaluation phase...')
  }

  if(D > 1){
    if (verbose == 3) {
      print("Using Latin Hypercube Sampling for initial design for D > 1.")
    }
    initial <- lhs::randomLHS(initial_design, D)
    initial <- t(apply(initial, 1, function(x) x*(upper - lower) + lower))
  }else{
    if (verbose == 3) {
      print("Using equally spaced spaced initial design for D = 1.")
    }
    initial <- matrix(rep(seq(from = (lower), to = (upper), length.out = initial_design), D),
                      nrow = initial_design, ncol = D, byrow = FALSE)
  }

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
      control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
    length_scale <- opt$par
    lik <- opt$value
  }
  vprint(3, paste("The new length.scale:", length_scale))
  vprint(3, paste("The new signal_var:", signal_var))

  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)

  # Initialize modal-related quantities.
  modal_max_diff <- Inf
  modal_f_old <- list(mode = NA, hessian = NA)
  modal_result <- data.frame(modal_max_diff = numeric(0), i = numeric(0))

  # Initialize the grid for the acquisition function, if provided.
  if (!is.null(opt.grid)) {
    grid_list <- replicate(D, seq(0, 1, length.out = opt.grid), simplify = FALSE)
    grid <- as.matrix(expand.grid(grid_list))
  }

  i <- 1
  while(i <= max_iter && modal_eps < modal_max_diff){
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
          control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
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
      initialize_UCB <- matrix(runif(D*optim.n, rep(0, D), rep((1), D)),
                               nrow = optim.n, ncol = D, byrow = TRUE)
      optimizer <- optimx::multistart(initialize_UCB, function(x)
        UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
            cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta),
        control = list(maxit = 100),
        lower = rep(0, D), upper = rep((1), D),
        method = 'L-BFGS-B')
      next_point <- unlist(unname(unique(optimizer[which.min(optimizer$value), 1:D])))
    } else {
      # Grid-based search for the acquisition function.
      af_values <- apply(grid, 1, function(x)
        -UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
             cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta))
      best_idx <- which.max(af_values)
      next_point <- grid[best_idx, ]
      # Take out the best point from the grid.
      grid <- grid[-best_idx, , drop = FALSE]
    }
    # if next_point is already covered, skip
    if(any(apply(xmat_trans, 1, function(row) all(abs(row - next_point) == 0)))) {
      next
      warning("Next point is already covered, skipping.")
    }
    # Update design matrices and evaluate the function.
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point_original <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point_original)
    y_new <- func(next_point_original)

    yvec <- c(yvec + rel, (y_new))
    rel <- mean(yvec)
    yvec <- yvec - rel

    if(verbose == 3){
      print(paste("Next point:", next_point_original))
      print(paste("Function value:", y_new))
    } else if(verbose == 2) {
      print(paste("Iteration:", i, "Next point:", next_point_original, "Function value:", y_new))
    }

    # Check convergence with modal every modal_iter_check iterations.
    if(i %% modal_iter_check == 0  && i >= modal_check_warmup){
      if(verbose == 3) print("Time to check modal difference!")
      surrogate <- function(xvalue, data_to_smooth) {
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D),
                   choice_cov = choice_cov, noise_var = noise_var)$mean
      }

      lf_design <- list(x = xmat_trans, y = yvec)
      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lf_design))
      # find current optimizer
      mode_point <- xmat_trans[which.max(yvec),]

      # find hessian at the current optimizer
      mode_hess <- numDeriv::hessian(fn_new, matrix(mode_point, nrow = 1, ncol = D))
      mode_hess <- (mode_hess + t(mode_hess))/2

      # compute the trace of the hessian
      post_sd <- 1/sqrt(sum(abs(eigen(mode_hess)$values)))

      # Compute Euclidean distances from the best point to all other design points.
      distances <- sqrt(rowSums((xmat_trans - matrix(rep(mode_point, nrow(xmat_trans)),
                                                     nrow = nrow(xmat_trans), byrow = TRUE))^2))
      # Exclude the zero distance (self-distance)
      nonzero_distances <- distances[distances > 0]

      # find average nearest neighbor distance
      nn_dists <- mean(nonzero_distances[order(nonzero_distances)[1:modal_k.nn]])

      #stop if nn_dists is less than certain percentage (acc) of post_sd
      mode <- nn_dists/post_sd

      modal_f_new <- list()
      # opt_surrogate <- optim(par = lf_design$x[which.max(yvec),],
      #                        fn = fn_new,
      #                        method = "L-BFGS-B",
      #                        lower = lower, upper = upper,
      #                        hessian = TRUE)
      modal_f_new$mode <- mode
      modal_f_new$hessian <- 1/post_sd

      modal_diff_new <- list(
        mode = modal_f_new$mode,
        hessian = abs(modal_f_new$hessian - modal_f_old$hessian)/abs(modal_f_old$hessian)
      )

      # if any of the modal differences is NaN or NA, set it to Inf
      modal_diff_new$mode <- ifelse(is.na(modal_diff_new$mode) | is.nan(modal_diff_new$mode), Inf, modal_diff_new$mode)
      modal_diff_new$hessian <- ifelse(is.na(modal_diff_new$hessian) | is.nan(modal_diff_new$hessian), Inf, modal_diff_new$hessian)

      if(verbose == 3) {
        print(paste("Modal rel-difference:", modal_diff_new$mode))
        print(paste("Hessian rel-difference in second moment:", modal_diff_new$hessian))
      }

      modal_max_diff <- max(modal_diff_new$mode, modal_diff_new$hessian)
      modal_diff <- modal_diff_new
      modal_f_old <- modal_f_new
      modal_result <- rbind(modal_result, data.frame(modal_max_diff = modal_max_diff, i = i))
    }
    i <- i + 1
  }

  if (modal_max_diff < modal_eps && modal_eps > 0) {
    if(verbose >= 2) print("Posterior surrogate converged based on modal criteria!")
  } else {
    if(verbose >= 2) {
      print(paste0("Maximum iterations reached! Maximum modal difference: ", modal_max_diff))
    }
  }

  if(verbose == 1) close(pb)

  return(list(result = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var,
              modal_result = modal_result))
}


BOSS_aghq <- function(func, update_step = 5, max_iter = 100, D = 1,
                      lower = rep(0, D), upper = rep(1, D),
                      noise_var = 1e-6,
                      AGHQ_k = 3, AGHQ_iter_check = 10,  AGHQ_check_warmup = 20, AGHQ_eps = 0.1, buffer = 1e-4,
                      initial_design = 5, delta = 0.01,
                      optim.n = 5, optim.max.iter = 1000,
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
    print('Initial evaluation phase...')
  }
  # Initial design: uniformly initial points.
  if(D > 1){
    if (verbose == 3) {
      print("Using Latin Hypercube Sampling for initial design for D > 1.")
    }
    initial <- lhs::randomLHS(initial_design, D)
    initial <- t(apply(initial, 1, function(x) x*(upper - lower - 2*buffer) + lower - buffer))
  }else{
    if (verbose == 3) {
      print("Using equally spaced spaced initial design for D = 1.")
    }
    initial <- matrix(rep(seq(from = (lower + buffer), to = (upper - buffer), length.out = initial_design), D),
                      nrow = initial_design, ncol = D, byrow = FALSE)
  }

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
      -compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   signal_var = signal_var, noise_var = noise_var))
    max_idx <- which.max(like_vec)
    length_scale <- length_scale_vec[max_idx]
    lik <- like_vec[max_idx]
  } else {
    opt <- optim(runif(1, 0.01, 0.99), function(l)
      compute_like(length_scale = l, y = yvec, x = xmat_trans,
                   signal_var = signal_var, noise_var = noise_var),
      control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
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
  AGHQ_result <- data.frame(AGHQ_diff = numeric(0), i = numeric(0))

  # Initialize the grid for the acquisition function, if provided.
  if (!is.null(opt.grid)) {
    grid_list <- replicate(D, seq(buffer, (1 - buffer), length.out = opt.grid), simplify = FALSE)
    grid <- as.matrix(expand.grid(grid_list))
  }

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
          control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
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
      af_values <- apply(grid, 1, function(x)
        -UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
             cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta))
      best_idx <- which.max(af_values)
      next_point <- grid[best_idx, ]
      # Take out the best point from the grid.
      grid <- grid[-best_idx, , drop = FALSE]
    }
    # if next_point is already covered, add a small perturbation
    if(any(apply(xmat_trans, 1, function(row) all(abs(row - next_point) == 0)))) {
      next_point <- next_point + runif(D, (-buffer/2), (buffer/2))
      warning("Next point is already covered, adding a small perturbation.")
    }
    # Update design matrices and evaluate the function.
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point_original <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point_original)
    y_new <- func(next_point_original)

    yvec <- c(yvec + rel, (y_new))
    rel <- mean(yvec)
    yvec <- yvec - rel

    if(verbose == 3){
      print(paste("Next point:", next_point_original))
      print(paste("Function value:", y_new))
    } else if(verbose == 2) {
      print(paste("Iteration:", i, "Next point:", next_point_original, "Function value:", y_new))
    }

    # Check convergence with AGHQ every AGHQ_iter_check iterations.
    if(i %% AGHQ_iter_check == 0  && i >= AGHQ_check_warmup){
      if(verbose == 3) print("Time to check AGHQ difference!")
      surrogate <- function(xvalue, data_to_smooth) {
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D),
                   choice_cov = choice_cov, noise_var = noise_var)$mean
      }
      # surrogate <- function(xvalue, data_to_smooth){
      #   data_to_smooth$y <- data_to_smooth$y - mean(data_to_smooth$y)
      #   predict(ss(x = as.numeric(data_to_smooth$x), y = data_to_smooth$y, df = length(unique(data_to_smooth$y)), m = 2, all.knots = TRUE), x = xvalue)$y
      # }
      lf_design <- list(x = xmat_trans, y = yvec)
      lg_design <- list(x = apply(lf_design$x, 2, qnorm),
                        y = lf_design$y + apply(dnorm(apply(lf_design$x, 2, qnorm), log = TRUE), 1, sum))
      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lg_design))

      AGHQ_f_new <- obtain_aghq(f = fn_new, k = AGHQ_k, startingvalue = lg_design$x[which.max(yvec),])$normalized_posterior$nodesandweights
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
      AGHQ_result <- rbind(AGHQ_result, data.frame(AGHQ_diff = AGHQ_diff, i = i))
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
              length_scale = length_scale, signal_var = signal_var,
              AGHQ_result = AGHQ_result))
}


BOSS_KS <- function(func, update_step = 5, max_iter = 100, D = 1,
                       lower = rep(0, D), upper = rep(1, D),
                       noise_var = 1e-6,
                       KS_iter_check = 10,  KS_check_warmup = 20,
                       KS_eps = 0.1,
                       initial_design = 5, delta = 0.01,
                       optim.n = 5, optim.max.iter = 1000,
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

  # Check if D = 1, if not, return an error saying KS is only for 1D.
  if(D > 1) {
    stop("KS is only implemented for D = 1. Please use criterion = `aghq` or `modal`.")
  }

  # Initialize matrices/vectors to store evaluations.
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()

  if (verbose == 3) {
    print('Initial evaluation phase...')
  }

  if (verbose == 3) {
    print("Using equally spaced spaced initial design for D = 1.")
  }
  initial <- matrix(rep(seq(from = (lower), to = (upper), length.out = initial_design), D),
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
      control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
    length_scale <- opt$par
    lik <- opt$value
  }
  vprint(3, paste("The new length.scale:", length_scale))
  vprint(3, paste("The new signal_var:", signal_var))

  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)

  # Initialize KS-related quantities.
  KS_diff <- Inf
  px <- NULL # current posterior
  qx <- NULL # last posterior
  KS_result <- data.frame(i = NULL, KS = NULL)

  # Initialize the grid for the acquisition function, if provided.
  if (!is.null(opt.grid)) {
    grid_list <- replicate(D, seq(0, 1, length.out = opt.grid), simplify = FALSE)
    grid <- as.matrix(expand.grid(grid_list))
  }

  i <- 1
  while(i <= max_iter && KS_eps < KS_diff){
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
          control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
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
      opt.grid <- 100
      # Multi-start local optimization.
      initialize_UCB <- matrix(runif(D*optim.n, rep(0, D), rep((1), D)),
                               nrow = optim.n, ncol = D, byrow = TRUE)
      optimizer <- optimx::multistart(initialize_UCB, function(x)
        UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
            cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta),
        control = list(maxit = 100),
        lower = rep(0, D), upper = rep((1), D),
        method = 'L-BFGS-B')
      next_point <- unlist(unname(unique(optimizer[which.min(optimizer$value), 1:D])))
    } else {
      # Grid-based search for the acquisition function.
      af_values <- apply(grid, 1, function(x)
        -UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
             cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta))
      best_idx <- which.max(af_values)
      next_point <- grid[best_idx, ]
      # Take out the best point from the grid.
      grid <- grid[-best_idx, , drop = FALSE]
    }
    # if next_point is already covered, skip
    if(any(apply(xmat_trans, 1, function(row) all(abs(row - next_point) == 0)))) {
      next
      warning("Next point is already covered, skipping.")
    }
    # Update design matrices and evaluate the function.
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point_original <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point_original)
    y_new <- func(next_point_original)

    yvec <- c(yvec + rel, (y_new))
    rel <- mean(yvec)
    yvec <- yvec - rel

    if(verbose == 3){
      print(paste("Next point:", next_point_original))
      print(paste("Function value:", y_new))
    } else if(verbose == 2) {
      print(paste("Iteration:", i, "Next point:", next_point_original, "Function value:", y_new))
    }

    # Check convergence with KS every KS_iter_check iterations.
    if(i %% KS_iter_check == 0  && i >= KS_check_warmup){
      if(verbose == 3) print("Time to check KS difference!")
      surrogate <- function(xvalue, data_to_smooth) {
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D),
                   choice_cov = choice_cov, noise_var = noise_var)$mean
      }

      lf_design <- list(x = xmat_trans, y = yvec)
      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lf_design))

      # The grid used to compute the KS distance
      KS_grid <- (seq(
        from = lower,
        to = upper,
        length.out = opt.grid
      ) - lower) / (upper - lower)

      fn_vals <- sapply(KS_grid, function(x) fn_new(x))
      fn_vals <- fn_vals - max(fn_vals)

      # compute px by normalizing the function values
      dx <- diff(KS_grid)
      # Compute the trapezoidal areas and sum them up
      integral_approx <- sum(0.5 * (exp(fn_vals)[-1] + exp(fn_vals)[-length(fn_vals)]) * dx)
      px <- exp(fn_vals) / integral_approx
      if(is.null(qx)){
        qx = rep(Inf, length(px))
      }
      KS_diff <- Compute_KS(x = KS_grid, qx = qx, px = px)
      if(is.na(KS_diff) || is.nan(KS_diff)) {
        KS_diff <- Inf
      }
      if(verbose == 3) {
        print(paste("KS:", KS_diff))
      }
      KS_result <- rbind(KS_result, data.frame(i = i, KS = KS_diff))
      qx <- px
    }
    i <- i + 1
  }

  if (KS_diff < KS_eps && KS_eps > 0) {
    if(verbose >= 2) print("Posterior surrogate converged based on KS criteria!")
  } else {
    if(verbose >= 2) {
      print(paste0("Maximum iterations reached! KS difference: ", KS_diff))
    }
  }

  if(verbose == 1) close(pb)

  return(list(result = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var,
              KS_result = KS_result))
}



BOSS_KL <- function(func, update_step = 5, max_iter = 100, D = 1,
                    lower = rep(0, D), upper = rep(1, D),
                    noise_var = 1e-6,
                    KL_iter_check = 10,  KL_check_warmup = 20,
                    KL_eps = 0.1,
                    initial_design = 5, delta = 0.01,
                    optim.n = 5, optim.max.iter = 1000,
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

  # Check if D = 1, if not, return an error saying KL is only for 1D.
  if(D > 1) {
    stop("KL is only implemented for D = 1. Please use criterion = `aghq` or `modal`.")
  }

  # Initialize matrices/vectors to store evaluations.
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()

  if (verbose == 3) {
    print('Initial evaluation phase...')
  }

  if (verbose == 3) {
    print("Using equally spaced spaced initial design for D = 1.")
  }
  initial <- matrix(rep(seq(from = (lower), to = (upper), length.out = initial_design), D),
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
      control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
    length_scale <- opt$par
    lik <- opt$value
  }
  vprint(3, paste("The new length.scale:", length_scale))
  vprint(3, paste("The new signal_var:", signal_var))

  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)

  # Initialize KL-related quantities.
  KL_diff <- Inf
  px <- NULL # current posterior
  qx <- NULL # last posterior
  KL_result <- data.frame(i = NULL, KL = NULL)

  # Initialize the grid for the acquisition function, if provided.
  if (!is.null(opt.grid)) {
    grid_list <- replicate(D, seq(0, 1, length.out = opt.grid), simplify = FALSE)
    grid <- as.matrix(expand.grid(grid_list))
  }

  i <- 1
  while(i <= max_iter && KL_eps < KL_diff){
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
          control = list(maxit = optim.max.iter), lower = 0.01, upper = 0.99, method = 'L-BFGS-B')
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
      opt.grid <- 100
      # Multi-start local optimization.
      initialize_UCB <- matrix(runif(D*optim.n, rep(0, D), rep((1), D)),
                               nrow = optim.n, ncol = D, byrow = TRUE)
      optimizer <- optimx::multistart(initialize_UCB, function(x)
        UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
            cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta),
        control = list(maxit = 100),
        lower = rep(0, D), upper = rep((1), D),
        method = 'L-BFGS-B')
      next_point <- unlist(unname(unique(optimizer[which.min(optimizer$value), 1:D])))
    } else {
      # Grid-based search for the acquisition function.
      af_values <- apply(grid, 1, function(x)
        -UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
             cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta))
      best_idx <- which.max(af_values)
      next_point <- grid[best_idx, ]
      # Take out the best point from the grid.
      grid <- grid[-best_idx, , drop = FALSE]
    }
    # if next_point is already covered, skip
    if(any(apply(xmat_trans, 1, function(row) all(abs(row - next_point) == 0)))) {
      next
      warning("Next point is already covered, skipping.")
    }
    # Update design matrices and evaluate the function.
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point_original <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point_original)
    y_new <- func(next_point_original)

    yvec <- c(yvec + rel, (y_new))
    rel <- mean(yvec)
    yvec <- yvec - rel

    if(verbose == 3){
      print(paste("Next point:", next_point_original))
      print(paste("Function value:", y_new))
    } else if(verbose == 2) {
      print(paste("Iteration:", i, "Next point:", next_point_original, "Function value:", y_new))
    }

    # Check convergence with KL every KL_iter_check iterations.
    if(i %% KL_iter_check == 0  && i >= KL_check_warmup){
      if(verbose == 3) print("Time to check KL difference!")
      surrogate <- function(xvalue, data_to_smooth) {
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D),
                   choice_cov = choice_cov, noise_var = noise_var)$mean
      }

      lf_design <- list(x = xmat_trans, y = yvec)
      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lf_design))

      # The grid used to compute the KL distance
      KL_grid <- (seq(
        from = lower,
        to = upper,
        length.out = opt.grid
      ) - lower) / (upper - lower)

      fn_vals <- sapply(KL_grid, function(x) fn_new(x))
      fn_vals <- fn_vals - max(fn_vals)

      # compute px by normalizing the function values
      dx <- diff(KL_grid)
      # Compute the trapezoidal areas and sum them up
      integral_approx <- sum(0.5 * (exp(fn_vals)[-1] + exp(fn_vals)[-length(fn_vals)]) * dx)
      px <- exp(fn_vals) / integral_approx
      if(is.null(qx)){
        qx = rep(Inf, length(px))
      }
      KL_diff <- Compute_KL(x = KL_grid, qx = qx, px = px)
      if(is.na(KL_diff) || is.nan(KL_diff)) {
        KL_diff <- Inf
      }
      if(verbose == 3) {
        print(paste("KL:", KL_diff))
      }
      KL_result <- rbind(KL_result, data.frame(i = i, KL = KL_diff))
      qx <- px
    }
    i <- i + 1
  }

  if (KL_diff < KL_eps && KL_eps > 0) {
    if(verbose >= 2) print("Posterior surrogate converged based on KL criteria!")
  } else {
    if(verbose >= 2) {
      print(paste0("Maximum iterations reached! KL difference: ", KL_diff))
    }
  }

  if(verbose == 1) close(pb)

  return(list(result = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var,
              KL_result = KL_result))
}

