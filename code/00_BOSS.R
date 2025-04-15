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
compute_like <- function(length_scale, x, y, signal_var, noise_var, D, prior_l_mean, prior_l_sd){
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
  like <- as.numeric((t(y) %*% Q %*% y) / 2 - log_det_Q -
                       dlnorm(length_scale, prior_l_mean + log(D)/2, prior_l_sd, log = TRUE))

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
BO_adap_optim_AGHQ <- function(func, update_step = 5, max_iter = 100, D = 1,
                               lower = rep(0, D), upper = rep(1, D),
                               noise_var = 1e-6, prior_l_mean = log(1), prior_l_sd = 1,
                               AGHQ_iter_check = 10, AGHQ_eps = 0.1, buffer = 1e-3, # to avoid aghq-transformation to infinity
                               AGHQ_k = 3,
                               initial_design = 5, delta = 0.01){

  # Check if dimensions of lower_bounds and upper_bounds match
  if(length(lower) != D | length(upper) != D) {
    stop("lower_bounds and upper_bounds must have the same length as function input dimension")
  }

  # Initialize xvec and yvec to store evaluations
  xmat <- c()
  xmat_trans <- c()
  yvec <- c()

  print('Initial random evaluation phase...')
  initial <- matrix(runif(initial_design * D, lower, upper),
                    nrow = initial_design, ncol = D, byrow = T)

  for (i in 1:nrow(initial)) {
    # Add the initial point to xvec and evaluate the function
    xmat_trans <- rbind(xmat_trans, (initial[i,] - lower)/(upper - lower))
    xmat <- rbind(xmat, initial[i,])
    yvec <- c(yvec, func(initial[i,]))
  }
  # Assign the reference value
  rel <- mean(yvec)
  yvec <- yvec - rel

  num_initial <- initial_design

  signal_var = var(yvec)
  opt <- optim(runif(1, 0.01, 0.9), function(l) compute_like(length_scale = l, y = yvec, x = xmat_trans,
                                                             signal_var = signal_var, noise_var = noise_var,
                                                             D, prior_l_mean, prior_l_sd),
               control = list(maxit = 100), lower = 0.01, upper = 0.9, method = 'L-BFGS-B')
  length_scale <- opt$par
  lik <- opt$value

  print(paste("The new length.scale:", length_scale))
  print(paste("The new signal_var:", signal_var))

  choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)

  # Perform Bayesian Optimization

  # Initialize AGHQ here
  # Compute initial AGHQ points
  # initialize AGHQ_f_moments_old, AGHQ_range_old
  AGHQ_diff <- Inf
  AGHQ_f_moments_old <- list(first = rep(Inf, D), second = matrix(Inf, nrow = D, ncol = D))
  AGHQ_range_old <- matrix(Inf, nrow = 2, ncol = D)
  i <- 1
  while(i <= max_iter & AGHQ_eps < AGHQ_diff){
    print(paste("Iteration:", i))
    newdata <- list(x = xmat_trans, x_original = xmat, y = yvec)

    if(i %% update_step == 0){
      print(paste("Time to update the parameters!"))
      signal_var = var(newdata$y)
      opt <- optim(runif(1, 0.01, 0.9), function(l) compute_like(length_scale = l, y = newdata$y, x = newdata$x,
                                                                 signal_var = signal_var, noise_var = noise_var,
                                                                 D, prior_l_mean, prior_l_sd),
                   control = list(maxit = 100), lower = 0.01, upper = 0.9, method = 'L-BFGS-B')
      length_scale <- opt$par
      lik <- opt$value

      print(paste("The new length.scale:", length_scale))
      print(paste("The new signal_var:", signal_var))
      choice_cov <- square_exp_cov_generator_nd(length_scale = length_scale, signal_var = signal_var)
    }

    # Update the GP
    print('Maximize Acquisition Function')
    initialize_UCB <- matrix(runif(D, lower, upper),
                             nrow = 1, ncol = D, byrow = T)
    initialize_UCB <- t((t(initialize_UCB) - lower)/(upper - lower))
    next_point <- optim(initialize_UCB, function(x) UCB(x = matrix(x, nrow = 1, ncol = D), data = newdata,
                                                        cov = choice_cov, nv = noise_var, D = i + num_initial, d = delta),
                        control = list(maxit = 100), lower = rep(buffer, D),
                        upper = rep((1 - buffer), D), method = 'L-BFGS-B')$par

    # Select the next point to evaluate
    xmat_trans <- rbind(xmat_trans, next_point)
    next_point <- next_point*(upper - lower) + lower
    xmat <- rbind(xmat, next_point)
    yvec <- c(yvec + rel, (func(next_point)))
    rel <- mean(yvec)
    yvec <- yvec - rel

    # Print some information for debugging
    print(paste("Next point:", next_point))
    print(paste("Function value:", yvec[i + 1]))

    if(i %% AGHQ_iter_check == 0){
      print(paste("Time to check AGHQ difference!"))

      # surrogate function to be used in AGHQ
      surrogate <- function(xvalue, data_to_smooth){
        predict_gp(data = data_to_smooth, x_pred = matrix(xvalue, ncol = D), choice_cov = choice_cov, noise_var = noise_var)$mean
      }

      # The set of design points
      lf_design <- list(x = xmat_trans, y = yvec)
      # The set of design points after quantile transformation
      lg_design <- list(x = apply(lf_design$x, 2, qnorm),
                        y = lf_design$y + apply(dnorm(apply(lf_design$x, 2, qnorm), log = TRUE), 1, sum)
      )

      fn_new <- function(y) as.numeric(surrogate(xvalue = y, data_to_smooth = lg_design))

      # Update the AGHQ difference here
      AGHQ_f_new <- (obtain_aghq(f = fn_new, k = AGHQ_k, startingvalue = rep(0,D))$normalized_posterior$nodesandweights)

      # AGHQ_diff update
      AGHQ_f_new$prob <- (AGHQ_f_new$weights * exp(AGHQ_f_new$logpost_normalized))

      # Compare the two AGHQ difference in terms of the moments
      AGHQ_f_moments_new <- list()
      AGHQ_f_moments_new$first <- colSums(AGHQ_f_new[,1:D, drop = FALSE] * AGHQ_f_new$prob)
      AGHQ_f_moments_new$second <- t(as.matrix(AGHQ_f_new[,1:D, drop = FALSE])) %*% as.matrix(AGHQ_f_new[,1:D, drop = FALSE] * AGHQ_f_new$prob) - as.matrix(AGHQ_f_moments_new$first, nrow = D) %*% t(as.matrix(AGHQ_f_moments_new$first, nrow = D))

      difference_in_moments <- list()
      difference_in_moments$first <- sqrt(sum((AGHQ_f_moments_new$first - AGHQ_f_moments_old$first)^2))/sqrt(sum(AGHQ_f_moments_new$first^2))
      difference_in_moments$second <- norm(as.matrix(AGHQ_f_moments_new$second - AGHQ_f_moments_old$second), type = "F")/norm(as.matrix(AGHQ_f_moments_new$second), type = "F")

      # When k = 1, the second moment is always 0
      if(AGHQ_k == 1){
        difference_in_moments$second <- 0
      }

      # Compare the two AGHQ difference in ranges
      AGHQ_range_new <- apply(AGHQ_f_new[,1:D, drop = FALSE], 2, range)
      AGHQ_range_diff <- max(abs(AGHQ_range_new - AGHQ_range_old)/abs(AGHQ_range_new))
      AGHQ_diff_new <- max(difference_in_moments$first, difference_in_moments$second, AGHQ_range_diff)

      # print the AGHQ difference in each moment and the AGHQ range
      print(paste("AGHQ rel-difference in first moment:", difference_in_moments$first))
      print(paste("AGHQ rel-difference in second moment:", difference_in_moments$second))
      print(paste("AGHQ rel-difference in range:", AGHQ_range_diff))

      # update old AGHQ_f_moments, AGHQ_diff, AGHQ_range
      AGHQ_f_moments_old <- AGHQ_f_moments_new
      AGHQ_diff <- AGHQ_diff_new
      AGHQ_range_old <- AGHQ_range_new
    }
    i <- i + 1
  }

  if(AGHQ_diff < AGHQ_eps){
    print("Posterior surrogate converged based on AGHQ criteria!")
  }
  else{
    print("Maximum iterations reached and posterior surrogate did not converge based on AGHQ criteria. Adjust your expectation on final results!")
  }

  # Return the result
  return(list(result = list(x = xmat_trans, x_original = xmat, y = yvec + rel),
              length_scale = length_scale, signal_var = signal_var))
}
