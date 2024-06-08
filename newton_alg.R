library(data.table)

# NEWTON-RAPHSON ALGORITHM

first_deritative_theta_k <- function(T_k, delta_k, theta_k, lambda_k){
  d_k <- sum(delta_k)
  
  d_k * 1 / theta_k - d_k * log(lambda_k) + 
    sum(delta_k * log(T_k)) + 
    log(lambda_k) * lambda_k ^ (-theta_k) * sum(T_k ^ theta_k) - 
    lambda_k ^ (-theta_k) * sum(log(T_k) * T_k ^ theta_k)
}

first_deritative_lambda_k <- function(T_k, delta_k, theta_k, lambda_k){
  d_k <- sum(delta_k)
  
  - d_k * theta_k / lambda_k + theta_k * lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k)
}

second_deritative_theta_k <- function(T_k, delta_k, theta_k, lambda_k){
  d_k <- sum(delta_k)
  
  - d_k / theta_k ^ 2 -
    log(lambda_k) ^ 2 * lambda_k ^ (-theta_k) * sum(T_k ^ theta_k) +
    2 * log(lambda_k) * lambda_k ^ (-theta_k) * sum(log(T_k) * T_k ^ theta_k) -
    lambda_k ^ (-theta_k) * sum(log(T_k) ^ 2 * T_k ^ theta_k)
}

second_deritative_lambda_k <- function(T_k, delta_k, theta_k, lambda_k){
  d_k <- sum(delta_k)
  
  d_k * theta_k / lambda_k ^ 2 - theta_k * (theta_k + 1) * lambda_k ^ (-theta_k - 2) * sum(T_k ^ theta_k)
}

mix_deritative_k <- function(T_k, delta_k, theta_k, lambda_k){
  d_k <- sum(delta_k)
  
  -d_k / lambda_k +
    lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k) -
    log(lambda_k) * theta_k * lambda_k ^ (-theta_k - 1) * sum(T_k ^ theta_k) +
    theta_k * lambda_k ^ (-theta_k - 1) * sum(log(T_k) * T_k ^ theta_k)
}

grad_k <- function(T_k, delta_k, theta_k, lambda_k){
  dx_theta_k <- first_deritative_theta_k(T_k, delta_k, theta_k, lambda_k)
  dx_lambda_k <- first_deritative_lambda_k(T_k, delta_k, theta_k, lambda_k)
  
  c(dx_theta_k, dx_lambda_k)
}

Hessian_k <- function(T_k, delta_k, theta_k, lambda_k){
  H_k <- matrix(0, ncol = 2, nrow = 2)
  H_k[1,1] <- second_deritative_theta_k(T_k, delta_k, theta_k, lambda_k)
  H_k[2,2] <- second_deritative_lambda_k(T_k, delta_k, theta_k, lambda_k)
  H_k[1,2] <- mix_deritative_k(T_k, delta_k, theta_k, lambda_k)
  H_k[2,1] <- H_k[1,2]
  
  H_k
}

Newton_k_step <- function(theta_before_k, lambda_before_k, Hess_before_k, grad_before_k){
  x_before_k <- c(theta_before_k, lambda_before_k)
  x_before_k - as.vector(solve(Hess_before_k) %*% grad_before_k)
}

Newton_result <- function(T_dt, theta_0, lambda_0){
  T_v <- T_dt[, get(names(T_dt)[1])]
  delta <-  T_dt[, get(names(T_dt)[2])]
  x_k <- c(theta_0, lambda_0)
  Hess_k <- Hessian_k(T_v, delta, x_k[1], x_k[2])
  grd_k <- grad_k(T_v, delta, x_k[1], x_k[2])
  x_k2 <- Newton_k_step(x_k[1], x_k[2], Hess_k, grd_k)
  
  while((abs(x_k2[1] - x_k[1]) > 0.000005 | abs(x_k2[2] - x_k[2]) > 0.000005)){
    x_k <- c(x_k2[1], x_k2[2])
    Hess_k <- Hessian_k(T_v, delta, x_k[1], x_k[2])
    if(NaN %in% Hess_k) break
    grd_k <- grad_k(T_v, delta, x_k[1], x_k[2])
    if(NaN %in% grd_k) break
    x_k2 <- try(Newton_k_step(x_k[1], x_k[2], Hess_k, grd_k), silent = TRUE)
    if(isTRUE(is.character(x_k2[1]))) break
  }
  
  if(isTRUE(is.character(x_k2[1])) | isTRUE(is.character(x_k[1])) ){
    theta_k2 <- theta_0
    lambda_k2 <- lambda_0
  }
  else{
    theta_k2 <- x_k2[1]
    lambda_k2 <- x_k2[2]
  }
  
  c(theta_k2, lambda_k2)
}

# START POINTS OF ALGORITHM

thetha_hat_approx <- function(n, T_vals){
  function(x) n / x + sum(log(T_vals)) - (n / sum(T_vals ^ x)) * sum(log(T_vals) * (T_vals) ^ x)
}
mle_theta_full_sample <- function(n, T_vals) pracma::bisect(thetha_hat_approx(n, T_vals), 0, 10)$root

mle_lambda_full_sample <- function(n, T_vals, theta_hat_est) (1 / n * (sum(T_vals ^ theta_hat_est))) ^ (1 / theta_hat_est)


# PERCENTAGE LEVEL OF CENSURING

m_find <- function(x, pi){
  function(m){
    mean(punif(x, min = 0, max = m)) - pi
  }
}


