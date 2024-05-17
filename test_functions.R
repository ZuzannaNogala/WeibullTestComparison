# TEST NEWTON-RAPHSON ALGORITHM

MLE_params_estimator_Weibull <- function(X, perc_censored = 0){
  if(perc_censored == 0){
    C <- X
    T_vec <- X
    delta <- rep(1, length(X))
    T_v_dt <- data.table("T_val" = T_vec, "delta" = delta)
    theta_hat <- mle_theta_full_sample(X)
    lambda_hat <- mle_lambda_full_sample(T_vec, theta_hat_est)
  }
  
  else{
    m_thres <- pracma::bisect(m_find(X, perc_censored), 0, 10000)$root
    C <- runif(n, 0, m_thres) 
    
    T_vec <- sapply(1:length(X), function(i) min(X[i], C[i]))
    delta <- T_vec == X
    T_v_dt <- data.table("T_val" = T_vec, "delta" = delta)
    
    theta_hat_0 <- mle_theta_full_sample(X)
    lambda_hat_0 <- mle_lambda_full_sample(T_vec, theta_hat_0)
    
    params <- Newton_result(T_v_dt, theta_hat_0, lambda_hat_0) 
    theta_hat <- params[1]
    lambda_hat <- params[2]
  }
  
  return(c(theta_hat, lambda_hat))
}

X <- rweibull(n, 10, 9)
is.character(MLE_params_estimator_Weibull(X, 0.1))
is.character(MLE_params_estimator_Weibull(X, 0.4))
MLE_params_estimator_Weibull(X)