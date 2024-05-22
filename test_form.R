transform_to_Y <- function(T_vec_dt, theta_hat, lambda_hat){
  T_vec_dt[, .(Y_val = theta_hat * (log(T_val) - log(lambda_hat)) , delta)]
}

sort_Y_vals <- function(Y_vec_dt){
  Y_values <- Y_vec_dt[, Y_val]
  names(Y_values) <- Y_vec_dt[, delta]
  sort_Y_values <- sort(Y_values)
  if(sum(Y_vec_dt[, delta]) != length(Y_values)){
    sort_delta <- sort_delta <- as.logical(names(sort_Y_values))
  }
  else{
    sort_delta <- rep(TRUE, length(Y_values))
  }
  dt <- data.table(sort_Y_val = sort_Y_values, delta = sort_delta)
  dt
}


sort_T_vals <- function(T_vec_dt){
  T_values <- T_vec_dt[, T_val]
  names(T_values) <-T_vec_dt[, delta]
  sort_T_values <- sort(T_values)
  sort_delta <- as.logical(names(sort_T_values))
  
  dt <- data.table(sort_T_val = sort_T_values, delta = sort_delta)
  dt
}

sort_C_vals_bootstrap <- function(C_vec_dt){
  C_values <- C_vec_dt[, C_val]
  names(C_values) <- C_vec_dt[, delta]
  sort_C_values <- sort(C_values)
  sort_delta <- as.logical(names(sort_C_values))
  
  dt <- data.table(sort_C_val = sort_C_values, delta = sort_delta)
  dt
}


Est_kaplan_meier <- function(t, n, sort_vec_dt){
  Vec <- sort_vec_dt[, get(names(sort_vec_dt)[1])]
  delta_Vec <- sort_vec_dt[, get(names(sort_vec_dt)[2])]
  if(Vec[1] >= t) return(0)
  if(Vec[n] < t) return(1)
  else{
    k <- which.max(Vec >= t)
    k_before <- k - 1
    expr1 <- prod(((n - 1:k_before) / (n + 1 - 1:k_before)) ^ delta_Vec[1:k_before])
    return(1 - expr1)
  }
}

# sapply(sort(theta_hat_est * ((log(X) - log(lambda_hat_est)))),
#        function(t) Est_kaplan_meier(t, 100, sort_Y_v_dt))

# fun_quantile_Est_kaplan_meier <- function(n, prob, sort_vec_dt){
#   function(t) Est_kaplan_meier(t, n, sort_vec_dt) - prob
# }
# 
# find_quantile_kaplan_meier <- function(n, prob, sort_vec_dt){
#   fun <- fun_quantile_Est_kaplan_meier(n, prob, sort_vec_dt)
#   a <-  min(sort_vec_dt[, get(names(sort_vec_dt)[1])])
#   b <-  max(sort_vec_dt[, get(names(sort_vec_dt)[1])]) + 1
#   quantile <- pracma::bisect(fun, a, b)$root
#   quantile
# }
# 
# non_param_C_bootstrap <- function(n, C, delta_list, theta_hat, lambda_hat){
#   probs <- runif(n)
#   C_transf_v_dt <- data.table("C_val" = theta_hat * (log(C) - log(lambda_hat)),
#                               "delta" = unlist(delta_list))
#   sort_C_transf_v_dt <- sort_C_vals_bootstrap(C_transf_v_dt)
#   C_bootstrap <- sapply(probs, function(p) find_quantile_kaplan_meier(n, p, sort_C_transf_v_dt))
#   C_stars <- exp(C_bootstrap / theta_hat + log(lambda_hat))
#   
#   C_stars
# }

Delta_k <- function(k, n, sort_Y_vec_dt){
  delta_vec <- sort_Y_vec_dt[, delta]
  if(k == 1) delta_vec[1] / n
  if(k == n) prod(((n - 1 : (n - 1)) / (n - 1 : (n - 1) + 1)) ^ delta_vec[1 : (n - 1)])
  else (delta_vec[k] / (n - k + 1)) * prod(((n - 1 : (k - 1)) / (n - 1 : (k - 1) + 1)) ^ delta_vec[1 : (k - 1)])
}

S_n_1 <- function(a, n, Y_dt){
  Delta <- sapply(1:n, function(k) Delta_k(k, n, sort_Y_vals(Y_dt)))
  
  Y_val <- Y_dt[, Y_val]
  S_j <- sapply(1:n, function(i){
    Y_j <- Y_val[i]
    S_1 <- exp((- (Y_j - Y_val) ^ 2) / (4 * a))
    S_2 <- - ((Y_j - Y_val) ^ 2 - 2 * a) / (4 * a ^ 2)
    S_3 <- 2 * (1 - exp(Y_j)) * (Y_j - Y_val) / (2 * a) + (1 - exp(Y_j)) * (1 - exp(Y_val))
    
    sum(Delta * (S_1 * (S_2 + S_3)))
  })
  
  n * sqrt(pi / a) * sum(Delta * S_j)
}

S_n_2 <- function(a, n, Y_dt){
  Delta <- sapply(1:n, function(k) Delta_k(k, n, sort_Y_vals(Y_dt)))
  Y_val <- Y_dt[, Y_val]
  S_j <- sapply(1:n, function(i){
    Y_j <- Y_val[i]
    S_1 <- (3 * (Y_j - Y_val) ^ 2 - a ^ 2) * (-4 * a)
    S_2 <- (Y_j - Y_val) ^ 2 + a ^ 2
    S_3 <- 8 * a * (Y_j - Y_val) * (1 - exp(Y_j))
    S_4 <- 2 * a * (1 - exp(Y_j)) * (1 - exp(Y_val))
    
    sum(Delta * (S_1 / (S_2) ^ 3 + S_3 / (S_2) ^ 2 + S_4 / S_2))
  })
  
  n * sum(Delta * S_j)
}


KS_test <- function(Y_vec_dt, n){
  #Y_values <- Y_vec_dt[, Y_val]
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[, sort_Y_val]
  Est_KM_sort_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  expr1 <- sapply(sort_Y_values, function(y) 1 - exp(-exp(y)))
  
  S1 <- max(Est_KM_sort_Y - expr1)
  S2 <- max(expr1 - Est_KM_sort_Y)
  
  max(S1, S2)
}

CM_test <- function(T_vec_dt, X_vec, theta_hat, lambda_hat, n){ # add sort 
  X_t <- theta_hat * (log(sort(X_vec)) - log(lambda_hat))
  d <- sum(T_vec_dt[, delta])
  sort_Y_vec_dt <- sort_Y_vals(transform_to_Y(T_vec_dt, theta_hat, lambda_hat))
  Est_KM_X_t <- sapply(X_t, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  
  if(d == n){
    expr1 <- sapply(2:n, function(i) X_t[i] - X_t[i - 1])
    expr2 <- sapply(2:n , function(i) X_t[i] + X_t[i - 1])
    S <- sum(Est_KM_X_t[1:(n-1)] * expr1 * (Est_KM_X_t[1:(n-1)] - expr2))
  }
  else{
    expr1 <- sapply(2:(d+1), function(i) X_t[i] - X_t[i - 1])
    expr2 <- sapply(2:(d+1), function(i) X_t[i] + X_t[i - 1])
    S <- sum(Est_KM_X_t[1:d] * expr1 * (Est_KM_X_t[1:d] - expr2))
  }
  
  n / 3 + n * S
}

LS_test <- function(Y_vec_dt, n){
  #Y_values <- Y_vec_dt[, Y_val]
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[, sort_Y_val]
  Est_KM_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  max_S1 <- 1:n / n - Est_KM_Y
  max_S2 <- Est_KM_Y - 0:(n-1) / n
  max_S <- sapply(1:n, function(j) max(max_S1[j], max_S2[j]))
  
  expr1 <- sqrt(Est_KM_Y * (1 - Est_KM_Y))
  
  1 / sqrt(n) * sum(ifelse(max_S / expr1 == Inf, 0, max_S / expr1))
}
