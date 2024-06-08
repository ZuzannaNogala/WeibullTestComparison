library(data.table)

# LEAST SQUARE ESTIMATORS TO LS TEST
library(EWGoF)

transform_to_Y <- function(T_vec_dt, theta_hat, lambda_hat){
  T_vec_dt[, .(Y_val = theta_hat * (log(T_val) - log(lambda_hat)) , delta)]
}

sort_Y_vals <- function(Y_vec_dt){
  Y_values <- Y_vec_dt[, Y_val]
  names(Y_values) <- Y_vec_dt[, delta]
  sort_Y_values <- sort(Y_values)
  if(sum(Y_vec_dt[, delta]) != length(Y_values)){
    sort_delta <- as.logical(names(sort_Y_values))
  }
  else{
    sort_delta <- rep(TRUE, length(Y_values))
  }
  dt <- data.table(sort_Y_val = sort_Y_values, delta = sort_delta)
  dt
}

Est_kaplan_meier <- function(t, n, sort_vec_dt){ # lack of >= if num1, lack k -1  
  Vec <- sort_vec_dt[, get(names(sort_vec_dt)[1])]
  delta_Vec <- sort_vec_dt[, get(names(sort_vec_dt)[2])]
  if(Vec[1] > t) return(0)
  if(Vec[n] < t) return(1)
  else{
    k <- which.max(Vec >= t)
    k_before <- k 
    expr1 <- prod(((n - 1:k_before) / (n + 1 - 1:k_before)) ^ delta_Vec[1:k_before])
    return(1 - expr1)
  }
}

Delta_k <- function(k, n, sort_Y_vec_dt){
  delta_vec <- sort_Y_vec_dt[, delta]
  if(k == 1) delta_vec[1] / n
  if(k == n) prod(((n - 1 : (n - 1)) / (n - 1 : (n - 1) + 1)) ^ delta_vec[1 : (n - 1)])
  else (delta_vec[k] / (n - k + 1)) * prod(((n - 1 : (k - 1)) / (n - 1 : (k - 1) + 1)) ^ delta_vec[1 : (k - 1)])
}

# NEW TESTS

S_n_1 <- function(a, n, Y_dt){
  Y_val <- sort(Y_dt[, Y_val])

  Delta <- sapply(1:n, function(k) Delta_k(k, n, sort_Y_vals(Y_dt)))
  
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
  Y_val <- sort(Y_dt[, Y_val])
 
  Delta <- sapply(1:n, function(k) Delta_k(k, n, sort_Y_vals(Y_dt)))

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

# OTHER TESTS FOR CENSURED DATA

KS_test <- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[delta == 1, sort_Y_val]
  Est_KM_sort_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  expr1 <- 1 - exp(-exp(sort_Y_values))
    
  S1 <- max(Est_KM_sort_Y - expr1)
  S2 <- max(expr1 - Est_KM_sort_Y)

  max(S1, S2)
}

CM_test <- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[, sort_Y_val]
  num_uncens <- sort_Y_vec_dt[, which(delta == 1)]
  d <- sort_Y_vec_dt[delta == 1, .N]

  Est_KM_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  
  expr1 <- sort_Y_values[num_uncens][2:d] - sort_Y_values[num_uncens][1:(d-1)]
  expr2 <- sort_Y_values[num_uncens][2:d]^2 - sort_Y_values[num_uncens][1:(d-1)]^2
  
  expr_full <- Est_KM_Y[num_uncens][2:d]^2 * expr1 - Est_KM_Y[num_uncens][2:d] * expr2

  sum(expr_full) * n + n / 3
}

LS_test <- function(Y_vec_dt, n){
  sort_Y_vec_dt <- sort_Y_vals(Y_vec_dt)
  sort_Y_values <- sort_Y_vec_dt[, sort_Y_val]
  Est_KM_Y <- sapply(sort_Y_values, function(t) Est_kaplan_meier(t, n, sort_Y_vec_dt))
  
  d <- length(sort_Y_values)
  
  max_S1 <- sort_Y_vec_dt[, .("sort_Y_val" = 1:n / n * delta, delta)][delta == 1, sort_Y_val] - Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]]
  max_S2 <- Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]] - sort_Y_vec_dt[, .("sort_Y_val" = 0:(n-1) / n * delta, delta)][delta == 1, sort_Y_val]
  max_S <- pmax(max_S1, max_S2)
  
  expr1 <- sqrt(Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]] * (1 - Est_KM_Y[sort_Y_vec_dt[, which(delta == 1)]]))
  
  1 / sqrt(n) * sum(ifelse(max_S / expr1 == Inf, 0, max_S / expr1))
}

# OTHER TESTS FOR FULL SAMPLE, WITHOUT CENSORING

KS_test_full <- function(X_vec, theta_hat, lambda_hat, n){
  U <- 1 - exp( -exp(theta_hat * (log(X_vec) - log(lambda_hat))))
  U_sort <- sort(U)
  expr1 <- max(1:n / n - U_sort)
  expr2 <- max(U_sort - 0:(n - 1) / n)
  
  sqrt(n) * max(expr1, expr2)
}

CM_test_full <- function(X_vec, theta_hat, lambda_hat, n){
  U <- 1 - exp( -exp(theta_hat * (log(X_vec) - log(lambda_hat))))
  U_sort <- sort(U)
  
  sum((U_sort - (2 * (1:n) - 1) / (2 * n))^2) + 1 / (12 * n)
}

LS_test_full <- function(X_vec, n){
  lambda_hat_LS <- EWGoF::LSEst(X_vec)$eta
  theta_hat_LS <- EWGoF::LSEst(X_vec)$beta
  
  U <- 1 - exp( -exp(theta_hat_LS * (log(X_vec) - log(lambda_hat_LS))))
  U_sort <- sort(U)
  
  expr_1 <- (1:n) / n - U_sort
  expr_2 <- U_sort - (0:(n - 1)) / n
  
  1 / sqrt(n) * sum(pmax(expr_1, expr_2) / sqrt(U_sort * (1 - U_sort)))
}

