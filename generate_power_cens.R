library(data.table)
library(foreach)
library(doParallel)

# DISTRIBIUTION PACKAGES
library(RelDists)
library(invgamma)
library(rmutil)
library(ggamma)


bootstrap_sample_test_value <- function(n, X, C, theta_hat, lambda_hat){
  X_bstr <- rweibull(n, theta_hat, lambda_hat)
  
  C_bstr <- sample(C, n, replace = TRUE)
  
  T_bstr_vec <- sapply(1:n, function(i) min(X_bstr[i], C_bstr[i]))
  delta_bstr <- T_bstr_vec == X_bstr
  
  T_bstr_v_dt <- data.table("T_val" = T_bstr_vec, "delta" = delta_bstr)
  
  theta_hat_0 <- mle_theta_full_sample(n, T_bstr_v_dt[, T_val])
  lambda_hat_0 <- mle_lambda_full_sample(n, T_bstr_v_dt[, T_val], theta_hat_0)
  
  params_bstr <- Newton_result(T_bstr_v_dt, theta_hat_0, lambda_hat_0)
  theta_hat_bstr <- params_bstr[1]
  lambda_hat_bstr <- params_bstr[2]
  
  Y_bstr_v_dt <- transform_to_Y(T_bstr_v_dt, theta_hat_bstr, lambda_hat_bstr)
  
  # diff a params! 
  
  list("S_n_1_1" =  S_n_1(a = 0.75, n, Y_bstr_v_dt),
       "S_n_1_2" =  S_n_1(a = 1, n, Y_bstr_v_dt),
       "S_n_1_5" =  S_n_1(a = 2, n, Y_bstr_v_dt),
       "S_n_2_1" =  S_n_2(a = 0.75, n, Y_bstr_v_dt),
       "S_n_2_2" =  S_n_2(a = 1, n, Y_bstr_v_dt),
       "S_n_2_5" =  S_n_2(a = 2, n, Y_bstr_v_dt),
       "KS" = KS_test(Y_bstr_v_dt, n),
       "CM" = CM_test(Y_bstr_v_dt, n),
       "LS" = LS_test(Y_bstr_v_dt, n))
}

H1_wins_check <- function(X, cens_perc, n){
  m_thres <- pracma::bisect(m_find(X, cens_perc), 0, 10000)$root 
  C <- runif(n, 0, m_thres)
  
  T_v_dt <- data.table("T_val" = pmin(X,C), "delta" = pmin(X,C) == X)
  
  theta_hat_0 <- mle_theta_full_sample(n, T_bstr_v_dt[, T_val])
  lambda_hat_0 <- mle_lambda_full_sample(n, T_bstr_v_dt[, T_val], theta_hat_0)
  
  params_adjust <- Newton_result(T_bstr_v_dt, theta_hat_0, lambda_hat_0)
  theta_hat <- params_adjust[1]
  lambda_hat <- params_adjust[2]
  
  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)
  
  values_of_statistics <- c(S_n_1(a = 0.75, n, Y_v_dt), S_n_1(a = 1, n, Y_v_dt), S_n_1(a = 2, n, Y_v_dt),
                            S_n_2(a = 0.75, n, Y_v_dt), S_n_2(a = 1, n, Y_v_dt), S_n_2(a = 2, n, Y_v_dt), 
                            KS_test(Y_v_dt, n), CM_test(Y_v_dt, n),
                            LS_test(Y_v_dt, n)) 
  
  values_of_statistics_bstr <- sapply(1:1000, function(i) bootstrap_sample_test_value(n, X, C, theta_hat, lambda_hat))
  
  quantiles <- apply(values_of_statistics_bstr, MARGIN = 1, function(vals) sort(unlist(vals))[1000 * (0.9)])
  
  values_of_statistics > quantiles
} 

m <- 1000
n <- 100
cens_perc <- 0.1

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

computer_values_to_powers_cens10 <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{
  W_1 <- rweibull(n, 0.5, 1)
  W_2 <- rweibull(n, 1, 1)
  W_3 <- rweibull(n, 1.5, 1)
  
  # H1, IHR
  G_1 <- rgamma(n, 2, 1)
  G_2 <- rgamma(n, 3, 1)
  AW_1 <- RelDists::rAddW(n, mu = 10, sigma = 1 / 0.02, nu = 1, tau = 5.2)
  while(0 %in% AW_1){
    AW_1 <- RelDists::rAddW(n, mu = 10, sigma = 1 / 0.02, nu = 1, tau = 5.2)
  }
  
  # H1, UBT
  LN_1 <- rlnorm(n, 0, 0.8)
  IG_1 <- invgamma::rinvgamma(n, 3, 1)
  EW_1 <- RelDists::rEW(n, nu = 4, mu = 1/ 12, sigma = 0.6)
  
  
  # H1, DHR
  G_3 <- rgamma(n, 0.2, 1)
  H_1 <- rmutil::rhjorth(n, s = 0, m = 1, f = 1) #f theta, s beta m delta
  AW_2 <- RelDists::rAddW(n, mu = 2, sigma = 1 / 20, nu = 1, tau = 0.1)
  while(0 %in% AW_2){
    AW_2 <- RelDists::rAddW(n, mu = 2, sigma = 1 / 20, nu = 1, tau = 0.1)
  }
  
  # H1, BT
  GG_1 <- ggamma::rggamma(n, k = 0.1, a = 1, b = 4)
  GG_2 <- ggamma::rggamma(n, k = 0.2, a = 1, b = 3)
  B_2 <- rbeta(n, shape1 = 0.1, shape2 = 3)
  
  data.frame("W_1" = H1_wins_check(W_1, cens_perc, n),
             "W_2" = H1_wins_check(W_2, cens_perc, n),
             "W_3" = H1_wins_check(W_3, cens_perc, n),
             "G_1" = H1_wins_check(G_1, cens_perc, n),
             "G_2" = H1_wins_check(G_2, cens_perc, n),
             "AW_1" = H1_wins_check(AW_1, cens_perc, n),
             "LN_1" = H1_wins_check(LN_1, cens_perc, n),
             "IG_1" = H1_wins_check(IG_1, cens_perc, n),
             "EW_1" = H1_wins_check(EW_1, cens_perc, n),
             "G_3" = H1_wins_check(G_3, cens_perc, n),
             "H_1" = H1_wins_check(H_1, cens_perc, n),
             "AW_2" = H1_wins_check(AW_2, cens_perc, n),
             "GG_1" = H1_wins_check(GG_1, cens_perc, n),
             "GG_2" = H1_wins_check(GG_2, cens_perc, n),
             "B_2" = H1_wins_check(B_2, cens_perc, n)
  )
}
stopImplicitCluster()


computePower <- function(distName, StatTestName){
  list_of_col_value <- lapply(computer_values_to_powers_cens10, function(df) df[distName]) 
  
  list_of_statistic_value <- unlist(lapply(list_of_col_value, function(df) df[StatTestName, ]))
  
  power <- mean(list_of_statistic_value)
  power
}

distNames <- colnames(computer_values_to_powers_cens10[[1]])
statsNames <- rownames(computer_values_to_powers_cens10[[1]])

powers_dt_cens10 <- sapply(distNames, function(colName) lapply(statsNames, function(rowName){
  computePower(distName = colName, StatTestName = rowName)
}))

rownames(powers_dt_cens10) <- statsNames


