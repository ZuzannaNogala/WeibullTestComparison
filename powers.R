sampling_uncensured <- function(X, rowNameQDF = "Uncensured_50"){
  n <- length(X)
  quantiles <- unlist(Quantiles_DF[rowNameQDF, ])
  
  T_v_dt <- data.table("T_val" = X, "delta" = 1)
  
  theta_hat <- mle_theta_full_sample(n, T_v_dt[, T_val])
  lambda_hat <- mle_lambda_full_sample(n, T_v_dt[, T_val], theta_hat)
  
  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)
  sort_Y_v_dt <- sort_Y_vals(Y_v_dt)
  
  c(S_n_1(a = 1, n, Y_v_dt), S_n_1(a = 2, n, Y_v_dt), S_n_1(a = 5, n, Y_v_dt),
    S_n_2(a = 1, n, Y_v_dt), S_n_2(a = 2, n, Y_v_dt), S_n_2(a = 5, n, Y_v_dt), 
    KS_test(Y_v_dt, n),CM_test(T_v_dt, X, theta_hat, lambda_hat, n),
    LS_test(Y_v_dt, n)) > quantiles
}

sampling_cenusured <- function(X, perc_censored = 0.1, rowNameQDF = "Censured_10perc_50"){
  n <- length(X)
  
  quantiles <- unlist(Quantiles_DF[rowNameQDF, ])
  
  m_thres <- pracma::bisect(m_find(X, perc_censored), 0, 10000)$root
  C <- runif(n, 0, m_thres)
  
  T_vec <- sapply(1:n, function(i) min(X[i], C[i]))
  delta <- T_vec == X
  T_v_dt <- data.table("T_val" = T_vec, "delta" = delta)
  
  theta_hat_0 <- mle_theta_full_sample(n, T_v_dt[, T_val])
  lambda_hat_0 <- mle_lambda_full_sample(n, T_v_dt[, T_val], theta_hat_0)
  
  params <- Newton_result(T_v_dt, theta_hat_0, lambda_hat_0) 
  theta_hat <- params[1]
  lambda_hat <- params[2]
  
  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)
  
  c(S_n_1(a = 1, n, Y_v_dt), S_n_1(a = 2, n, Y_v_dt), S_n_1(a = 5, n, Y_v_dt),
    S_n_2(a = 1, n, Y_v_dt), S_n_2(a = 2, n, Y_v_dt), S_n_2(a = 5, n, Y_v_dt), 
    KS_test(Y_v_dt, n),CM_test(T_v_dt, X, theta_hat, lambda_hat, n),
    LS_test(Y_v_dt, n)) > quantiles
}

sampling_uncensured(rweibull(50, 5, 2), "Uncensured_100")
sampling_cenusred(rgamma(50, 4, 1))




sapply(rownames(Quantiles_DF), function(i) sampling_uncensured(W_1, i))
lapply(rownames(Quantiles_DF)[3:4], function(i) sapply(c(0.1, 0.2), function(perc) sampling_cenusured(W_1, perc_censored = perc, i))) 

m <- 20

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

X_vals <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{
  W_1 <- rweibull(50, 0.5, 1)
  W_2 <- rweibull(50, 1, 1)
  W_3 <- rweibull(50, 1.5, 1)
  
  G_1 <- rgamma(50, 2, 1)
  G_2 <- rgamma(50, 3, 1)
  AW_1 <- RelDists::rAddW(50, 10, 1 / 0.02, 1, 5.2)
  
  LN_1 <- rlnorm(50, 0, 0.8)
  IG_1 <- invgamma::rinvgamma(50, 3, 1)
  IS_1 <- rmutil::rinvgauss(50, 1, 0.25)
  
  G_3 <- rgamma(50, 0.2, 1)
  H_1 <- rmutil::rhjorth(50, 1, 0, 1)
  EW_1 <- RelDists::dEW(50, 1/ 0.01, 0.95, 0.1)
  
  GG_1 <- ggamma::rggamma()
  
  
  data.frame("W_1" = sampling_uncensured(W_1),
             "W_2" = sampling_uncensured(W_2),
             "W_3" = sampling_uncensured(W_3),
             "W_1_10perc" = sampling_cenusured(W_1, perc_censored = 0.1, "Censured_10perc_50"),
             "W_2_10perc" = sampling_cenusured(W_2, perc_censored = 0.1, "Censured_10perc_50"),
             "W_3_10perc" = sampling_cenusured(W_3, perc_censored = 0.1, "Censured_10perc_50"),
             "W_1_20perc" = sampling_cenusured(W_1, perc_censored = 0.2, "Censured_20perc_50"),
             "W_2_20perc" = sampling_cenusured(W_2, perc_censored = 0.2, "Censured_20perc_50"),
             "W_3_20perc" = sampling_cenusured(W_3, perc_censored = 0.2, "Censured_20perc_50"),
             )
}
stopImplicitCluster()

X <- sapply(1:20, function(i) rweibull(50, 0.5, 1))
apply(X, MARGIN = 2, sampling_cenusured)

X <- sapply(1:100, function(i) rweibull(50, 1, 1))
X <- sapply(1:100, function(i) rweibull(50, 1.5, 1))


hist(rmutil::rinvgauss(1000, 1, 0.25), xlim = c(0, 5), breaks = 100)

?RelDists::dEW()
