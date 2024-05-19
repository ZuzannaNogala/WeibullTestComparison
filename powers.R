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
  
  c(S_n_1(a = 0.75, n, Y_v_dt), S_n_1(a = 1, n, Y_v_dt), S_n_1(a = 2, n, Y_v_dt),
    S_n_2(a = 0.75, n, Y_v_dt), S_n_2(a = 1, n, Y_v_dt), S_n_2(a = 2, n, Y_v_dt), 
    KS_test(Y_v_dt, n),CM_test(T_v_dt, X, theta_hat, lambda_hat, n),
    LS_test(Y_v_dt, n)) > quantiles
}

m <- 800

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

Sample_vals_2 <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{
  # H0
  W_1 <- rweibull(50, 0.5, 1)
  W_2 <- rweibull(50, 1, 1)
  W_3 <- rweibull(50, 1.5, 1)
  
  # H1, IHR
  G_1 <- rgamma(50, 2, 1)
  G_2 <- rgamma(50, 3, 1)
  AW_1 <- RelDists::rAddW(50, mu = 10, sigma = 1 / 0.02, nu = 1, tau = 5.2)
  while(0 %in% AW_1){
    AW_1 <- RelDists::rAddW(50, mu = 10, sigma = 1 / 0.02, nu = 1, tau = 5.2)
  }
  
  # H1, UBT
  LN_1 <- rlnorm(50, 0, 0.8)
  IG_1 <- invgamma::rinvgamma(50, 3, 1)
  #IG_2 <- invgamma::rinvgamma(50, 5, 2)
  EW_1 <- RelDists::rEW(50, nu = 4, mu = 1/ 12, sigma = 0.6)
  
  
  # H1, DHR
  G_3 <- rgamma(50, 0.2, 1)
  H_1 <- rmutil::rhjorth(50, s = 0, m = 1, f = 1) #f theta, s beta m delta
  AW_2 <- RelDists::rAddW(50, mu = 2, sigma = 1 / 20, nu = 1, tau = 0.1)
  while(0 %in% AW_2){
    AW_2 <- RelDists::rAddW(50, mu = 2, sigma = 1 / 20, nu = 1, tau = 0.1)
  }
  
  # H1, BT
  GG_1 <- ggamma::rggamma(50, k = 0.1, a = 1, b = 4)
  GG_2 <- ggamma::rggamma(50, k = 0.2, a = 1, b = 3)
  #EW_2 <- RelDists::rEW(50, nu = 0.1, mu = 1/ 100, sigma = 5)
  B_2 <- rbeta(50, shape1 = 0.1, shape2 = 3)
  
  
  data.frame("W_1" = sampling_uncensured(W_1), # H0 ok
             "W_2" = sampling_uncensured(W_2),
             "W_3" = sampling_uncensured(W_3),
             "W_1_10perc" = sampling_cenusured(W_1, perc_censored = 0.1, "Censured_10perc_50"),
             "W_2_10perc" = sampling_cenusured(W_2, perc_censored = 0.1, "Censured_10perc_50"),
             "W_3_10perc" = sampling_cenusured(W_3, perc_censored = 0.1, "Censured_10perc_50"),
             "W_1_20perc" = sampling_cenusured(W_1, perc_censored = 0.2, "Censured_20perc_50"),
             "W_2_20perc" = sampling_cenusured(W_2, perc_censored = 0.2, "Censured_20perc_50"),
             "W_3_20perc" = sampling_cenusured(W_3, perc_censored = 0.2, "Censured_20perc_50"),
             "IHR_G_1" = sampling_uncensured(G_1), # H1 ok, IHR
             "IHR_G_2" = sampling_uncensured(G_2),
             "IHR_AW_1" = sampling_uncensured(AW_1),
             "IHR_G_1_10perc" = sampling_cenusured(G_1, perc_censored = 0.1, "Censured_10perc_50"),
             "IHR_G_2_10perc" = sampling_cenusured(G_2, perc_censored = 0.1, "Censured_10perc_50"),
             "IHR_AW_1_10perc" = sampling_cenusured(AW_1, perc_censored = 0.1, "Censured_10perc_50"),
             "IHR_G_1_20perc" = sampling_cenusured(G_1, perc_censored = 0.2, "Censured_20perc_50"),
             "IHR_G_2_20perc" = sampling_cenusured(G_2, perc_censored = 0.2, "Censured_20perc_50"),
             "IHR_AW_1_20perc" = sampling_cenusured(AW_1, perc_censored = 0.2, "Censured_20perc_50"),
             "UBT_LN_1" = sampling_uncensured(LN_1), # H1 ok, UBT
             "UBT_IG_1" = sampling_uncensured(IG_1),
             "UBT_EW_1" = sampling_uncensured(EW_1),
             "UBT_LN_1_10perc" = sampling_cenusured(LN_1, perc_censored = 0.1, "Censured_10perc_50"),
             "UBT_IG_1_10perc" = sampling_cenusured(IG_1, perc_censored = 0.1, "Censured_10perc_50"),
             "UBT_EW_1_10perc" = sampling_cenusured(EW_1, perc_censored = 0.1, "Censured_10perc_50"),
             "UBT_LN_1_20perc" = sampling_cenusured(LN_1, perc_censored = 0.2, "Censured_20perc_50"),
             "UBT_IG_1_20perc" = sampling_cenusured(IG_1, perc_censored = 0.2, "Censured_20perc_50"),
             "UBT_EW_1_20perc" = sampling_cenusured(EW_1, perc_censored = 0.2, "Censured_20perc_50"),
             "DHR_G_3" = sampling_uncensured(G_3), # H1 ok, DHR
             "DHR_H_2" = sampling_uncensured(H_1),
             "DHR_AW_2" = sampling_uncensured(AW_2),
             "DHR_G_3_10perc" = sampling_cenusured(G_3, perc_censored = 0.1, "Censured_10perc_50"),
             "DHR_H_2_10perc" = sampling_cenusured(H_1, perc_censored = 0.1, "Censured_10perc_50"),
             "DHR_AW_2_10perc" = sampling_cenusured(AW_2, perc_censored = 0.1, "Censured_10perc_50"),
             "DHR_G_3_20perc" = sampling_cenusured(G_3, perc_censored = 0.2, "Censured_20perc_50"),
             "DHR_H_2_20perc" = sampling_cenusured(H_1, perc_censored = 0.2, "Censured_20perc_50"),
             "DHR_AW_2_20perc" = sampling_cenusured(AW_2, perc_censored = 0.2, "Censured_20perc_50"),
             "BT_GG_1" = sampling_uncensured(GG_1), # H1 ok, BT
             "BT_GG_2" = sampling_uncensured(GG_2),
             "BT_B_2" = sampling_uncensured(B_2),
             "BT_GG_1_10perc" = sampling_cenusured(GG_1, perc_censored = 0.1, "Censured_10perc_50"),
             "BT_GG_2_10perc" = sampling_cenusured(GG_2, perc_censored = 0.1, "Censured_10perc_50"),
             "BT_B_2_10perc" = sampling_cenusured(B_2, perc_censored = 0.1, "Censured_10perc_50"),
             "BT_GG_1_20perc" = sampling_cenusured(GG_1, perc_censored = 0.2, "Censured_20perc_50"),
             "BT_GG_2_20perc" = sampling_cenusured(GG_2, perc_censored = 0.2, "Censured_20perc_50"),
             "BT_B_2_20perc" = sampling_cenusured(B_2, perc_censored = 0.2, "Censured_20perc_50"))
}
stopImplicitCluster()

#Sample_vals <- c(Sample_vals_2, Sample_vals)

save(Sample_vals_2, file = "/Users/zuza/Desktop/studia/licencjat/SimulationData/licencjat_R/Sample_vals_2.rda")

Sample_vals <- c(Sample_vals, Sample_vals_2)

load("/Users/maniek/Desktop/licencjat_R/Sample_vals.rda")
load("/Users/maniek/Desktop/licencjat_R/list_of.rda")

load("/Users/zuza/Desktop/studia/licencjat/SimulationData/licencjat_R/Sample_vals.rda")

colNames_vec <- colnames(Sample_vals[[1]])

computePower <- function(distName, StatTestName){
  list_of_col_value <- lapply(Sample_vals, function(df) df[distName])
  list_of_statistic_value <- unlist(lapply(list_of_col_value, function(df) df[StatTestName, ]))
  
  power <- mean(list_of_statistic_value)
  power
}


Power_df <- sapply(colNames_vec,
                   function(colName) lapply(colnames(Quantiles_DF), 
                                            function(rowName) computePower(distName = colName, StatTestName = rowName)))

rownames(Power_df) <- colnames(Quantiles_DF)

Power_df <- t(Power_df)

Power_df[1:3, ] # Weibull, full sample

Power_df[4:9, ] # Weibull, censured 10% and 20%

Power_df[10:12, ] # IHR, full sample

Power_df[13:18, ] # IHR, censured 10% and 20%

Power_df[19:21, ] # UBT, full sample

Power_df[22:27, ] # UBT, censured 10% and 20%

Power_df[28:30, ] # DHR, full sample

Power_df[31:36, ] # DHR, censured 10% and 20%

Power_df[37:39, ] # BT, full sample

Power_df[40:45, ] # BT, censured 10% and 20%
