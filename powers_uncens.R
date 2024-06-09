library(data.table)
library(foreach)
library(doParallel)
#library(xtable)

# DISTRIBIUTION PACKAGES
library(RelDists)
library(invgamma)
library(rmutil)
library(ggamma)

sampling_uncensured <- function(X, rowNameQDF = "Uncensured_100"){
  n <- length(X)
  quantiles <- unlist(Quantiles_DF[rowNameQDF, ])
  
  T_v_dt <- data.table("T_val" = X, "delta" = 1)
  
  theta_hat <- mle_theta_full_sample(n, T_v_dt[, T_val])
  lambda_hat <- mle_lambda_full_sample(n, T_v_dt[, T_val], theta_hat)
  
  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)
  
  c(S_n_1(a = 1, n, Y_v_dt), S_n_1(a = 2, n, Y_v_dt), S_n_1(a = 5, n, Y_v_dt),
    S_n_2(a = 1, n, Y_v_dt), S_n_2(a = 2, n, Y_v_dt), S_n_2(a = 5, n, Y_v_dt), 
    KS_test_full(X, theta_hat, lambda_hat, n),CM_test_full(X, theta_hat, lambda_hat, n),
    LS_test_full(X, n)) > quantiles
}

m <-10000

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

list_of_df_amount_of_H1_wins <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{
  # H0
  W_1 <- rweibull(10, 0.5, 1)
  W_2 <- rweibull(10, 1, 1)
  W_3 <- rweibull(10, 1.5, 1)
  
  # H1, IHR
  G_1 <- rgamma(10, 2, 1)
  G_2 <- rgamma(10, 3, 1)
  AW_1 <- RelDists::rAddW(10, mu = 10, sigma = 1 / 0.02, nu = 1, tau = 5.2)
  while(0 %in% AW_1){
    AW_1 <- RelDists::rAddW(10, mu = 10, sigma = 1 / 0.02, nu = 1, tau = 5.2)
  }
  
  # H1, UBT
  LN_1 <- rlnorm(10, 0, 0.8)
  IG_1 <- invgamma::rinvgamma(10, 3, 1)
  #IG_2 <- invgamma::rinvgamma(100, 5, 2)
  EW_1 <- RelDists::rEW(10, nu = 4, mu = 1/ 12, sigma = 0.6)
  
  
  # H1, DHR
  G_3 <- rgamma(10, 0.2, 1)
  H_1 <- rmutil::rhjorth(10, s = 0, m = 1, f = 1) #f theta, s beta m delta
  AW_2 <- RelDists::rAddW(10, mu = 2, sigma = 1 / 20, nu = 1, tau = 0.1)
  while(0 %in% AW_2){
    AW_2 <- RelDists::rAddW(10, mu = 2, sigma = 1 / 20, nu = 1, tau = 0.1)
  }
  
  # H1, BT
  GG_1 <- ggamma::rggamma(10, k = 0.1, a = 1, b = 4)
  GG_2 <- ggamma::rggamma(10, k = 0.2, a = 1, b = 3)
  B_2 <- rbeta(10, shape1 = 0.1, shape2 = 3)
  
  
  data.frame("W_1" = sampling_uncensured(W_1, "Uncensured_10"), # H0 ok
             "W_2" = sampling_uncensured(W_2, "Uncensured_10"),
             "W_3" = sampling_uncensured(W_3, "Uncensured_10"),
             "IHR_G_1" = sampling_uncensured(G_1, "Uncensured_10"), # H1 ok, IHR
             "IHR_G_2" = sampling_uncensured(G_2, "Uncensured_10"),
             "IHR_AW_1" = sampling_uncensured(AW_1, "Uncensured_10"),
             "UBT_IG_1" = sampling_uncensured(IG_1, "Uncensured_10"),
             "UBT_EW_1" = sampling_uncensured(EW_1, "Uncensured_10"),
             "DHR_G_3" = sampling_uncensured(G_3, "Uncensured_10"), # H1 ok, DHR
             "DHR_H_2" = sampling_uncensured(H_1, "Uncensured_10"),
             "DHR_AW_2" = sampling_uncensured(AW_2, "Uncensured_10"),
             "BT_GG_1" = sampling_uncensured(GG_1, "Uncensured_10"), # H1 ok, BT
             "BT_GG_2" = sampling_uncensured(GG_2, "Uncensured_10"),
             "BT_B_2" = sampling_uncensured(B_2, "Uncensured_10"))
}
stopImplicitCluster()



computePower <- function(list_of_vals, distName, StatTestName){
  list_of_col_value <- lapply(list_of_vals, function(df) df[distName]) 
  list_of_statistic_value <- unlist(lapply(list_of_col_value, function(df) df[StatTestName, ]))
  
  power <- mean(list_of_statistic_value)
  power
}

colNames_vec <- colnames(list_of_df_amount_of_H1_wins[[1]])



Power_df <- sapply(colNames_vec,
                   function(colName){
                     lapply(colnames(Quantiles_DF10), 
                            function(rowName){
                              computePower(list_of_df_amount_of_H1_wins,
                                           distName = colName, 
                                           StatTestName = rowName)})})

rownames(Power_df) <- colnames(Quantiles_DF)

Power_df<- t(Power_df)

print(xtable(Power_df, type = "latex"), file = "...")


#hist(rmutil::rhjorth(10000, s = 0, m = 1, f = 1), breaks = 100)
#hist(rweibull(10000, 2,1.4), breaks = 100)
#load("/Users/zuza/Downloads/Sample_vals_10n.rda")