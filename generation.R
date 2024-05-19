library(data.table)
library(foreach)
library(doParallel)

compute_test_values_uncens <- function(n){
  X <- rweibull(n, 1.5, 1)
  
  T_v_dt <- data.table("T_val" = X, "delta" = 1)
  
  theta_hat <- mle_theta_full_sample(n, T_v_dt[, T_val])
  lambda_hat <- mle_lambda_full_sample(n, T_v_dt[, T_val], theta_hat)
  
  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)
  
  list("S_n_1_1" =  S_n_1(a = 1, n, Y_v_dt),
       "S_n_1_2" =  S_n_1(a = 2, n, Y_v_dt),
       "S_n_1_5" =  S_n_1(a = 5, n, Y_v_dt),
       "S_n_2_1" =  S_n_2(a = 1, n, Y_v_dt),
       "S_n_2_2" =  S_n_2(a = 2, n, Y_v_dt),
       "S_n_2_5" =  S_n_2(a = 5, n, Y_v_dt),
       "KS" = KS_test(Y_v_dt, n),
       "CM" = CM_test(T_v_dt, X, theta_hat, lambda_hat, n),
       "LS" = LS_test(Y_v_dt, n))
}

bootstrap_sample_test_value <- function(n, X, C, theta_hat, lambda_hat){
  X_bstr <- rweibull(n, theta_hat, lambda_hat)
  
  #C_bstr <- non_param_C_bootstrap(n, C, list(X == C), theta_hat, lambda_hat)
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
       "CM" = CM_test(T_bstr_v_dt, X_bstr, theta_hat_bstr, lambda_hat_bstr, n),
       "LS" = LS_test(Y_bstr_v_dt, n))
}


compute_statistics_cens <- function(n, perc_censored){
  X <- rweibull(n, 1.5, 1)
  m_thres <- pracma::bisect(m_find(X, perc_censored), 0, 10000)$root # znaleznie takiego m, by w próbie było pi=0.3 ocenzurowanych obs
  C <- runif(n, 0, m_thres)
  
  T_vec <- sapply(1:n, function(i) min(X[i], C[i]))
  delta <- T_vec == X
  T_v_dt <- data.table("T_val" = T_vec, "delta" = delta)
  
  theta_hat_0 <- mle_theta_full_sample(n, T_v_dt[, T_val])
  lambda_hat_0 <- mle_lambda_full_sample(n, T_v_dt[, T_val], theta_hat_0)
  lambda_hat_0
  
  params <- Newton_result(T_v_dt, theta_hat_0, lambda_hat_0) 
  theta_hat <- params[1]
  lambda_hat <- params[2]
  
  Y_v_dt <- transform_to_Y(T_v_dt, theta_hat, lambda_hat)
  
  list("X" = X, 
       "C" = C, 
       "theta_hat" = theta_hat,
       "lambda_hat" = lambda_hat,
       "Y_v_dt" = Y_v_dt)
}



Cens_statistics_compute_list <- lapply(c(50, 100), 
                                       function(n) lapply(c(0.1, 0.2), function(perc) compute_statistics_cens(n, perc)))
names(Cens_statistics_compute_list) <- c(50, 100)
names(Cens_statistics_compute_list[["50"]]) <- c("10%", "20%")
names(Cens_statistics_compute_list[["100"]]) <- c("10%", "20%")

saveRDS(Cens_statistics_compute_list, file = "/Users/maniek/Desktop/licencjat_R/Cens_statistics_compute_list.RDS")
m <- 5000

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

list_of_test_vals <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{
  
  data.table("uncensured_50" = compute_test_values_uncens(50),
             "uncensured_100" = compute_test_values_uncens(100),
             "censured_10perc_50" = bootstrap_sample_test_value(50, Cens_statistics_compute_list[["50"]][["10%"]]$X,
                                                                Cens_statistics_compute_list[["50"]][["10%"]]$C,
                                                                Cens_statistics_compute_list[["50"]][["10%"]]$theta_hat,
                                                                Cens_statistics_compute_list[["50"]][["10%"]]$lambda_hat),
             "censured_20perc_50" = bootstrap_sample_test_value(50, Cens_statistics_compute_list[["50"]][["20%"]]$X,
                                                                Cens_statistics_compute_list[["50"]][["20%"]]$C,
                                                                Cens_statistics_compute_list[["50"]][["20%"]]$theta_hat,
                                                                Cens_statistics_compute_list[["50"]][["20%"]]$lambda_hat),
             "censured_10perc_100" = bootstrap_sample_test_value(100, Cens_statistics_compute_list[["100"]][["10%"]]$X,
                                                                 Cens_statistics_compute_list[["100"]][["10%"]]$C,
                                                                 Cens_statistics_compute_list[["100"]][["10%"]]$theta_hat,
                                                                 Cens_statistics_compute_list[["100"]][["10%"]]$lambda_hat),
             "censured_20perc_100" = bootstrap_sample_test_value(100, Cens_statistics_compute_list[["100"]][["20%"]]$X,
                                                                 Cens_statistics_compute_list[["100"]][["20%"]]$C,
                                                                 Cens_statistics_compute_list[["100"]][["20%"]]$theta_hat,
                                                                 Cens_statistics_compute_list[["100"]][["20%"]]$lambda_hat))
}
stopImplicitCluster()

save(list_of_test_vals, file = "/Users/maniek/Desktop/licencjat_R/list_of_test_vals.rda")
load("/Users/zuza/Desktop/studia/licencjat/SimulationData/licencjat_R/list_of_test_vals.rda")


getStatisticTestsValue <- function(i, strColName){
  unlist(lapply(list_of_test_vals, function(dt) dt[i, get(strColName)]))
}

getDataTableForSample <- function(strColName){
  dt <- as.data.table(sapply(1:9, function(i) getStatisticTestsValue(i, strColName)))
  names(dt) <- c("S_n_1_1", "S_n_1_2", "S_n_1_5", "S_n_2_1", "S_n_2_2", "S_n_2_5", "KS", "CM", "LS")
  dt
}

getQuantileStatisticTest <- function(dt, strTestName, m = 10000, alpha = 0.05){
  test_vals <- dt[, get(strTestName)]
  sort(test_vals)[floor(m * (1 - alpha))]
}

getAllQuantiles <- function(dt, m, alpha = 0.05){
  sapply(names(dt), function(name) getQuantileStatisticTest(dt, name, m, alpha))
}

getDataFrameQuantiles <- function(m, alpha = 0.05){
  
  list_of_dt <- list(Uncensured_50_dt, Uncensured_100_dt, Censured_10perc_50_dt, 
                     Censured_20perc_50_dt, Censured_10perc_100_dt, Censured_20perc_100_dt)
  
  list_to_dt <- lapply(list_of_dt, function(dt) getAllQuantiles(dt, m, alpha))
  
  df <- as.data.frame(t(data.frame(list_to_dt[[1]], list_to_dt[[2]], list_to_dt[[3]], list_to_dt[[4]], list_to_dt[[5]], list_to_dt[[6]])))
  rownames(df) <- c(deparse(substitute(Uncensured_50)), deparse(substitute(Uncensured_100)),
                    deparse(substitute(Censured_10perc_50)),deparse(substitute(Censured_20perc_50)),
                    deparse(substitute(Censured_10perc_100)),deparse(substitute(Censured_20perc_100)))
  df
}

Uncensured_50_dt <- getDataTableForSample("uncensured_50")
Uncensured_100_dt <-  getDataTableForSample("uncensured_100")
Censured_10perc_50_dt <- getDataTableForSample("censured_10perc_50")
Censured_20perc_50_dt <- getDataTableForSample("censured_20perc_50")
Censured_10perc_100_dt <- getDataTableForSample("censured_10perc_100")
Censured_20perc_100_dt <- getDataTableForSample("censured_20perc_100")

Quantiles_DF <- getDataFrameQuantiles(5000)
