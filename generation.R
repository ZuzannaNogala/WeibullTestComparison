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
       "KS" = KS_test_full(X, theta_hat, lambda_hat, n),
       "CM" = CM_test_full(X, theta_hat, lambda_hat, n),
       "LS" = LS_test_full(X, n))
}

m <- 10000

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

list_of_values_to_quantiles <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{
  
  data.table("uncensured_10" = compute_test_values_uncens(10),
             "uncensured_50" = compute_test_values_uncens(50),
             "uncensured_100" = compute_test_values_uncens(100))
  
}
stopImplicitCluster()

# CREATE QUANTILES DATA.FRAME FROM GENERATED DATA 

getStatisticTestsValue <- function(list_of_vals, i, strColName){
   unlist(lapply(list_of_vals, function(dt) dt[i, get(strColName)]))
}

getDataTableForSample <- function(list_of_vals, strColName){
   dt <- as.data.table(sapply(1:9, function(i) getStatisticTestsValue(list_of_vals, i, strColName)))
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


getDataFrameQuantiles <- function(list_of_dt, m, alpha = 0.05){
  list_to_dt <- lapply(list_of_dt, function(dt) getAllQuantiles(dt, m, alpha))
  
  df <- t(do.call(cbind, lapply(list_to_dt, data.frame)))
  rownames(df) <- names(list_of_dt)
  df
}

list_of_dt <- list("Uncensured_10" = getDataTableForSample(list_of_values_to_quantiles, "uncensured_10"),
                   "Uncensured_50" = getDataTableForSample(list_of_values_to_quantiles, "uncensured_50"),
                   "Uncensured_100" = getDataTableForSample(list_of_values_to_quantiles, "uncensured_100"))

Quantiles_DF <- getDataFrameQuantiles(list_of_dt, m, alpha = 0.1)
