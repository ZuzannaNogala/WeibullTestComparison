# POWER COMPUTING FOR FULL SAMPLE, DIFFRENT DISTRIBUTIONS
# CHECK IF RESULTS AS SIMMILAR AS IN ORIGINAL WORK

library(data.table)
library(foreach)
library(doParallel)
# library(xtable)


m <- 10000

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

list_of_df_amount_of_H1_wins <- foreach(i = 1:m, .packages = c("foreach", "data.table")) %dopar%{

  W_1 <- rweibull(100, 0.5, 1)
  W_2 <- rweibull(100, 1, 1)
  W_3 <- rweibull(100, 1.5, 1)
  
  G_1 <- rgamma(100, 2, 1)
  G_2 <- rgamma(100, 3, 1)
  
  LN_1 <- rlnorm(100, 0, 0.5)
  LN_2 <- rlnorm(100, 0, 1)
  
  Chi_1 <- rchisq(100, 8)
  Chi_2 <- rchisq(100, 10)
  
  B_2 <- rbeta(100, shape1 = 1, shape2 = 1)
  B_1 <- rbeta(100, shape1 = 0.5, shape2 = 1)
  
  
  data.frame("W_1" = sampling_uncensured(W_1, "Uncensured_100"), # H0 ok
             "W_2" = sampling_uncensured(W_2, "Uncensured_100"),
             "W_3" = sampling_uncensured(W_3, "Uncensured_100"),
             "G_1" = sampling_uncensured(G_1, "Uncensured_100"), # H1 ok
             "G_2" = sampling_uncensured(G_2, "Uncensured_100"),
             "LN_1" = sampling_uncensured(LN_1, "Uncensured_100"), 
             "LN_2" = sampling_uncensured(LN_2, "Uncensured_100"),
             "Chi_1" = sampling_uncensured(Chi_1, "Uncensured_100"), 
             "Chi_2" = sampling_uncensured(Chi_2, "Uncensured_100"),
             "B_1" = sampling_uncensured(B_1, "Uncensured_100"), 
             "B_2" = sampling_uncensured(B_2, "Uncensured_100"))
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
                     lapply(colnames(Quantiles_DF), 
                            function(rowName){
                              computePower(list_of_df_amount_of_H1_wins,
                                           distName = colName, 
                                           StatTestName = rowName)})})

rownames(Power_df) <- colnames(Quantiles_DF)

Power_df<- t(Power_df)

print(xtable(Power_df, type = "latex"), file = "...")

#load("/Users/zuza/Downloads/bibibibu/Sample_vals_another_2.rda")