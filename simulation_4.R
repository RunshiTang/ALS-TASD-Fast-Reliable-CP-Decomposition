setwd("~/density_CP/04102025")

set.seed(123)

dimensions <- c(10, 12, 15)
R_list = c(2)

#condition_number = c(1,1000,10000,100000)

result_df = expand.grid(R = R_list, iter = 1:50, sigma = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1))
result_df$rand_seed = sample.int(1000000, length(result_df$R), replace = TRUE)


library(doParallel)
myCluster <- makeCluster(64, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)
result_df <- foreach(i = 1:nrow(result_df), .combine = 'rbind') %dopar% {
  library(rTensor)
  source("helpers.R")
  
  set.seed(result_df$rand_seed[i])
  
  sigma = result_df$sigma[i]
  R = result_df$R[i]
  random_tensor <- generate_random_cp_tensor_condition_number_1(dimensions, R)
  
  noise_tensor = rand_tensor(modes = as.tensor(random_tensor[["tensor"]])@modes)
  
  result_random_initialized_cp = cp(as.tensor(random_tensor[["tensor"]] + noise_tensor@data * sigma), R)
  
  result_our_initialization = h_d_new(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  result_our_full_approach = our_method_new(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  result_sim_diag = sim_diag(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  result_r1als = r1als(random_tensor$tensor + noise_tensor@data * sigma, R, nu = 0.1)
  
  result_r1als_2 = r1als_2(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  result_ico = ico(random_tensor$tensor + noise_tensor@data * sigma, 
                   R,
                   a0 = cpca(random_tensor$tensor + noise_tensor@data * sigma, 
                             R = R)$U)
  
  if (any(is.na(result_r1als$U))){
    r1als_loss = NA
  }else{
    r1als_loss = my_loss(random_tensor$U, result_r1als$U)
  }
  
  rlist = c(R = R, 
            sigma = sigma,
            random_initialized_cp_loss = my_loss(random_tensor$U, result_random_initialized_cp$U),
            our_full_approach_loss = my_loss(random_tensor$U, result_our_full_approach$cp_result$U),
            our_initialization_loss = my_loss(random_tensor$U, result_our_initialization$U),
            sim_diag_loss = my_loss(random_tensor$U, result_sim_diag$U),
            r1als_loss = r1als_loss,
            r1als_loss_2 = my_loss(random_tensor$U, result_r1als_2$U),
            ico_loss = my_loss(random_tensor$U, result_ico$U)
  )
  
  return(rlist)
  
}
stopCluster(myCluster)

result_df_save = result_df
write.csv(result_df_save, "result_df_save2.csv", row.names = F)

library(tidyverse)

result_df = as.data.frame(result_df) %>% pivot_longer(cols = 3:9, names_to = "type", values_to = "loss")

# result_df = result_df %>% group_by(R, sigma, type) %>% summarise(loss_mean = mean(loss),
#                                                                  loss_sd = sd(loss))




result_df$sigma = as.factor(result_df$sigma)

# R_list = unique(result_df$R)

result_df$type <- factor(result_df$type,
                        levels = c("our_full_approach_loss",
                                   "ico_loss", 
                                   "random_initialized_cp_loss", 
                                   "r1als_loss",
                                   "r1als_loss_2",
                                   "our_initialization_loss",
                                   "sim_diag_loss"),
                        labels = c("TASD-ALS", "CPCA-ICO", "Random ALS", "R1-ALS-1", "R1-ALS-2", 
                                   "TASD", "SimDiag")
)



p = ggplot(result_df %>% filter(type %in% c("TASD-ALS", "CPCA-ICO", 
                                            "Random ALS", "R1-ALS-1", "R1-ALS-2", 
                                            "TASD", "SimDiag"))) + 
  geom_boxplot(aes(x = sigma, y = loss, fill = type), outlier.size = 1,
               outlier.alpha = 1, lwd = 0.2) + 
  xlab("Sigma")+
  labs(fill = "Method")

ggsave("plot_4.pdf", 
       p, 
       width = 9, height = 4)
