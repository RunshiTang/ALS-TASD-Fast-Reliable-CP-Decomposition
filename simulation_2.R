setwd("~/density_CP/04102025")

set.seed(0)
pp = 30
dimensions <- c(pp, pp, pp)
R_list = c(3)

#condition_number = c(1,1000,10000,100000)

alpha = seq(0.25,1, 0.05)



result_df = expand.grid(R = R_list, iter = 1:100, alpha = alpha)
result_df$rand_seed = sample.int(1000000, length(result_df$R), replace = TRUE)


library(doParallel)
myCluster <- makeCluster(64, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)
result_df <- foreach(i = 1:nrow(result_df), .combine = 'rbind') %dopar% {
  library(rTensor)
  source("helpers.R")
  
  set.seed(result_df$rand_seed[i])
  
  alpha = result_df$alpha[i]
  sigma = 15/pp^alpha
  
  R = result_df$R[i]
  random_tensor <- generate_random_cp_tensor(dimensions, R)
  
  noise_tensor = rand_tensor(modes = as.tensor(random_tensor[["tensor"]])@modes)
  
  # result_random_initialized_cp = cp(as.tensor(random_tensor[["tensor"]] + noise_tensor@data * sigma), R)
  # 
  # result_our_initialization = h_d(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  result_our_full_approach = our_method_new(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  # result_sim_diag = sim_diag(random_tensor$tensor + noise_tensor@data * sigma, R)
  
  rlist = c(R = R, 
            alpha = alpha,
            sigma= sigma,
            #random_initialized_cp_loss = my_loss(random_tensor$U, result_random_initialized_cp$U),
            #our_initialization_loss = my_loss(random_tensor$U, result_our_initialization$U),
            #sim_diag_loss = my_loss(random_tensor$U, result_sim_diag$U),
            our_full_approach_loss = my_loss(random_tensor$U, result_our_full_approach$cp_result$U)
  )
  
  return(rlist)
  
}
stopCluster(myCluster)

result_df_save = result_df
write.csv(result_df_save, "result_df_save_2.csv", row.names = F)

library(tidyverse)
#result_df = result_df_save
result_df = as.data.frame(result_df) %>% group_by(alpha) %>%
  mutate(loss_mean = mean(our_full_approach_loss),
            loss = our_full_approach_loss)

result_df_2 = as.data.frame(result_df) %>% group_by(alpha) %>%
  summarise(loss_median = median(our_full_approach_loss),
            loss_sd = sd(our_full_approach_loss))


p = ggplot() +
  geom_line(data =  result_df_2, aes(x = alpha, y = loss_median)) +
  geom_point(data =  result_df, aes(x = alpha, y = loss),
              alpha = 0.3) +
  labs(
    x = "Alpha",
    y = "loss"
  ) + 
  scale_y_log10()+
  geom_vline(xintercept = 0.75, linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 0.75, y = min(result_df_2$loss_median), label = "alpha = 0.75", 
           vjust = 2, hjust = 0.5)+
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 0.5, y = min(result_df_2$loss_median), label = "alpha = 0.5", 
           vjust = 2, hjust = 0.5)

ggsave("plot_2.pdf", 
       p, 
       width = 7, height = 3)
