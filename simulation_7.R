
# compare new approach with exhaustive search

setwd("~/density_CP/04102025")

set.seed(0)


d = 3

pp = 15
dimensions <- rep(pp, d)
R_list = c(2,3,4, 5, 6)
sigma = c(0.1,0.5, 1)


result_df = expand.grid(R = R_list, iter = 1:50, sigma = sigma)
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
  random_tensor <- generate_random_cp_tensor_high_d(dimensions, R)
  
  noise_tensor = rand_tensor(modes = as.tensor(random_tensor[["tensor"]])@modes)
  
  time1 = Sys.time()
  result_our_initialization = h_d(random_tensor$tensor + noise_tensor@data * sigma, R)
  result_our_initialization_time_cost = Sys.time() - time1
  
  time1 = Sys.time()
  result_our_initialization_new = h_d_new(random_tensor$tensor + noise_tensor@data * sigma, R)
  result_our_initialization_new_time_cost = Sys.time() - time1
  
  obs_norm = fnorm(as.tensor(random_tensor$tensor + noise_tensor@data * sigma))
  
  # loss_per_iteration = rep(NA, 10)
  # loss_per_iteration[1] = result_our_full_approach$initialization_loss_X#/obs_norm
  # loss_per_iteration[2:(1+length(result_our_full_approach$cp_result$all_resids))] = result_our_full_approach$cp_result$all_resids#/obs_norm
  # 
  
  rlist = c(R = R, 
            #sigma= sigma,
            sigma =sigma,
            our_initialization_loss = my_loss(random_tensor$U, result_our_initialization$U),
            our_initialization_new_loss = my_loss(random_tensor$U, result_our_initialization_new$U),
            result_our_initialization_time_cost = result_our_initialization_time_cost,
            result_our_initialization_new_time_cost = result_our_initialization_new_time_cost
  )
  
  return(rlist)
  
}
stopCluster(myCluster)
result_df = as.data.frame(result_df)
result_df$d = d
result_df_final = result_df



library(patchwork)

#threshold = 0.1

# result_df = result_df %>% group_by(sigma, d) %>%
#   mutate(congverge_percent = sum(final_loss < threshold)/n())



result_df = as.data.frame(result_df)[,1:4] %>% pivot_longer(cols = 3:4, names_to = "type", values_to = "loss")


result_df$type <- factor(result_df$type,
                         levels = c("our_initialization_loss", 
                                    "our_initialization_new_loss"),
                         labels = c("Exhaustively Search", "TASD-Part II")
)

plot_list = list()
for (sigma_tmp in unique(result_df$sigma)) {
  p1 = ggplot(result_df %>% filter(sigma == sigma_tmp)) + 
    geom_boxplot(aes(x = as.factor(R), y = loss, fill = type)) + 
    ggtitle(paste0("sigma = ", sigma_tmp)) + 
    xlab("R")+
    labs(fill = "Method")+
    scale_y_log10()
  plot_list[[length(plot_list) + 1]] = p1
}

ggsave("plot_7.pdf", 
       ggpubr::ggarrange(plotlist = plot_list, ncol = 3, nrow = 1, common.legend = T), 
       width = 9, height = 4)


















