setwd("~/density_CP/04102025")

set.seed(0)

threshold = 0.05

d = c(3, 5)

pp = 5

R_list = c(3)
sigma = c(0.01)
cor_list = seq(0,0.9, 0.05)



result_df = expand.grid(d = d, R = R_list, iter = 1:300, sigma = sigma, cor = cor_list)
result_df$rand_seed = sample.int(1000000, length(result_df$R), replace = TRUE)


library(doParallel)
myCluster <- makeCluster(64, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)
result_df <- foreach(i = 1:nrow(result_df), .combine = 'rbind') %dopar% {
  library(rTensor)
  source("helpers.R")
  
  set.seed(result_df$rand_seed[i])
  
  cor = result_df$cor[i]
  sigma = result_df$sigma[i]
  d = result_df$d[i]
  R = result_df$R[i]
  
  dimensions <- rep(pp, d)
  
  random_tensor <- generate_random_cp_tensor_high_d_with_cor(dimensions, R, cor = cor)
  
  noise_tensor = rand_tensor(modes = as.tensor(random_tensor[["tensor"]])@modes)
  
  result_our_full_approach = our_method_new(random_tensor$tensor + noise_tensor@data * sigma, R, tol = 0, max_iter = 20,
                                            test_conv = T)
  
  loss_per_iter = c()
  for (j in 1:length(result_our_full_approach$U_all_list)) {
    loss_per_iter = c(loss_per_iter, my_loss(result_our_full_approach$U_all_list[[j]], random_tensor$U))
  }
  
  num_iter = min(which(loss_per_iter<threshold))
  
  rlist = c(R = R, 
            d = d,
            sigma =sigma,
            num_iter = num_iter,
            cor_coef = cor
  )
  
  return(rlist)
  
}
stopCluster(myCluster)
result_df = as.data.frame(result_df)
result_df_final = result_df



write.csv(result_df, "result_df_save_6_1.csv", row.names = F)




# result_df = result_df_final %>% select(-R)%>% filter(num_iter < Inf) %>%
#   group_by(d, cor_coef, sigma) %>%
#   summarise(iter_mean = mean(num_iter)) 
# 
# p1 = ggplot(result_df) + 
#   geom_point(aes(x = as.factor(cor_coef), y = iter_mean, colour = as.factor(sigma))) + 
#   xlab("Xi")+
#   ylab("Proportion")+
#   labs(fill = "Number of iterations")+
#   ggtitle("d = 3")+ 
#   guides(fill = guide_legend(nrow = 1))+
#   scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)])


library(patchwork)

result_df = result_df_final %>% select(-R)%>% filter(sigma == 0.01) %>%
  group_by(d, cor_coef, sigma) %>%
  summarise(`0 iteration` = mean(num_iter == 1),
            `1 iteration` = mean(num_iter == 2),
            `2 iteration` = mean(num_iter == 3),
            `3 iteration` = mean(num_iter == 4),
            `4+ iteration` = mean(num_iter > 4 & num_iter < Inf),
            `never` = mean(num_iter == Inf)) %>%
  pivot_longer(4:9, names_to = "num_iter", values_to = "proportion")


p1 = ggplot(result_df %>% filter(d==3)) + 
  geom_bar(aes(x = as.factor(cor_coef), y = proportion, fill = as.factor(num_iter)),
           stat = "identity") + 
  xlab("Xi")+
  ylab("Proportion")+
  labs(fill = "")+
  ggtitle("d = 3")+ 
  guides(fill = guide_legend(nrow = 1))+
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)])

p2 = ggplot(result_df %>% filter(d==5)) + 
  geom_bar(aes(x = as.factor(cor_coef), y = proportion, fill = as.factor(num_iter)),
           stat = "identity") + 
  xlab("Xi")+
  ylab("Proportion")+
  labs(fill = "")+
  ggtitle("d = 5")+ 
  guides(fill = guide_legend(nrow = 1))+
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)])

# p3 = ggplot(result_df %>% filter(d==7)) + 
#   geom_bar(aes(x = as.factor(cor_coef), y = proportion, fill = as.factor(num_iter)),
#            stat = "identity") + 
#   xlab("Xi")+
#   ylab("Proportion")+
#   labs(fill = "")+
#   ggtitle("d = 7")+ 
#   guides(fill = guide_legend(nrow = 1))+
#   scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)])




ggsave("plot_6_1.pdf", 
       ggpubr::ggarrange(p1,p2, ncol = 2, nrow = 1, common.legend = T)+
         scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]), 
       width = 9, height = 3)



