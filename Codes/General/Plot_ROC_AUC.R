
#BSJ_updrs_3 model
calculate_auc_roc <- function(model_list, metric_name){
  train_roc <- list()
  train_auc_final <- NULL 
  
  for (i in 1:length(model_list)) {
    train_roc[[i]] <- model_list[[i]]$final_model$roc
    train_auc <- ci.auc(model_list[[i]]$final_model$roc)
    train_auc <- data.frame(higher_ci = train_auc[3], auc = train_auc[2], lower_ci = train_auc[1])
    train_auc <- train_auc %>%
      mutate(metric = metric_name, name = names(model_list[i]))
    
    train_auc_final <- rbind(train_auc_final, train_auc)
  }
 return(list( Train_AUC = train_auc_final, Train_ROC = train_roc))
}

#Separating lists for OFF and ON scores
models_UPDRS_3 <- list(models$circ_u3, models$circG_u3, models$Gene_u3)
names(models_UPDRS_3) <- c("CircRNAs", "CircRNAs + Target Genes", "Target Genes")  

models_UPDRS_3_on <- list(models$circ_u3_on, models$circG_u3_on, models$Gene_u3_on)
names(models_UPDRS_3_on) <- c("CircRNAs", "CircRNAs + Target Genes", "Target Genes")  

#get output 
updrs3_res <- calculate_auc_roc(models_UPDRS_3, "UPDRS_III")
updrs3_on_res <- calculate_auc_roc(models_UPDRS_3_on, "UPDRS_III_(ON)")


#Plot ROC curve 
roc_list_off_train <- updrs3_res$Train_ROC
names(roc_list_off_train) <- c("circRNAs", "CircRNAs + Target Genes", "Target Genes")

updrs3_off <- ggroc(roc_list_off_train) +
  theme(legend.title = element_text(size = 12))+
  guides(col=guide_legend(title="Predictors"))+
  geom_abline(intercept = 1, slope = 1, lty =2, color = "gray")+
  scale_color_discrete(labels = names(roc_list_off_train))+
  ggtitle("UPDRS III (OFF) Model")

#Plot ROC curve 
roc_list_on_train <- updrs3_on_res$Train_ROC
names(roc_list_on_train) <- c("circRNAs", "CircRNAs + Target Genes", "Target Genes")

updrs3_on <- ggroc(roc_list_on_train) +
  theme(legend.title = element_text(size = 12))+
  guides(col=guide_legend(title="Predictors"))+
  geom_abline(intercept = 1, slope = 1, lty =2, color = "gray")+
  scale_color_discrete(labels = names(roc_list_on_train))+
  ggtitle("UPDRS III (ON) Model")

#AUC graph for updrs 3
train_auc_u3 <- bind_rows(updrs3_res$Train_AUC, updrs3_on_res$Train_AUC)

train_plot <- ggplot(train_auc_u3, aes(y = name, x = auc, colour = name)) +
  facet_wrap(~ metric, scales = "free_y") +
  geom_point(show.legend = FALSE) +
  geom_pointrange(aes(xmin = lower_ci, xmax = higher_ci)) +
  scale_colour_brewer(palette = "Set2") +
  labs(x = "AUC", y = "Predictor", colour = "Predictor") +
  theme(legend.position = 'bottom')+ 
  ggtitle("PPMI Cohort")

train_plot

#write.csv(train_auc_u3, file = "Logistic_regression/auc_final_results.csv")

cowplot::plot_grid(updrs3_off, updrs3_on,
                   labels = c("a", "b"),
                   ncol = 2,
                   label_size = 12)
