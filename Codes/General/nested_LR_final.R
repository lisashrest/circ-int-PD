library(dplyr)
library(plyr)
library(tidyverse)
library(pROC)
library(glmnet)
library(caret)
library(nestedcv)

setwd("data/")

#write the function for updrs 3 OFF and ON progression
updrs_3_progression_model <- function(mat, len){
  #Making model to test the progress of updrs 3 scores
  #Take the required columns only
  set.seed(1234)
  updrs3_final_mat <- mat[!is.na(mat$change_updrs3), ]
  final_mat <- updrs3_final_mat[, c(1:len, ncol(updrs3_final_mat)-2)]
  
  #Extract predictor variables
  X_train <- model.matrix(Progression_updrs3 ~. , data = final_mat)[, -ncol(final_mat)]
  
  #Extract the response variables
  y_train <- as.factor(final_mat$Progression_updrs3)
  
  #Random Sampling
  out <- randomsample(y_train, X_train)
  y_rs <- out$y
  X_rs <- out$x
  
  #Check no of observations
  print(dim(X_rs))
  
  #Fit regularized multivariate logistic regression model using cross-validation
  cv_fit <- nestedcv::nestcv.glmnet(x = X_rs, y = y_rs,
                                    n_outer_folds = 10, n_inner_folds = 10,                         
                                    outer_train_predict = TRUE,
                                    family = "binomial", alphaSet = seq(0, 1, 0.05), 
                                    finalCV = TRUE, cv.cores = 4)
  
  return(list(final_model = cv_fit))
  
}


updrs_3_progression_on_model <- function(mat, len){
  
  #Making model to test the progress of updrs 3 scores
  #Take the required columns only
  set.seed(1234)
  updrs3_final_mat <- mat[!is.na(mat$change_updrs3_on), ]
  final_mat <- updrs3_final_mat[, c(1:len, ncol(updrs3_final_mat)-1)]
  
  #Extract predictor variables
  X_train <- model.matrix(Progression_updrs3_on ~. , data = final_mat)[, -ncol(final_mat)]
  
  #Extract the response variables
  y_train <- as.factor(final_mat$Progression_updrs3_on)
  
  #Random Sampling
  out <- randomsample(y_train, X_train)
  y_rs <- out$y
  X_rs <- out$x
  
  #Check no of observations
  print(dim(X_rs))
  
  #Fit regularized multivariate logistic regression model using cross-validation
  cv_fit <- nestedcv::nestcv.glmnet(x = X_rs, y = y_rs,
                                    n_outer_folds = 10, n_inner_folds = 10, 
                                    outer_train_predict = TRUE,
                                    family = "binomial", alphaSet = seq(0, 1, 0.05), 
                                    finalCV = TRUE, cv.cores = 4)
  
  return(list(final_model = cv_fit))
  
}

#import clinical metadata from PPMI
details <- read.csv("meta/sample_details_1.csv", stringsAsFactors = FALSE, header = TRUE)
details$filename <- paste0("P1_", details$PATNO, "_", details$EVENT_ID)
rownames(details) <- details$filename

#Creating Object List 
details_up <- subset(details, details$genetic == "LRRK2-/SNCA-/GBA-")
details_up$filename <- paste0("P1_", details_up$PATNO, "_", details_up$EVENT_ID)

#Import the list of genes of interest
circ_imp <- read.csv("meta/circ_of_interest.csv", header = TRUE, row.names = 1)

#loading normalized reads 
norm_cts <- circ_imp[, details_up$filename]
head(norm_cts)

#Obtain the read counts for the host genes 
hs_cts <- read.csv("meta/genes_of_interest.csv", header =TRUE, row.names = 1)
hs_cts <- hs_cts %>%
  rownames_to_column("ID")%>%
  separate(col = ID, into = c("ENSEMBL", "Gene"), sep = "\\|")

rownames(hs_cts) <- hs_cts$Gene

#For OFF state
#import the updrs scores data
updrs3 <- read.csv("meta/change_scores.csv", header = TRUE)

updrs3$Y2 <- "V06"
updrs3$filename <- paste0("P1_", updrs3$PATNO, "_", updrs3$Y2)
head(updrs3)

#extract IDs for samples from PD participants at Year 2 and Baseline
meta_pd <- details_up %>%
  dplyr::filter(APPRDX == 1) %>%
  dplyr::filter(EVENT_ID == "BL" | EVENT_ID == "V06")

head(meta_pd)

#Separate metadata for samples collected at Baseline
meta_pd_bl <- meta_pd %>%
  dplyr::filter(EVENT_ID == "BL")

#Separate metadata for samples collected at Year 2 
meta_pd_y2 <- meta_pd %>%
  dplyr::filter(EVENT_ID == "V06")

#Get normalized values for samples at Baseline
normalized_bl <- norm_cts[, colnames(norm_cts)%in%meta_pd_bl$filename]
head(normalized_bl)

#Get normalized values for samples at Year 2
normalized_y2 <- norm_cts[, colnames(norm_cts)%in%meta_pd_y2$filename]

#Get change in normalized expression
normalized_final <- normalized_y2 - normalized_bl

#Get normalized values for samples at Baseline
hs_bl <- hs_cts[, colnames(hs_cts)%in%meta_pd_bl$filename]
head(hs_bl)
hs_y2 <- hs_cts[, colnames(hs_cts)%in%meta_pd_y2$filename]

#Get change in normalized expression
hs_final <- hs_y2 - hs_bl

host_gene_final <- hs_final %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("filename")

circ_final <- normalized_final %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("filename") 
circ_final <- circ_final[, c(2:7,1)]

head(circ_final)

#--------------------------------------BSJ Only ----------------------------------

final_mat <- inner_join(circ_final, updrs3, by = "filename")
final_mat$Progression_updrs3 <- factor(final_mat$Progression_updrs3,
                                       levels = c("Unchanged", "Progressed"))

final_mat$Progression_updrs3_on <- factor(final_mat$Progression_updrs3_on, 
                                          levels = c("Unchanged", "Progressed"))

#Add the number of circRNAs in the model
len <- 6

#updrs3
updrs3_BSJ <- updrs_3_progression_model(final_mat, len)

#updrs3_on
updrs3_BSJ_on <- updrs_3_progression_on_model(final_mat, len)

#----------------------------------------- BSJ + Gene ----------------------------------
col <- "filename"

final_matrix <- cbind(circ_final[, !(names(circ_final)%in%col)], host_gene_final[, !(names(host_gene_final)%in%col)])
final_matrix$filename <- circ_final$filename
head(final_matrix)

f_matrix <- inner_join(final_matrix, updrs3, by = "filename")
f_matrix$Progression_updrs3 <- factor(f_matrix$Progression_updrs3,
                                      levels = c("Unchanged", "Progressed"))

f_matrix$Progression_updrs3_on <- factor(f_matrix$Progression_updrs3_on, 
                                         levels = c("Unchanged", "Progressed"))

#Add the number of circRNAs in the model
len <- 17

#updrs3
updrs3_BSJ_Gene <- updrs_3_progression_model(f_matrix, len)

#updrs3_on
updrs3_BSJ_Gene_on <- updrs_3_progression_on_model(f_matrix, len)

#----------------------------------------- Host Gene Only  ----------------------------------
final_matrix_g <- host_gene_final[, -1]
final_matrix_g$filename <- host_gene_final$filename

final_matrix_g <- inner_join(final_matrix_g, updrs3, by = "filename")
final_matrix_g$Progression_updrs3 <- factor(final_matrix_g$Progression_updrs3,
                                            levels = c("Unchanged", "Progressed"))

final_matrix_g$Progression_updrs3_on <- factor(final_matrix_g$Progression_updrs3_on, 
                                               levels = c("Unchanged", "Progressed"))

#Add the number of circRNAs in the model
len <- 11

#updrs3
updrs3_Gene <- updrs_3_progression_model(final_matrix_g, len)

#updrs3_on
updrs3_Gene_on <- updrs_3_progression_on_model(final_matrix_g, len)

#Export model redults 
models <- list(circ_u3 = updrs3_BSJ,
               circG_u3 = updrs3_BSJ_Gene,
               Gene_u3 = updrs3_Gene,
               circ_u3_on = updrs3_BSJ_on,
               circG_u3_on = updrs3_BSJ_Gene_on,
               Gene_u3_on = updrs3_Gene_on)

setwd("circRNA-sponging_pipeline")
saveRDS(models, file = "Models_final.rds")
