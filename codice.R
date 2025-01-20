

# -----------------------------------------
# Predictive Modeling and Network Analysis 
# of Gene Expression in Breast Cancer
# -----------------------------------------


# set the seed
set.seed(123)

# import the data
breast <- read.csv("Breast_GSE70947.csv", header=TRUE)
head(breast, n=c(6,4))

# remove the column with the samples names
breast <- breast[, -1]

# factorize the response
levels(as.factor(breast[,1]))
breast[,1] <- as.factor(breast[,1])

# recast the response variable type as 1 for tumor, 0 for control
breast[,1] <- ifelse(breast[,1] == "breast_adenocarcinoma", 1, 0)
dim(breast)
head(breast, n=c(6,5))

# chech if there are missing values
any(apply(breast, 2, is.na))
# FALSE


### Split the dataset into train and test

# get the indices
idx_train <- sample(1:nrow(breast), round(nrow(breast)*(3/4)))
idx_test <- setdiff(1:nrow(breast), idx_train)

y_train <- breast[idx_train, 1]
y_test <- breast[idx_test, 1]
breast <- breast[, -1]
# head(breast, n=c(6,4))

x_train <- as.matrix(breast[idx_train, ])
# head(x_train, n=c(6,4))

x_test <- as.matrix(breast[idx_test, ])
# head(x_test, n=c(6,4))


table(y_train)
# balanced


# free some memory
rm(breast)


#---------------------------------------
# Logistic Regression with LASSO penalty
# --------------------------------------

library(glmnet)


### Cross Validation for lambda

cv.fit <- cv.glmnet(x_train, y_train, family="binomial")

par(mfrow=c(1,2))
plot(cv.fit)
plot(cv.fit$glmnet.fit, xvar='lambda')
abline(v=log(cv.fit$lambda.min), lty=2)
abline(v=log(cv.fit$lambda.1se), lty=2)
par(mfrow=c(1,1))

cv.fit$lambda.min
# 0.01217137
cv.fit$lambda.1se
# 0.0746824


# Best:
log_lasso <- glmnet(x_train, y_train, lambda=cv.fit$lambda.min,
                    family='binomial')

# coefficients
coef.log_lasso <- as.matrix(coef(log_lasso))
coef.log_lasso <- subset(coef.log_lasso, coef.log_lasso != 0)  
length(coef.log_lasso)
# 72


# prediction error
pred.log_lasso <- as.numeric(predict(log_lasso, newx=x_test, type="class"))
err.log_lasso <- 1 - sum(diag(table(pred.log_lasso, y_test)))/NROW(y_test)
err.log_lasso
# 0.1111111

# compute other metrics
library(caret)

conf_matrix1 <- confusionMatrix(as.factor(pred.log_lasso), as.factor(y_test), 
                               positive = "1")

conf_matrix1$table
conf_matrix1$byClass["Sensitivity"]
# 0.8
conf_matrix1$byClass["Specificity"]  
# 1


#-------------------------------------------------
# SVM with LASSO penalization - squared hinge loss
# ------------------------------------------------

library(gcdnet)

set.seed(123)

### Cross validation
cv.svm1 <- cv.gcdnet(x_train, y_train, method="sqsvm", standardize=TRUE)

par(mfrow=c(1,2))
plot(cv.svm1)
plot(cv.svm1$gcdnet.fit, xvar='lambda')
abline(v=log(cv.svm1$lambda.min), lty=2)
abline(v=log(cv.svm1$lambda.1se), lty=2)
par(mfrow=c(1,1))

cv.svm1$lambda.min
# 0.0586
cv.svm1$lambda.1se
# 0.1709


# Best:
lasso_svm.sq <- gcdnet(x_train, y_train, method="sqsvm", 
                       lambda=cv.svm1$lambda.min, standardize=TRUE)

# coefficients
coef.lasso_svm.sq <- as.matrix(coef(lasso_svm.sq))
coef.lasso_svm.sq <- subset(coef.lasso_svm.sq, coef.lasso_svm.sq !=0)     
length(coef.lasso_svm.sq)
# 84

# prediction error
pred.lasso_svm.sq <- predict(lasso_svm.sq, x_test)	
# values: -1/1 => convert in 0/1
pred.lasso_svm.sq[pred.lasso_svm.sq == -1] <- 0

(err.lasso_svm.sq <- 1 - 
    sum(diag(table(pred.lasso_svm.sq, y_test)))/NROW(y_test))
# 0.111111


conf_matrix2 <- confusionMatrix(as.factor(pred.lasso_svm.sq), as.factor(y_test), 
                                positive = "1")
conf_matrix2$table
conf_matrix2$byClass["Sensitivity"] 
# 0.8
conf_matrix2$byClass["Specificity"]  
# 1


# free memory
detach(package:gcdnet)


#-----------------------------------------
# SVM with LASSO penalization - hinge loss
# ----------------------------------------

library(sparseSVM)

set.seed(123)

### Cross validation
cv.svm2 <- cv.sparseSVM(x_train, y_train, preprocess='standardize')

cv.svm2$lambda.min
# 0.1124

#1se rule
rule1se <- min(cv.svm2$lambda[cv.svm2$cve <= min(cv.svm2$cve) + 
                                cv.svm2$cvse[cv.svm2$min]])
rule1se
# 0.0356

par(mfrow=c(1,2))
plot(cv.svm2)
abline(v=log(cv.svm2$lambda.min), lty=2, lwd=1, col='black')
abline(v=log(rule1se), col='green', lty=2, lwd=1)
plot(cv.svm2$fit, xvar='lambda')
abline(v=log(cv.svm2$lambda.min), lty=2)
abline(v=log(rule1se), lty=2, col='green')
par(mfrow=c(1,1))

# coefficients
coef.lasso_svm.h <- coef(cv.svm2, lambda=cv.svm2$lambda.min)
coef.lasso_svm.h <- subset(coef.lasso_svm.h, coef.lasso_svm.h != 0) 
length(coef.lasso_svm.h)
# 117

# prediction error
pred.lasso_svm_h <- predict(cv.svm2, x_test, lambda=cv.svm2$lambda.min)	
(err.lasso_svm_h <- 1 - 
    sum(diag(table(pred.lasso_svm_h, y_test)))/NROW(y_test))
# 0.097222


conf_matrix3 <- confusionMatrix(as.factor(pred.lasso_svm_h), as.factor(y_test), 
                                positive = "1")
conf_matrix3$table
conf_matrix3$byClass["Sensitivity"] 
# 0.825
conf_matrix3$byClass["Specificity"]  
# 1


# free memory
detach(package:sparseSVM)



#----------------
# Relaxed LASSO
# --------------

# find the genes selected by the LASSO Logistic Regression:
names_selected_genes <- as.vector(rownames(coef.log_lasso))
matches <- which(colnames(x_train) %in% names_selected_genes)

relaxed_log_lasso <- glm(y_train ~ ., family='binomial',
                         data=as.data.frame(scale(x_train[, matches])))
# we have probalems, the algorithm does not converge
# because there is collinearity between explanatory variables


# the model does not work
summary(relaxed_log_lasso)


#### Check the correlations between the explanatory variables
library(corrplot)

cor_matrix <- cor(x_train[, matches])

# heatmap
corrplot(cor_matrix, method="color", tl.col=NA, tl.pos='n')


### Filter strong correlations

# threshold
thr <- 0.69

to_keep <- apply(abs(cor_matrix) > thr & abs(cor_matrix) < 1, 1, any)
filtered_cor_matrix <- cor_matrix[to_keep, to_keep]

corrplot(filtered_cor_matrix, tl.col=NA)
# some variables are strongly correlated to each others

# check the correlation of these variables also with the response
sapply(matches[to_keep], function(x) cor(x_train[,x], y_train))
# they are also strong correlated to the response



#-------------
# Elastic Net
# ------------

# generate fixed folds in order to make the results deterministic
set.seed(123)
nfolds <- 10
myfolds <- sample(rep(1:nfolds, length.out=nrow(x_train)))

# sequence for alpha
alpha <- seq(0, 1, length=6)

mods.elasticNet <- sapply(alpha, function(x) cv.glmnet(x_train, y_train, 
                                                       family="binomial",
                                                       alpha=x, 
                                                       foldid=myfolds,
                                                       nfolds=10))

best <- apply(mods.elasticNet, 2, function(x) min(x$cvm))
(alphaBest <- alpha[which.min(best)])
# 0.8
# so we search again for the best alpha in the neighbor of 0.8

alpha <- c(seq(alphaBest-0.1, alphaBest+0.1, length=10), alphaBest)
mods.elasticNet <- sapply(alpha, function(x) cv.glmnet(x_train, y_train,
                                                       family = "binomial",
                                                       alpha=x,
                                                       foldid=myfolds,
                                                       nfolds=10))

best <- apply(mods.elasticNet, 2, function(x) min(x$cvm))
(alphaBest <- alpha[which.min(best)])
# 0.7444444

# Cross Validation for lambda using the best alpha
cv.elasticNet <- cv.glmnet(x_train, y_train, family="binomial",
                           alpha=alphaBest)

par(mfrow=c(1,2))
plot(cv.elasticNet)
plot(cv.elasticNet$glmnet.fit, xvar='lambda')
abline(v=log(cv.elasticNet$lambda.min), lty=2)
abline(v=log(cv.elasticNet$lambda.1se), lty=2)
par(mfrow=c(1,1))

cv.elasticNet$lambda.min
# 0.0112
cv.elasticNet$lambda.1se
# 0.063


# Best Elastic Net:
elasticNet <- glmnet(x_train, y_train, lambda=cv.elasticNet$lambda.min,
                     alpha=alphaBest, family='binomial')

# coefficients
coef.elasticNet <- as.matrix(coef(elasticNet))
coef.elasticNet <- subset(coef.elasticNet, coef.elasticNet != 0)  
length(coef.elasticNet)
# 104

# prediction error
pred.elasticNet <- as.numeric(predict(elasticNet, newx=x_test, type="class"))
err.elasticNet <- 1 - sum(diag(table(pred.elasticNet, y_test)))/NROW(y_test)
err.elasticNet
# 0.09722

conf_matrix4 <- confusionMatrix(as.factor(pred.elasticNet), as.factor(y_test), 
                                positive = "1")
conf_matrix4$table
conf_matrix4$byClass["Sensitivity"] 
# 0.825
conf_matrix4$byClass["Specificity"]  
# 1


# --------------
# Random Forest
# --------------

library(ranger)

# hyperparameter grid search
hyper_grid <- expand.grid(
  mtry=seq(60, 200, by=20),
  ntrees=seq(200, 800, by=100),
  node_size=seq(1, 13, by=3),
  OOB_error=0)


# the following for cycle takes around 10 minutes...  
for(i in 1:nrow(hyper_grid)) {
  cat('performing evaluation', i, '/', nrow(hyper_grid), '\n')
  rf_temp <- ranger(x=x_train, y=as.factor(y_train), 
                    num.trees=hyper_grid$ntrees[i],
                    mtry=hyper_grid$mtry[i],
                    min.node.size=hyper_grid$node_size[i],
                    seed=123)
  
  hyper_grid$OOB_error[i] <- rf_temp$prediction.error
}


# results:
head(sort_by(hyper_grid, hyper_grid$OOB_error), 8)

(best_hyperparam <- sort_by(hyper_grid, hyper_grid$OOB_error)[1,])
# best: mtry=180, ntrees=400, node_size=10 => OOB_error=0.1059908

# refit the model
rf <- ranger(x=x_train, y=as.factor(y_train), 
             num.trees=best_hyperparam$ntrees,
             mtry=best_hyperparam$mtry,
             min.node.size=best_hyperparam$node_size,
             importance='impurity',
             seed=123)


# The importance of each variable is measured by the mean decrease 
# in Gini node impurity when the variable is used to make a split in a tree.
# => If a variable is significant for the model, splitting on it will 
# correspond to a greater reduction in node impurity.
# => The higher the decrease, the more important the variable is.
rf_imp_variables <- rf$variable.importance
rf_imp_variables <- subset(rf_imp_variables, rf_imp_variables>0)
length(rf_imp_variables)
# 3317

# variable importance plot
library(tidyverse)
rf$variable.importance %>% 
  broom::tidy() %>%
  dplyr::arrange(desc(x)) %>%
  dplyr::top_n(25) %>%
  ggplot(aes(reorder(names, x), x)) +
  geom_col() +
  coord_flip() +
  ggtitle("Top 25 important variables") +
  xlab('') + ylab('Gini Impurity decrease')


# predictions
pred.rf <- predict(rf, data=x_test)$predictions

(err.rf <- 1 - 
    sum(diag(table(pred.rf, y_test)))/NROW(y_test))
# 0.08333333

conf_matrix5 <- confusionMatrix(as.factor(pred.rf), as.factor(y_test), 
                                positive='1')
conf_matrix5$table
conf_matrix5$byClass["Sensitivity"] 
# 0.875
conf_matrix5$byClass["Specificity"]  
# 0.96875



# -------------------------------------------------------
# Consider the 20 most important variables for each model
# -------------------------------------------------------

# LASSO logistic regression
tibble(variable=rownames(coef.log_lasso), coef=coef.log_lasso) %>%  
  filter(abs(coef) > 0.35) %>% 
  arrange(desc(abs(coef))) %>% 
  print(n=20)

### NOTE: one gene can be in more than one category or in none of them!
gene_function_log_lasso <- c(6, 3, 7, 3, 2, 2, 1)
names(gene_function_log_lasso) <- c(
  "Cell Proliferation",
  "Cellular Stress Response",
  "Gene Expression Regulation",
  "Regulation of Cell Growth",
  "Cellular Metabolism",
  "Genomic Markers and Undefined Functions",
  "Unknown or Poorly Defined Biological Function"
)


# Elastic Net logistic regression
tibble(variable=rownames(coef.elasticNet), coef=coef.elasticNet) %>%  
  filter(abs(coef) > 0.35) %>% 
  arrange(desc(abs(coef))) %>%
  print(n=20)

gene_function_elnet <- c(5, 3, 8, 3, 2, 2, 2)
names(gene_function_elnet) <- c(
  "Cell Proliferation",
  "Cellular Stress Response",
  "Gene Expression Regulation",
  "Regulation of Cell Growth",
  "Cellular Metabolism",
  "Genomic Markers and Undefined Functions",
  "Unknown or Poorly Defined Biological Function"
)


# LASSO SVM square hinge loss
tibble(variable=rownames(coef.lasso_svm.sq), coef=coef.lasso_svm.sq) %>%  
  filter(abs(coef) > 0.10) %>% 
  arrange(desc(abs(coef))) %>% 
  print(n=20)

gene_function_SVMsq <- c(5, 0, 9, 0, 4, 0, 2)
names(gene_function_SVMsq) <- c(
  "Cell Proliferation",
  "Cellular Stress Response",
  "Gene Expression Regulation",
  "Regulation of Cell Growth",
  "Cellular Metabolism",
  "Genomic Markers and Undefined Functions",
  "Unknown or Poorly Defined Biological Function"
)


# LASSO SVM hinge loss
tibble(variable=names(coef.lasso_svm.h), coef=coef.lasso_svm.h) %>%  
  filter(abs(coef) > 0.05) %>% 
  arrange(desc(abs(coef))) %>%
  print(n=20)

gene_function_SVM_h <- c(3, 2, 7, 0, 4, 0, 2)
names(gene_function_SVM_h) <- c(
  "Cell Proliferation",
  "Cellular Stress Response,  stimula,",
  "Gene Expression Regulation",
  "Regulation of Cell Growth",
  "Cellular Metabolism",
  "Genomic Markers, and Undefined Functions",
  "Unknown or Poorly Defined Biological Function"
)


# Random Forest
tibble(variable=names(rf_imp_variables), importance=rf_imp_variables) %>%  
  filter(abs(importance) > 0.35) %>% 
  arrange(desc(abs(importance))) %>%
  print(n=20)

gene_function_RF<- c(6, 2, 8, 0, 2, 0, 2)
names(gene_function_RF) <- c(
  "Cell Proliferation",
  "Cellular Stress Response,  stimula,",
  "Gene Expression Regulation",
  "Regulation of Cell Growth",
  "Cellular Metabolism",
  "Genomic Markers, and Undefined Functions",
  "Unknown or Poorly Defined Biological Function"
)


top20_genes_function <- as.matrix(rbind(gene_function_log_lasso,
                                        gene_function_elnet,
                                        gene_function_SVMsq,
                                        gene_function_SVM_h,
                                        gene_function_RF))
View(top20_genes_function)

x11()
par(mfrow=c(3,2))
col <- c("red", "yellow", "green", "blue", "skyblue", "orange", "black")

barplot(top20_genes_function[1,], col=col, main="Logistic Lasso", names=FALSE, ylim= c(0,10))
barplot(top20_genes_function[2,], col=col, main="Elastic Net", names=FALSE, ylim= c(0,10))
barplot(top20_genes_function[3,], col=col, main="Lasso SVM - sq hinge", names=FALSE, ylim= c(0,10))
barplot(top20_genes_function[4,], col=col, main="Lasso SVM - hinge", names=FALSE, ylim= c(0,10))
barplot(top20_genes_function[5,], col=col, main="Random Forest", names=FALSE, ylim= c(0,10))
plot(1, type = "n", xlab = "", ylab = "", axes=FALSE)
legend("center", legend = names(gene_function_RF), fill = col, cex=1, bty= "n")
par(mfrow=c(1,1))


# ------------------
# Models comparison
# -----------------

df <- matrix(data=round(c(err.log_lasso, err.elasticNet, err.lasso_svm.sq, 
                          err.lasso_svm_h, err.rf, conf_matrix1$byClass["Sensitivity"], 
                          conf_matrix2$byClass["Sensitivity"], conf_matrix3$byClass["Sensitivity"],
                        conf_matrix4$byClass["Sensitivity"], conf_matrix5$byClass["Sensitivity"], 
                        conf_matrix1$byClass["Specificity"], conf_matrix2$byClass["Specificity"],
                        conf_matrix3$byClass["Specificity"], conf_matrix4$byClass["Specificity"],
                        conf_matrix5$byClass["Specificity"]), 
                        3), ncol=3, nrow=5)

colnames(df) <- c('Test Error', 'Sensitivity', 'Specificity')
rownames(df) <- c('Lasso Logistic Reg.', 'Elastic Net Logistic Reg.', 'Penalized s.h. SVM',
                  'Penalized h. SVM', 'Random Forest')
df


#------------------------------------
# Graphical Lasso - Normal vs Cancer
# -----------------------------------

# Correlation Graphs

# We want to create two graphs representing the correlations between the variables 
# selected with LASSO logistic regression with cv.fit$lambda.min, 
# one in normal patients, the other in patients with cancer
library(igraph)

# First, we need to separate the training set in normal patients and cancer patients:
rows_normal <- which(y_train == 0)
rows_cancer <- which(y_train == 1)
x_train0 <- x_train[rows_normal, matches]
x_train1 <- x_train[rows_cancer, matches]
# Correlation matrices
cor_matrix_normal <- cor(x_train0)
cor_matrix_cancer <- cor(x_train1)

# Now we convert them in Gaussian weights
sigma <- 0.5  # Gaussian transformation parameter
gaussian_weights_normal <- exp(-(1 - cor_matrix_normal)^2 / (2 * sigma^2))
gaussian_weights_cancer <- exp(-(1 - cor_matrix_cancer)^2 / (2 * sigma^2))
# And we use them as adjacency matrices for the graphs
# Note: strong correlation translates into big weight

# Partial Graphs
# Finding patterns and differences in the full graphs is difficult due 
# to the completeness of the graphs and the high number of correlations.
# Thus, we set a threshold, so not to visualize the weakest correlations.
threshold <- 0.3
adj_matrix_normal_part <- gaussian_weights_normal
adj_matrix_normal_part[abs(cor_matrix_normal) < threshold | cor_matrix_normal==1] <- 0

adj_matrix_cancer_part <- gaussian_weights_cancer
adj_matrix_cancer_part[abs(cor_matrix_cancer) < threshold | cor_matrix_cancer==1] <- 0

graph_normal_part <- graph_from_adjacency_matrix(adj_matrix_normal_part, mode = "undirected", weighted = TRUE)
graph_cancer_part <- graph_from_adjacency_matrix(adj_matrix_cancer_part, mode = "undirected", weighted = TRUE)

# Plot
# Create a fixed layout
set.seed(123)
fixed_layout <- layout_with_fr(graph_normal_part) 
# 1. We choose a fixed node layout because I want to compare the graphs
# 2. We use the function layout_with_fr() to obtain a graph where the stronger is the
# correlation between the variables, the shorter the distance between the nodes

# Graph for Normal Patients
plot(
  graph_normal_part,
  layout = fixed_layout,
  vertex.label = NA,
  vertex.size = 7,
  main = "Graph for Normal Patients"
)

# Graph for Cancer Patients
plot(
  graph_cancer_part,
  layout = fixed_layout,             
  vertex.label = NA,
  vertex.size = 7,
  main = "Graph for Cancer Patients"
)


# Partial graphs with heatmap
# To better visualize the correlations, we create two graphs where we use directly 
# the correlation value as edge weight and we introduce a heatmap.

# Define a function to map correlations (-1 to 1) to colors
get_edge_colors_direct <- function(correlations) {
  colors <- colorRampPalette(c("blue", "white", "red"))(100) # Blue to white to red
  indices <- round(((correlations + 1) / 2) * 99) + 1 # Map -1 to 1 into 0 to 1
  return(colors[indices])
}

# Threshold correlations to filter weak connections
adj_matrix_normal_direct <- cor_matrix_normal
adj_matrix_normal_direct[abs(cor_matrix_normal) < threshold | cor_matrix_normal == 1] <- 0

adj_matrix_cancer_direct <- cor_matrix_cancer
adj_matrix_cancer_direct[abs(cor_matrix_cancer) < threshold | cor_matrix_cancer == 1] <- 0

# Create graphs from adjacency matrices
graph_normal_direct <- graph_from_adjacency_matrix(adj_matrix_normal_direct, mode = "undirected", weighted = TRUE)
graph_cancer_direct <- graph_from_adjacency_matrix(adj_matrix_cancer_direct, mode = "undirected", weighted = TRUE)

# Get edge correlations and corresponding colors
edge_cor_normal <- E(graph_normal_direct)$weight
edge_colors_normal_direct <- get_edge_colors_direct(edge_cor_normal)

edge_cor_cancer <- E(graph_cancer_direct)$weight
edge_colors_cancer_direct <- get_edge_colors_direct(edge_cor_cancer)

# Connected components in the normal and the cancer graph
connected_components_0 <- components(graph_normal_direct)
connected_components_1 <- components(graph_cancer_direct)
# Number of connected components
connected_components_0$no
# 8
connected_components_1$no
# 12
connected_components_0$csize
# A big one (64 nodes) and the rest of the nodes are disconnected one from another
connected_components_1$csize
# Again, a big one (57 nodes), two small ones (2/3 nodes) and the rest is disconnected


# Plot
# Graph for normal patients
plot(
  graph_normal_direct,
  layout = fixed_layout,            
  vertex.label = NA,
  vertex.size = 7,
  vertex.color = "gray",
  edge.width = abs(edge_cor_normal) * 5, # Adjust width based on correlation strength
  edge.color = edge_colors_normal_direct, # Set edge colors based on correlation
  main = "Graph for Normal Patients (Direct Correlation Colors)"
)

# Graph for cancer patients
plot(
  graph_cancer_direct,
  layout = fixed_layout,             
  vertex.label = NA,
  vertex.size = 7,
  vertex.color="gray",
  edge.width = abs(edge_cor_cancer) * 5, # Adjust width based on correlation strength
  edge.color = edge_colors_cancer_direct, # Set edge colors based on correlation
  main = "Graph for Cancer Patients (Direct Correlation Colors)"
)



# Graphical Lasso

# Instead of working directly with correlations, we can also perform 
# graphical lasso. 

library(glasso)

# Create a function that fits graphical lasso given the covariance matrix S and lambda
# and creates the corresponding graph:
fit_graphical_lasso <- function(S, lambda) {
  #fit graphical lasso
  fit <- glasso(S, lambda)
  #get adjacency matrix from the precision matrix W of the fitted model
  adj_matrix <- as.matrix(fit$wi != 0)
  diag(adj_matrix) <- 0
  #create the graph
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  return(graph)
}

# Define S and lambda
S_normal <- cov(x_train0)
S_cancer <- cov(x_train1)
lambda_values <- c(0,0.1, 0.2, 0.3, 0.4, 0.5)

# Create a function to plot graphs for multiple lambda values
plot_graphical_lasso <- function(S, lambda_values, group_title) {
  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  for (lambda in lambda_values) {
    #fit glasso and create graph
    graph <- fit_graphical_lasso(S, lambda)
    #plot graph
    plot(
      graph,
      main = bquote(lambda == .(lambda)),
      vertex.label = NA,                  
      vertex.size = 7,                    
      edge.width = 1,                     
      layout = fixed_layout
    )
  }
  title(group_title, outer = TRUE, cex.main = 1.5, font.main = 2)
  par(mfrow = c(1, 1))
}

plot_graphical_lasso(S_normal, lambda_values=lambda_values, group_title = "Graphical Lasso for Normal Patients")

plot_graphical_lasso(S_cancer, lambda_values=lambda_values, group_title = "Graphical Lasso for Cancer Patients")

# By looking at these plots, we select lambda=0.2 and I comment the relative graphs:
lambda = 0.2
# Fit normal graph
fit_normal <- glasso(S_normal, lambda)
glasso_adj_matrix_0 <- as.matrix(fit_normal$wi != 0)
diag(glasso_adj_matrix_0) <- 0
glasso_normal <- graph_from_adjacency_matrix(glasso_adj_matrix_0, mode = "undirected", weighted = TRUE)
# Fit cancer graph
fit_cancer <- glasso(S_cancer, lambda)
glasso_adj_matrix_1 <- as.matrix(fit_cancer$wi != 0)
diag(glasso_adj_matrix_1) <- 0
glasso_cancer <- graph_from_adjacency_matrix(glasso_adj_matrix_1, mode = "undirected", weighted = TRUE)

# Rename
variable_names_0 <- rownames(S_normal)
variable_names_1 <- rownames(S_cancer)
V(glasso_normal)$names <- variable_names_0
V(glasso_cancer)$names <- variable_names_1
# Show the mappings
variable_mapping_0 <- data.frame(Index = 1:length(variable_names_0), Name = variable_names_0)
variable_mapping_1 <- data.frame(Index = 1:length(variable_names_1), Name = variable_names_1)
View(variable_mapping_0)
View(variable_mapping_1)

# Extract precision matrix estimate
precision_matrix_normal <- fit_normal$wi
precision_matrix_cancer <- fit_cancer$wi

# Get graph weights
edge_weights_normal <- precision_matrix_normal[as.matrix(as_edgelist(glasso_normal))]
edge_weights_cancer <- precision_matrix_cancer[as.matrix(as_edgelist(glasso_cancer))]

# Connected components in the normal and the cancer graph
connected_components_normal <- components(glasso_normal)
connected_components_cancer <- components(glasso_cancer)
# Number of connected components
num_components_0 <- connected_components_normal$no
num_components_0
# 52
num_components_1 <- connected_components_cancer$no
num_components_1
# 35
size_components_0 <- connected_components_normal$csize
size_components_0
#A big one (20 nodes)
size_components_1 <- connected_components_cancer$csize
size_components_1
#Again, a big one (37 nodes)

# Generate a color palette with one color per component (if the size is greater than 1)
library(randomcoloR)
set.seed(123)
colors_0 <- distinctColorPalette(num_components_0)
colors_1 <- distinctColorPalette(num_components_1)
connected_colors_0 <- rep("gray", num_components_0)
connected_colors_1 <- rep("gray", num_components_1)
connected_colors_0[size_components_0 > 1] <- colors_0[size_components_0 > 1]
connected_colors_1[size_components_1 > 1] <- colors_1[size_components_1 > 1]

# Plot
# Normal Glasso
plot(
  glasso_normal,
  layout = fixed_layout,            
  vertex.label = NA,
  vertex.size = 7,
  vertex.color = connected_colors_0[connected_components_normal$membership],
  edge.width = abs(edge_weights_normal) * 5,
  main = "Graphical Lasso for Normal Patients"
)


# Cancer Glasso
plot(
  glasso_cancer,
  layout = fixed_layout,            
  vertex.label = NA,
  vertex.size = 7,
  edge.width = abs(edge_weights_cancer) * 5,
  vertex.color = connected_colors_1[connected_components_cancer$membership],
  main = "Graphical Lasso for Cancer Patients"
)



# Differences in the precision matrix estimate

# Print the precision matrix estimates for cancer patients and normal patients
precision_matrix_normal
precision_matrix_cancer
# Exploit the sparsity
sparse_precision_normal <- Matrix(precision_matrix_normal, sparse = TRUE)
sparse_precision_cancer <- Matrix(precision_matrix_cancer, sparse = TRUE)
# Visualize just the non-zero components
summary(sparse_precision_normal)
summary(sparse_precision_cancer)
# 145 non zero elements - normal case, 219 non zero elements - cancer case
non_zero_normal <- as.data.frame(summary(sparse_precision_normal))
colnames(non_zero_normal) <- c("Row", "Column", "Value")
non_zero_cancer <- as.data.frame(summary(sparse_precision_cancer))
colnames(non_zero_cancer) <- c("Row", "Column", "Value")
View(non_zero_normal)
View(non_zero_cancer)

# Create a list of the non zero components that are in non_zero_cancer 
# but not in non_zero_normal and viceversa
non_zero_normal$Index <- paste(non_zero_normal$Row, non_zero_normal$Column, sep = ",")
non_zero_cancer$Index <- paste(non_zero_cancer$Row, non_zero_cancer$Column, sep = ",")
# Find indices that are in cancer but not in normal
diff_indices_0 <- setdiff(non_zero_cancer$Index, non_zero_normal$Index)
diff_indices_0
# 122 couples
# Viceversa, the indexes that are in normal but not in cancer
diff_indices_1 <- setdiff(non_zero_normal$Index, non_zero_cancer$Index)
diff_indices_1
# 48 couples

##########

