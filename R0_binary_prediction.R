library(caret)
library(parallel)
library(doParallel)

train <- 
  data_all_mF %>% 
  filter(flag == "train") %>% 
  select(-c(flag, imp, hist_type, ID))
test <- 
  data_all_mF %>% 
  filter(flag == "test") %>% 
  select(-c(flag, imp, hist_type, ID))

train$Class <- ifelse(train$residual_cancer == 0, "R0", "Rest") %>% as.factor()
test$Class <- ifelse(test$residual_cancer == 0, "R0", "Rest") %>% as.factor()

train <- train %>% select(-residual_cancer)
test <- test %>% select(-residual_cancer)

train$Class %>% table()
test$Class %>% table()
train$Class %>% class()

fourStats <- function (data, lev = levels(data$obs), model = NULL) {
  ## This code will get use the area under the ROC curve and the
  ## sensitivity and specificity values using the current candidate
  ## value of the probability threshold.
  out <- c(twoClassSummary(data, lev = levels(data$obs), model = NULL))
  
  ## The best possible model has sensitivity of 1 and specificity of 1. 
  ## How far are we from that value?
  coords <- matrix(c(1, 1, out["Spec"], out["Sens"]),
                   ncol = 2,
                   byrow = TRUE)
  colnames(coords) <- c("Spec", "Sens")
  rownames(coords) <- c("Best", "Current")
  c(out, Dist = dist(coords)[1])
}

fitControl <- trainControl(method = "repeatedcv", 
                           number = 10, repeats = 100, 
                           classProbs = TRUE, 
                           summaryFunction = fourStats, 
                           allowParallel = TRUE)

tGrid  <-  expand.grid(mtry= (1:6)*1)

cluster <- makeCluster(3)
registerDoParallel(cluster)

t <- proc.time()
set.seed(123456)
model.tune <- caret::train(Class ~ ., data = train, 
                           method = "rf", 
                           trControl = fitControl, 
                           tuneGrid = tGrid,
                           ntree = 1000, 
                           metric = "ROC")
                           # metric = "Dist", 
                           # maximize = FALSE)
proc.time() - t
# user  system elapsed 
# 6.109   0.913 392.058
stopCluster(cluster)
registerDoSEQ()

set.seed(123456)
rf_best <- caret::train(Class ~ ., data = train, 
                           method = "rf", 
                           trControl = trainControl(method="none"), 
                           tuneGrid = data.frame(mtry=1),
                           ntree = 5000)

# caret::confusionMatrix(data = test$Class, predict(rf_best, newdata = test))
table(test$Class, predict(rf_best, newdata = test))

pred <- predict(rf_best, newdata = test, type = "prob")
roc_0vsrest <- pROC::roc(test$Class, 1 - pred$R0)
auc(roc_0vsrest)
plot.roc(roc_0vsrest, col="red", add = TRUE)
ci(auc(roc_0vsrest))

plot(varImp(rf_best))

seed_list <- gtools::permutations(6, 6) %>% 
  apply(1, paste, collapse="") %>% 
  as.numeric()

df_varimp <- matrix(rep(NA, 33*100), ncol = 33) %>% as.data.frame
colnames(df_varimp) <- colnames(train %>% select(-Class))
df_varimp %>% head()
t <- proc.time()
for(k in 1:100){
  set.seed(seed_list[k])
  rf_best <- caret::train(Class ~ ., data = train, 
                          method = "rf", 
                          trControl = trainControl(method="none"), 
                          tuneGrid = data.frame(mtry=1),
                          ntree = 5000)
  if(k == 1){
    pred <- predict(rf_best, newdata = test, type = "prob")
    
    pred %>% mutate(ID = data_all_mF %>% 
                      filter(flag == "test") %>% 
                      .$ID) %>% 
      write.table("./data/rf_best_pred_prob_R0vsRest.txt", row.names = FALSE, quote = FALSE)
    roc_0vsrest <- pROC::roc(test$Class, 1 - pred$R0)
    auc(roc_0vsrest)
    pdf(file = paste("./plot/rf_best_ROC_R0vsRest_AUC_", auc(roc_0vsrest), ".pdf", sep = ""), useDingbats = FALSE)
    plot.roc(roc_0vsrest, col="red")
    dev.off()
  }
  df_varimp[k, ] <- varImp(rf_best, scale = FALSE)$importance %>% 
    unlist() %>% 
    as.numeric()
}
proc.time() - t

df_varimp %>% 
  write.table("./data/rf_best_varimp_R0vsRest.txt", row.names = FALSE, quote = FALSE)

df_varimp_long <- df_varimp %>% mutate(iter = 1:100) %>% 
  gather(key = var, val = varimp, -iter)

median_var_order <- df_varimp_long %>% 
  group_by(var) %>% 
  summarise(med = median(varimp)) %>% 
  arrange(-med) %>% 
  .$var

df_varimp_long$var <- factor(df_varimp_long$var, levels = median_var_order)

p <- df_varimp_long %>% 
  ggplot(aes(x = var, y = varimp, fill = var)) + 
  geom_boxplot(outlier.size=0.5) +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=10, angle = 90, hjust = 1)
  )
ggsave(file = "./plot/varimp_boxplot_R0vsRest_allvar.pdf", plot = p, useDingbats = FALSE)
  
p <- df_varimp_long %>% 
  filter(var %in% median_var_order[1:8]) %>% 
  ggplot(aes(x = var, y = varimp, fill = var)) + 
  geom_boxplot(outlier.size=0.5) +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=10, angle = 90, hjust = 1)
  )
ggsave(file = "./plot/varimp_boxplot_R0vsRest_top8var.pdf", plot = p, useDingbats = FALSE)
