# Predict Ovary residual cancer using Ordinal Forest
library(ordinalForest)
library(dplyr)
data_train <- read.table("../data/ovary_cancer_update171205_train.txt",sep="\t",header=T)
data_test <- read.table("../data/ovary_cancer_update171205_test.txt",sep="\t",header=T)

exp_vals <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")

use_exp_vals <- c()

for (exp_val in exp_vals){
  print(exp_val)
  if (length(na.omit(data_train[,exp_val])) >= 170){
    use_exp_vals <- c(use_exp_vals,exp_val)
  }
}

data_train_exp <- data_train[data_train$FIGO!="I",use_exp_vals]
data_test_exp <- data_test[data_test$FIGO!="I",use_exp_vals]

data_train_merged <- na.omit(data.frame(data_train_exp,residual_cancer=data_train[data_train$FIGO!="I","residual_cancer"]))
data_test_merged <- na.omit(data.frame(data_test_exp,residual_cancer=data_test[data_test$FIGO!="I","residual_cancer"]))

data_train_merged$residual_cancer <- factor(data_train_merged$residual_cancer, levels = c("0", "1", "2"))

ordforres <- ordfor(depvar="residual_cancer", data=data_train_merged, nsets=1000, ntreeperdiv=100,
                    ntreefinal=5000, perffunction="oneclass", classimp = "1")

ordforres
# Study variable importance values:
sort(ordforres$varimp, decreasing=TRUE)
# Take a closer look at the top variables:
boxplot(data_train_merged$Neu ~ data_train_merged$residual_cancer, horizontal=TRUE)

# Predict values of the ordinal target variable in the test dataset:
preds <- predict(ordforres, newdata=data_test_merged)
preds
# Compare predicted values with true values:
table(data.frame(true_values=data_test_merged$residual_cancer, predictions=preds$ypred))

# Predicted ordered categories
head(preds$ypred)

# Prediction probabilities
head(preds$classfreqtree)
# preds$classfreqtreeg

df <- as.data.frame(preds$classfreqtree)
colnames(df) <- c("R0", "R1", "R2")
df$id <- 1:dim(df)[1]
df <- gather(df, key = pred, value = prob, -id) %>% 
  arrange(id) # %>% head()

p <- ggplot(df,aes(as.factor(id),as.factor(pred)))+
  geom_tile(aes(fill = prob),colour="white")+
  scale_fill_gradient(low="white",high="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "ID", y = "Probability")
ggsave(file = pdf,plot=p,width=10,height=7)
