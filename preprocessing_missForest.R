library(tidyverse)
library(caret)
library(superheat)
library(missForest)
library(ggridges)

data_train <- read.table("../data/ovary_cancer_update171205_train2.txt",sep="\t",header=T)
data_test <- read.table("../data/ovary_cancer_update171205_test2.txt",sep="\t",header=T)

exp_vals <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")
# data_train %>% str()

use_exp_vals <- c()

for (exp_val in exp_vals){
  print(exp_val)
  if (length(na.omit(data_train[,exp_val])) >= 170){
    use_exp_vals <- c(use_exp_vals,exp_val)
  }
}

data_train <- data_train %>% 
  mutate(flag = rep("train", nrow(.)))

data_test <- data_test %>% 
  mutate(flag = rep("test", nrow(.)))

data_all <-
  data_train %>% rbind(., data_test) %>% 
  filter(FIGO!="I") %>% 
  select(ID, use_exp_vals, residual_cancer, hist_type, flag) %>% 
  # na.omit() %>% # To omit rows with NA value(s), use this line.
  filter(!is.na(residual_cancer))

subtype <- "all"
png(paste("./data_NA_pattern_", subtype, ".png", sep=""), height = 600, width = 600)
superheat::superheat(data_all %>% select(-c(ID, residual_cancer, hist_type, flag)) %>% as.matrix(), 
                     membership.rows = data_all$flag, 
                     # scale the matrix
                     scale = T,
                     bottom.label.text.angle = 90, 
                     # change color of missing values
                     heat.na.col = "red")
dev.off()

png(paste("./data_NA_pattern_hist_type_", subtype, ".png", sep=""), height = 600, width = 600)
superheat::superheat(data_all %>% select(-c(ID, residual_cancer, hist_type, flag)) %>% as.matrix(), 
                     membership.rows = data_all$hist_type, 
                     # scale the matrix
                     scale = T,
                     bottom.label.text.angle = 90, 
                     # change color of missing values
                     heat.na.col = "red")
dev.off()

# Check variables with a few unique values: 
data_all %>% select(-c(ID, residual_cancer, hist_type, flag)) %>% 
  summarise_all(funs(n_distinct(., na.rm = TRUE))) %>% 
  gather(key = var, value = n_unique) %>% 
  arrange(n_unique) %>% head()
# Convert variable(s) with only a few unique values
# into categorical variable(s) (e.g., factor).
data_all$Baso <- as.factor(data_all$Baso)

set.seed(123456)
res_mF <- missForest(data_all %>% select(-c(ID, residual_cancer, hist_type, flag)))
res_mF$ximp$Baso <- res_mF$ximp$Baso %>% as.character() %>% as.numeric()
data_all_mF <- res_mF$ximp

png(paste("./data_NA_imputed_", subtype, ".png", sep=""), height = 600, width = 600)
superheat::superheat(data_all_mF %>% as.matrix(), 
                     membership.rows = data_all$flag, 
                     # scale the matrix
                     scale = T,
                     bottom.label.text.angle = 90, 
                     # change color of missing values
                     heat.na.col = "red")
dev.off()

png(paste("./data_NA_imputed_hist_type_", subtype, ".png", sep=""), height = 600, width = 600)
superheat::superheat(data_all_mF %>% as.matrix(), 
                     membership.rows = data_all$hist_type, 
                     # scale the matrix
                     scale = T,
                     bottom.label.text.angle = 90, 
                     # change color of missing values
                     heat.na.col = "red")
dev.off()

data_all_mF$flag <- data_all$flag

tmp <- data_all %>% select(-c(ID, residual_cancer, hist_type, flag))
tmp$Baso <- as.numeric(tmp$Baso)
tmp <- tmp %>% scale(center = TRUE, scale = TRUE)
tmp2 <- data_all_mF %>% select(-flag) %>%
  scale(center = TRUE, scale = TRUE)
tmp2[!is.na(tmp)] <- NA
dat_concat <- rbind(tmp %>% as.data.frame(), tmp2 %>% as.data.frame())
dat_concat$flag <- rep(c("obs", "pred"), each = nrow(tmp))

p <- dat_concat %>% gather(key = var, value = value, -flag) %>% 
  na.omit() %>% 
  ggplot(aes(x=value, y=var, fill = flag)) +
  geom_density_ridges(alpha = .2) + 
  coord_cartesian(xlim = c(-3, 3)) + 
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape") + 
  labs(x = "Standardized value", y = "Features", title = "Distributions of missForest prediction and observed data") + 
  theme_bw()
ggsave(file = paste("./ggridegs_pred_obs_", subtype, ".pdf", sep = ""), plot = p, useDingbats = FALSE)

data_all <- data_all %>% mutate(imp = rowSums(is.na(.)))
data_all$imp <- ifelse(data_all$imp > 0, "imp", "obs")

data_all_mF <- data_all %>% 
  select(c(ID, residual_cancer, hist_type, imp)) %>% 
  cbind(data_all_mF, .)
data_all_mF$residual_cancer <- 
  factor(data_all_mF$residual_cancer, levels = c("0", "1", "2"))

# colnames(data_all_mF)

data_all_mF %>%
  write.table(file="./data_all_mF_180604.txt",sep="\t",quote=FALSE, row.names = FALSE)

data_train_merged_mF_all <- 
  data_all_mF %>% 
  filter(flag == "train") %>% 
  select(-c(flag, imp, hist_type, ID))
data_test_merged_mF_all <- 
  data_all_mF %>% 
  filter(flag == "test") %>% 
  select(-c(flag, imp, hist_type, ID))

data_train_merged_mF_serous <- 
  data_all_mF %>% 
  filter(flag == "train", hist_type == "serous") %>% 
  select(-c(flag, imp, hist_type, ID))
data_test_merged_mF_serous <- 
  data_all_mF %>% 
  filter(flag == "test", hist_type == "serous") %>% 
  select(-c(flag, imp, hist_type, ID))

# data_train_merged_mF_serous <- data_train_merged_mF[data_train_merged$hist_type == "serous", ]
# data_test_merged_mF_serous <- data_test_merged_mF[data_test_merged$hist_type == "serous", ]
table_imp <- 
  rbind(
    data_all_mF %>% 
      filter(flag == "train", imp == "obs") %>% 
      .$residual_cancer %>% 
      table() %>% 
      as.data.frame() %>% 
      mutate(flag = rep("Train_omitNA", 3)), 
    data_all_mF %>% 
      filter(flag == "train") %>% 
      .$residual_cancer %>% 
      table() %>% 
      as.data.frame() %>% 
      mutate(flag = rep("Train_missForest", 3)), 
    data_all_mF %>% 
      filter(flag == "test", imp == "obs") %>% 
      .$residual_cancer %>% 
      table() %>% 
      as.data.frame() %>% 
      mutate(flag = rep("Test_omitNA", 3)), 
    data_all_mF %>% 
      filter(flag == "test") %>% 
      .$residual_cancer %>% 
      table() %>% 
      as.data.frame() %>% 
      mutate(flag = rep("Test_missForest", 3))
  )

colnames(table_imp)[1] <- "class"
table_imp$class <- paste("R", table_imp$class, sep = "")
table_imp$flag <- factor(table_imp$flag, levels = c("Train_omitNA", "Train_imputed", "Test_omitNA", "Test_imputed"))

p <- table_imp %>% ggplot(aes(x = flag, y = Freq, fill = class)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() +
  labs(x = "Preprocessing type")
ggsave(file = paste("./barplot_class_freq_vs_preprocessing_", subtype, ".pdf", sep = ""), plot = p, useDingbats = FALSE)

# To save subsets of the data by spesifying subtype and train/test: 
# subtype = "all"
# data_train_merged_mF %>% 
#   write.table(file=paste("./data_train_merged_mF_", subtype,".txt",sep=""),sep="\t",quote=FALSE, row.names = FALSE)
# 
# data_test_merged_mF %>% 
#   write.table(file=paste("./data_test_merged_mF_", subtype,".txt",sep=""),sep="\t",quote=FALSE, row.names = FALSE)
# 
# subtype = "serous"
# data_train_merged_mF_serous %>% 
#   write.table(file=paste("./data_train_merged_mF_", subtype,".txt",sep=""),sep="\t",quote=FALSE, row.names = FALSE)
# 
# data_test_merged_mF_serous %>% 
#   write.table(file=paste("./data_test_merged_mF_", subtype,".txt",sep=""),sep="\t",quote=FALSE, row.names = FALSE)
