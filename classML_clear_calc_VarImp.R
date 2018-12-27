library(caret)
#library(caretEnsemble)
library(pROC)
library(gplots)
library(doMC)
library(party)
library(randomForest)
library(gbm)
library(dplyr)

registerDoMC(cores=18)

#nthreads = as.numeric(commandArgs(trailingOnly=TRUE)[1])
#method_list = c("gbm","rf","svmRadial","nnet","nb","cforest","glmnet","glm")
method_list = c("rf")

ROC_colors <- c("#C7243A","#3261AB","#007FB1","#009250","#EDAD0B")

fourStats <- function (data, lev = levels(data$obs), model = NULL) {
  ## This code will get use the area under the ROC curve and the
  ## sensitivity and specificity values using the current candidate
  ## value of the probability threshold.
  out <- c(twoClassSummary(data, lev = levels(data$obs), model = NULL))
 
  ## The best possible model has sensitivity of 1 and specifity of 1. 
  ## How far are we from that value?
  coords <- matrix(c(1, 1, out["Spec"], out["Sens"]), 
                   ncol = 2, 
                   byrow = TRUE)
  colnames(coords) <- c("Spec", "Sens")
  rownames(coords) <- c("Best", "Current")
  c(out, Dist = dist(coords)[1])
}


#fitControl <- trainControl(## 
#                           method = "boot",
#                           number = 100)

fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 100,
                           savePredictions = TRUE,
                           classProbs = TRUE,
                           summaryFunction = fourStats)

# caret_package_list <- "caret_package_list.txt"

# caret_package_info <- read.table(caret_package_list,sep='\t',header=T)


#method_list <- caret_package_info[caret_package_info$type == "Classification" | caret_package_info$type == "Classification, Regression",]$method

#method_list <- method_list[-which(method_list %in% "xgbLinear")]
#method_list <- method_list[-which(method_list %in% "xgbTree")]
#method_list <- method_list[-which(method_list %in% "gaussprLinear")]
#method_list <- method_list[-which(method_list %in% "smda")]
#method_list <- method_list[-which(method_list %in% "randomGLM")]
#method_list <- method_list[-which(method_list %in% "gaussprPoly")]
#method_list <- c("xgbLinear","xgbTree",as.vector(method_list))

project <- "ovary"
class <- "clear"
margin <- 1
#order <- c("benign","early","late")
#label <- c("benign","early","late")
#order <- c("benign","I","II","III","IV")
#label <- c("benign","I","II","III","IV")
#order <- c("clearcell","serous","endometrioid","mucinous")
#label <- c("clear\ncell","serous","endo\nmetrioid","mucinous")
order <- c("clear","non_clear")
label <- c("clear","non_clear")

k <- 2  #k-fold cross validation

filltable <- function(res){
	res2 <- res
	cnames <- colnames(res)
	l=1
	for (j in 1:nrow(res)){
		if (l <= ncol(res)){
			if (rownames(res)[j] != colnames(res)[l]){
				res2 <- cbind(res2,rep(0,nrow(res)))
				cnames <- c(cnames,rownames(res)[j])
			}else{
				l = l+1
			}
		}else{
			res2 <- cbind(res2,rep(0,nrow(res)))
			cnames <- c(cnames,rownames(res)[j])
		}
	}
	colnames(res2) <- cnames
	#print(colnames(res2))
	res2 <- res2[order,order]
	rownames(res2) <- label
	colnames(res2) <- label

	res2 <- as.table(res2)

	return(res2)
}

setcolor <- function(res){
	colorvec <- c("palegreen3")
	for (l in 1:nrow(res)-1){
		colorvec <- c(colorvec,rep("gray",nrow(res)),"palegreen3")
	}
	return(colorvec)
}

calc_accuracy <- function(res){
	tmp_acc = 0
	for (l in 1:nrow(res)){
		tmp_acc = tmp_acc + res[l,l]
	}
	tmp_acc = tmp_acc / sum(res)

	return(tmp_acc)
}

convert_class <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == "clearcell"){
			tmp_vec <- c(tmp_vec,"clear")
		} else if (val == "others" | val == "cyst" | val == "teratoma") {
			tmp_vec <- c(tmp_vec,NA)
		} else {
			tmp_vec <- c(tmp_vec,"non_clear")
		}
	}

	return(tmp_vec)
}


data_train <- read.table("../data/ovary_cancer_update170830_dummyvar_train.txt",sep="\t",header=T)
data_test <- read.table("../data/ovary_cancer_update170830_dummyvar_test.txt",sep="\t",header=T)

exp_vals <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")

#data_all <- data_all[data_all$hist_type!="others",]
#data_all <- data_all[data_all$hist_type!="cyst",]
#data_all <- data_all[data_all$hist_type!="teratoma",]

# for (i in 1:length(val_type)){
# 	if (val_type[i] == "log"){
# 		data_all[,i] <- log(data_all[,i])
# 	}
# }

use_exp_vals <- c()

for (exp_val in exp_vals){
	print(exp_val)
	if (length(na.omit(data_train[,exp_val])) >= 170){
		use_exp_vals <- c(use_exp_vals,exp_val)
	}
}

data_train_exp <- data_train[,use_exp_vals]
data_test_exp <- data_test[,use_exp_vals]

#dir.create(paste('./classMLv6_',class,sep=''))
#setwd(paste('./classMLv6_',class,sep=''))

files <- list.files()

data_train_merged <- na.omit(data.frame(data_train_exp,class=data_train[,"hist_type"]))
data_train_merged$class <- factor(convert_class(data_train_merged$class))
data_train_merged <- na.omit(data_train_merged)

data_test_merged <- na.omit(data.frame(data_test_exp,class=data_test[,"hist_type"]))
data_test_merged$class <- factor(convert_class(data_test_merged$class))
data_test_merged <- na.omit(data_test_merged)

accuracy_list <- c()
auc_list <- c()

list_seeds <- gtools::permutations(6, 6) %>% 
  apply(1, paste, collapse="") %>% 
  as.numeric()

for (method in method_list){

	var_imp <- use_exp_vars

	for (seed in list_seeds[1:100]){

		if (method != "glm"){

			set.seed(seed)

			e <- try(model <- train(class ~ ., data = data_train_merged, 
		         method = method, 
		         trControl = fitControl,
		         maximize = FALSE,
		         preProcess = c('center', 'scale'),
		         metric = "Dist"),silent=FALSE)
		} else {

			set.seed(seed)

			e <- try(model <- train(class ~ ., data = data_train_merged, 
		         method = "glm",
		         family = "binomial", 
		         trControl = fitControl,
		         maximize = FALSE,
		         preProcess = c('center', 'scale'),
		         metric = "Dist"),silent=FALSE)
		}

	if (class(e) != "try-error") {
			#print(res)
			e <- try(imp <- varImp(model,scale = FALSE))
			#rownames(imp) <- gsub(paste("imp",i,".",sep=""),"",rownames(imp))
			if (class(e) != "try-error") {
				var_imp <- cbind(var_imp,imp$importance)
			}
			# e <- try(imp <- varImp(model,scale = FALSE))
			# #rownames(imp) <- gsub(paste("imp",i,".",sep=""),"",rownames(imp))
			# if (class(e) != "try-error") {
			# 	imp_file <- paste(project,class,"_",method,"importance_raw.txt",sep='')
			# 	write.table(imp$importance,file=imp_file,sep="\t",quote=FALSE)
			# }
		}
	}

	imp_file <- paste(project,class,"_",method,"importance.txt",sep='')
	write.table(var_imp,file=imp_file,sep="\t",quote=FALSE,row.names=F)
}
