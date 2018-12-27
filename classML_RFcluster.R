library(caret)
#library(caretEnsemble)
library(pROC)
library(gplots)
#library(doMC)

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
                           repeats = 10,
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
class <- "MDScluster_early"
margin <- 1
#order <- c("benign","early","late")
#label <- c("benign","early","late")
#order <- c("benign","I","II","III","IV")
#label <- c("benign","I","II","III","IV")
#order <- c("clearcell","serous","endometrioid","mucinous")
#label <- c("clear\ncell","serous","endo\nmetrioid","mucinous")
order <- c("C2","C1")
label <- c("C2","C1")

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
		if (val == 1){
			tmp_vec <- c(tmp_vec,"C1")
		} else {
			tmp_vec <- c(tmp_vec,"C2")
		}
	}

	return(tmp_vec)
}

data_train <- read.table("../data/ovary_cancer_update170830_dummyvar_train.txt",sep="\t",header=T,row.names=1)
data_test <- read.table("../data/ovary_cancer_update170830_dummyvar_test.txt",sep="\t",header=T,row.names=1)

exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")
#exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")
#exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP")

#data_all <- data_all[data_all$hist_type!="others",]
#data_all <- data_all[data_all$hist_type!="cyst",]
#data_all <- data_all[data_all$hist_type!="teratoma",]

# for (i in 1:length(val_type)){
# 	if (val_type[i] == "log"){
# 		data_all[,i] <- log(data_all[,i])
# 	}
# }

use_exp_vars <- c()

for (exp_var in exp_vars){
	print(exp_var)
	if (length(na.omit(data_train[,exp_var])) >= 170){
		use_exp_vars <- c(use_exp_vars,exp_var)
	}
}

data_train_exp <- data_train[data_train$stage=="early",use_exp_vars]
data_test_exp <- data_test[data_train$stage=="early",use_exp_vars]

#dir.create(paste('./classMLv6_',class,sep=''))
#setwd(paste('./classMLv6_',class,sep=''))

MDSdata <- read.table('MDSdata_cluster.txt')

data_train_merged <- na.omit(data.frame(data_train_exp,class=MDSdata[rownames(data_train_exp),"clust"]))
data_train_merged$class <- factor(convert_class(data_train_merged$class))
data_train_merged <- na.omit(data_train_merged)

data_test_merged <- na.omit(data.frame(data_test_exp,class=MDSdata[rownames(data_test_exp),"clust"]))
data_test_merged$class <- factor(convert_class(data_test_merged$class))
data_test_merged <- na.omit(data_test_merged)

accuracy_list <- c()
auc_list <- c()

for (method in method_list){

	print(method)

	if (method != "glm"){

		set.seed(123456)

		e <- try(model <- train(class ~ ., data = data_train_merged, 
	         method = method, 
	         trControl = fitControl,
	         maximize = FALSE,
	         preProcess = c('center', 'scale'),
	         metric = "Dist"),silent=FALSE)
	} else {

		set.seed(123456)
		
		e <- try(model <- train(class ~ ., data = data_train_merged, 
	         method = "glm",
	         family = "binomial", 
	         trControl = fitControl,
	         maximize = FALSE,
	         preProcess = c('center', 'scale'),
	         metric = "Dist"),silent=FALSE)
	}
	
	if (class(e) != "try-error") {

		pred<-predict(model, data_test_merged, type="prob")
		roc.res <- roc(data_test_merged$class, pred$C1)
		auc_list <- c(auc_list,pROC::auc(roc.res))
		rocplot <- paste("ovary_ROCplot_",method,"_",k,"-fold_",class,".pdf",sep='')
		pdf(rocplot,useDingbats=FALSE)
		plot.roc(roc.res,col=ROC_colors[1])
		#sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=ROC_colors[i]))
		dev.off()
		# pred<-predict(model, data_test_merged)
		# res <- table(data_test_merged[,"class"],pred)
		# #print(res)
		# res2 <- filltable(res)
		# print(calc_accuracy(res2))
		res.coords <- coords(roc.res,"best", best.method="closest.topleft",ret=c("threshold","accuracy"))
		accuracy_list <- c(accuracy_list,res.coords["accuracy"])

		#print(res2)
		colorvec <- setcolor(res2)
		blnplot <- paste("ovary_balloonplot_",method,"_",k,"-fold_",class,".pdf",sep='')
		pdf(blnplot,useDingbats=FALSE)
		balloonplot(t(res2),show.zeros=TRUE,main=NULL,rowmar=margin,text.size=1,label.size=1+t(res2)/max(t(res2)),dotcolor=colorvec,show.margins=TRUE, xlab="prediction",ylab="category",cum.margins=FALSE)
		dev.off()

		#print(res)
		e <- try(imp <- varImp(model,scale = TRUE))
		#rownames(imp) <- gsub(paste("imp",i,".",sep=""),"",rownames(imp))
		if (class(e) != "try-error") {
			imp_file <- paste(project,class,"_",method,"importance.txt",sep='')
			write.table(imp$importance,file=imp_file,sep="\t",quote=FALSE)
		}
		e <- try(imp <- varImp(model,scale = FALSE))
		#rownames(imp) <- gsub(paste("imp",i,".",sep=""),"",rownames(imp))
		if (class(e) != "try-error") {
			imp_file <- paste(project,class,"_",method,"importance_raw.txt",sep='')
			write.table(imp$importance,file=imp_file,sep="\t",quote=FALSE)
		}
	} else {
		auc_list <- c(auc_list,NA)
		accuracy_list <- c(accuracy_list,NA)
	}
}
acc_result <- data.frame(method=method_list,accuracy=accuracy_list,auc=auc_list)
acc_result_file <- paste(class,"accuracy_",k,"-fold_validation.txt",sep="")
write.table(acc_result,file=acc_result_file,sep="\t",quote=FALSE)
