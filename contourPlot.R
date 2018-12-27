library(caret)
library(scales)

method_list = c("gbm","rf","svmRadial","nnet","nb","cforest","glmnet","glm")

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

fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10,
                           savePredictions = TRUE,
                           classProbs = TRUE,
                           summaryFunction = fourStats)

convert_class <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == "benign"){
			tmp_vec <- c(tmp_vec,"benign")
		} else {
			tmp_vec <- c(tmp_vec,"malignant")
		}
	}

	return(tmp_vec)
}

convert_class2 <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == "benign"){
			tmp_vec <- c(tmp_vec,0)
		} else {
			tmp_vec <- c(tmp_vec,1)
		}
	}

	return(tmp_vec)
}

data_train <- read.table("../data/ovary_cancer_update170830_dummyvar_train.txt",sep="\t",header=T)
data_test <- read.table("../data/ovary_cancer_update170830_dummyvar_test.txt",sep="\t",header=T)

#exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso","Blood_A","Blood_B","Blood_AB","Blood_O")
exp_vars <- c("age","CA125")
#exp_vars <- c("age","CA125","Fbg","LDH","Alb","CRP","APTT")

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

data_train_exp <- data_train[,use_exp_vars]
data_test_exp <- data_test[,use_exp_vars]

data_train_merged <- na.omit(data.frame(data_train_exp,class=data_train[,"stage"]))
data_train_merged$class <- factor(convert_class(data_train_merged$class))
data_train_merged <- na.omit(data_train_merged)

data_test_merged <- na.omit(data.frame(data_test_exp,class=data_test[,"stage"]))
data_test_merged$class <- factor(convert_class(data_test_merged$class))
data_test_merged <- na.omit(data_test_merged)

px<-seq(from=21,to=90,length.out=300)
py<-seq(from=0.5,to=5,length.out=300)
pgrid<-expand.grid(px,py)
names(pgrid)<-c("age","CA125")

for (method in method_list){

	print(method)

	if (method != "glm"){

		e <- try(model <- train(class ~ ., data = data_train_merged, 
	         method = method, 
	         trControl = fitControl,
	         maximize = FALSE,
	         preProcess = c('center', 'scale'),
	         metric = "Dist"),silent=FALSE)
	} else {
		e <- try(model <- train(class ~ ., data = data_train_merged, 
	         method = "glm",
	         family = "binomial", 
	         trControl = fitControl,
	         maximize = FALSE,
	         preProcess = c('center', 'scale'),
	         metric = "Dist"),silent=FALSE)
	}

	if (class(e) != "try-error") {

		pred.border <- predict(model,newdata=pgrid)
		pred.border <-convert_class2(pred.border)

		pdf(paste("ovary_cancer",method,"pred_border_train.pdf",sep="_"),height=5,width=5,useDingbats=FALSE)

		plot(data_train_merged[data_train_merged$class=="benign",c("age","CA125")],col=alpha("#3D3939",0.7),pch=16,cex=1,xlim=c(21,90),ylim=c(0,5),alpha=0.7)
		points(data_train_merged[data_train_merged$class=="malignant",c("age","CA125")],col=alpha("#D5221E",0.7),pch=16,cex=1)

		par(new=T)
		contour(px,py,array(pred.border,dim=c(length(px),length(py))),xlim=c(21,90),ylim=c(0,5),col=alpha("purple",0.8),lwd=1.5,drawlabels=F,levels=0.5)

		dev.off()

		pdf(paste("ovary_cancer",method,"pred_border_test.pdf",sep="_"),height=5,width=5,useDingbats=FALSE)

		plot(data_test_merged[data_test_merged$class=="benign",c("age","CA125")],col=alpha("#3D3939",0.7),pch=16,cex=1,xlim=c(21,90),ylim=c(0,5),alpha=0.7)
		points(data_test_merged[data_test_merged$class=="malignant",c("age","CA125")],col=alpha("#D5221E",0.7),pch=16,cex=1)

		par(new=T)
		contour(px,py,array(pred.border,dim=c(length(px),length(py))),xlim=c(21,90),ylim=c(0,5),col=alpha("purple",0.8),lwd=1.5,drawlabels=F,levels=0.5)

		dev.off()
	}
}


