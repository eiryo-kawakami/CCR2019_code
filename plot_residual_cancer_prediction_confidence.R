library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pROC)
library(gplots)

margin <- 1
order <- c("benign","I","II","III","IV")

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
	#res2 <- res2[order,order]
	# rownames(res2) <- label
	# colnames(res2) <- label

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



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

is_max <- function(df){
	res <- c()
	for (i in 1:nrow(df)){
		tmp_res <- c()
		for (j in 1:ncol(df)){
			if (df[i,j] == max(df[i,])){
				tmp_res <- c(tmp_res,2)
			} else {
				tmp_res <- c(tmp_res,1)
			}
		}
		res <- rbind(res,tmp_res)
	}
	return(res)
}

as_max <- function(df){
	res <- c()
	for (i in 1:nrow(df)){
		for (j in 1:ncol(df)){
			if (df[i,j] == max(df[i,])){
				res <- c(res,colnames(df)[j])
			} 
		}
	}
	return(res)
}


convert_residual_class1 <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == 0){
			tmp_vec <- c(tmp_vec,"R0")
		} else {
			tmp_vec <- c(tmp_vec,"R1")
		}
	}
	return(tmp_vec)
}

convert_residual_class2 <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == 2){
			tmp_vec <- c(tmp_vec,"R2")
		} else {
			tmp_vec <- c(tmp_vec,"R0")
		}
	}
	return(tmp_vec)
}

rowMax <- function(x, na.rm=TRUE) {
    apply(x, 1, max, na.rm=na.rm)
}

convert_class <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == 0){
			tmp_vec <- c(tmp_vec,"R0")
		} else if (val == 1) {
			tmp_vec <- c(tmp_vec,"R1")
		} else {
			tmp_vec <- c(tmp_vec,"R2")
		}
	}

	return(tmp_vec)
}



cols1 <- c(brewer.pal(5, "Set1")[c(1,2,3,5)],"#555555")
cols2 <- gg_color_hue(3)[c(3,2,1)]
cols3 <- c("#999999","#E69F00","#56B4E9")

df1 <- read.table('SMOTE_ordinalForest_prediction_probability.txt',sep=" ",header=T,row.names=4)
data_train <- read.table("../data/ovary_cancer_update171205_train.txt",sep="\t",header=T,row.names=1)
data_test <- read.table("../data/ovary_cancer_update171205_test.txt",sep="\t",header=T,row.names=1)

exp_vals <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")

use_exp_vals <- c()

for (exp_val in exp_vals){
  print(exp_val)
  if (length(na.omit(data_train[,exp_val])) >= 170){
    use_exp_vals <- c(use_exp_vals,exp_val)
  }
}

data_test_exp <- data_test[data_test$FIGO!="I",use_exp_vals]
data_test_merged <- data.frame(data_test_exp,data_test[data_test$FIGO!="I",c("residual_cancer","FIGO","hist_type")])

df_merge <- data.frame(df1,data_test_merged[as.vector(rownames(df1)),c("FIGO","hist_type","residual_cancer")])
df_merge$hist_type <- factor(df_merge$hist_type,levels=c("serous","clearcell","endometrioid","mucinous","others"))

df_merge$resid <- convert_class(df_merge$residual_cancer)

df_merge <- df_merge[order(df_merge$hist_type),]
df_merge <- df_merge[order(df_merge$FIGO),]
df_max <- as.data.frame(is_max(df_merge[,c(1:3)]))
df_max$pred <- as_max(df_merge[,c(1:3)])
colnames(df_max) <- c("R0_is_max","R1_is_max","R2_is_max","pred")
df_merge <- data.frame(df_merge,df_max)

res <- table(df_merge$resid,df_merge$pred)
res2 <- filltable(res)

colorvec <- setcolor(res2)
blnplot <- paste("ovary_balloonplot_ordinalForest_2-fold_residual_cancer_all.pdf",sep='')
pdf(blnplot,useDingbats=F)
balloonplot(t(res2),show.zeros=TRUE,main=NULL,rowmar=margin,text.size=1,label.size=1+t(res2)/max(t(res2)),dotcolor=colorvec,show.margins=TRUE, xlab="prediction",ylab="category",cum.margins=FALSE)
dev.off()


df_merge$resid1 <- factor(convert_residual_class1(df_merge$residual_cancer),levels=c("R0","R1"))
df_merge$resid2 <- factor(convert_residual_class2(df_merge$residual_cancer),levels=c("R0","R2"))

roc.res <- roc(df_merge$resid1, df_merge$R1+df_merge$R2)
#auc_list <- c(auc_list,pROC::auc(roc.res))
rocplot <- paste("ovary_ROCplot_ordinalForest_2-fold_residual_cancer_prediction1.pdf",sep='')
pdf(rocplot,useDingbats=FALSE)
plot.roc(roc.res,col=cols3[1])
#sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=ROC_colors[i]))
dev.off()

roc.res <- roc(df_merge$resid2, df_merge$R2)
#auc_list <- c(auc_list,pROC::auc(roc.res))
rocplot <- paste("ovary_ROCplot_ordinalForest_2-fold_residual_cancer_prediction2.pdf",sep='')
pdf(rocplot,useDingbats=FALSE)
plot.roc(roc.res,col=cols3[2])
#sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=ROC_colors[i]))
dev.off()


height <- - (df_merge$R0*log2(df_merge$R0) + df_merge$R1*log2(df_merge$R1) + df_merge$R2*log2(df_merge$R2))

df_merge$R0 <- df_merge$R0 * (log2(3) - height)
df_merge$R1 <- df_merge$R1 * (log2(3) - height)
df_merge$R2 <- df_merge$R2 * (log2(3) - height)

df_merge$confidence <- log2(3) - height
df_merge$max_confidence <- rowMax(df_merge[,c("R0","R1","R2")])

df_merge$residual_cancer <- as.factor(df_merge$residual_cancer)

ggplot(data=df_merge,aes(x=hist_type,y=max_confidence,color=hist_type))+
geom_boxplot(color="black",outlier.shape = NA)+
geom_jitter(position=position_jitter(0.2),size=1,shape=20)+
scale_color_manual(values=cols1)+
theme(
	panel.background = element_rect(fill = "white", colour="black"),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.text.x = element_text(angle = 45, hjust = 1),
	legend.position = "none"
)
ggsave(paste("ovary_cancer_residual_prediction_confidence_jitterplot.pdf",sep="_"),width=2.0,height=2.5,useDingbats=FALSE)


df_high_conf <- df_merge[df_merge$max_confidence>=0.2,]
df_low_conf <- df_merge[df_merge$max_confidence<0.20,]

res <- table(df_high_conf$resid,df_high_conf$pred)
res2 <- filltable(res)

colorvec <- setcolor(res2)
blnplot <- paste("ovary_balloonplot_ordinalForest_2-fold_residual_cancer_highconf.pdf",sep='')
pdf(blnplot,useDingbats=F)
balloonplot(t(res2),show.zeros=TRUE,main=NULL,rowmar=margin,text.size=1,label.size=1+t(res2)/max(t(res2)),dotcolor=colorvec,show.margins=TRUE, xlab="prediction",ylab="category",cum.margins=FALSE)
dev.off()
