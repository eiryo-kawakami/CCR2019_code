library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pROC)

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


cols1 <- c(brewer.pal(5, "Set1")[c(1,2,3,5)],"#555555")
cols2 <- gg_color_hue(3)[c(3,2,1)]
cols3 <- c("#999999","#E69F00","#56B4E9")

df1 <- read.table('ovary_residual_cancer_test_prediction1_rf.txt',sep="\t",header=T,row.names=1)
df2 <- read.table('ovary_residual_cancer_test_prediction2_rf.txt',sep="\t",header=T,row.names=1)

data_train <- read.table("../data/ovary_cancer_update171205_train2.txt",sep="\t",header=T)
data_test <- read.table("../data/ovary_cancer_update171205_test2.txt",sep="\t",header=T)

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

df_merge <- data.frame(R0=df1[,1],R1=df1[,2]-df2[,2],R2=df2[,2],data_test_merged[as.vector(rownames(df1)),c("FIGO","hist_type","residual_cancer")])
df_merge$hist_type <- factor(df_merge$hist_type,levels=c("serous","clearcell","endometrioid","mucinous","others"))

df_merge <- df_merge[order(df_merge$hist_type),]
df_merge <- df_merge[order(df_merge$FIGO),]
df_max <- is_max(df_merge[,c(1:2)])
colnames(df_max) <- c("R0_is_max","R1_is_max")
df_merge <- data.frame(df_merge,df_max)

df_merge$resid1 <- factor(convert_residual_class1(df_merge$residual_cancer),levels=c("R0","R1"))
#df_merge$resid2 <- factor(convert_residual_class2(df_merge$residual_cancer),levels=c("R0","R2"))

# roc.res <- roc(df_merge$resid1, df_merge$R1+df_merge$R2)
# #auc_list <- c(auc_list,pROC::auc(roc.res))
# rocplot <- paste("ovary_ROCplot_ordinalForest_2-fold_residual_cancer_prediction1.pdf",sep='')
# pdf(rocplot,useDingbats=FALSE)
# plot.roc(roc.res,col=cols3[1])
# #sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=ROC_colors[i]))
# dev.off()

# roc.res <- roc(df_merge$resid2, df_merge$R2)
# #auc_list <- c(auc_list,pROC::auc(roc.res))
# rocplot <- paste("ovary_ROCplot_ordinalForest_2-fold_residual_cancer_prediction2.pdf",sep='')
# pdf(rocplot,useDingbats=FALSE)
# plot.roc(roc.res,col=cols3[2])
# #sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=ROC_colors[i]))
# dev.off()


height <- - (df_merge$R0*log2(df_merge$R0) + df_merge$R1*log2(df_merge$R1))

df_merge$R0 <- df_merge$R0 * (log2(2) - height)
df_merge$R1 <- df_merge$R1 * (log2(2) - height)
#df_merge$R2 <- df_merge$R2 * (log2(3) - height)

data <- melt(as.matrix(df_merge[,c(1:2)]))
data$Var1 <- factor(data$Var1,levels=rownames(df_merge))

pdf <- "ovary_residual_cancer_test_prediction_rf.pdf"

p <- ggplot(df_merge)+
geom_line(aes(x=(1:nrow(df_merge)),y=R0),colour=cols3[1])+
geom_line(aes(x=(1:nrow(df_merge)),y=R1),colour=cols3[2])+
#geom_line(aes(x=(1:nrow(df_merge)),y=R2),colour=cols3[3])+
geom_point(aes(x=(1:nrow(df_merge)),y=R0),shape=16,alpha=0.7,colour=cols3[1])+
geom_point(aes(x=(1:nrow(df_merge)),y=R1),shape=16,alpha=0.7,colour=cols3[2])+
#geom_point(aes(x=(1:nrow(df_merge)),y=R2),shape=16,alpha=0.7,colour=cols3[3])+
#scale_fill_gradient(low="white",high="red")+
theme(
	panel.background = element_rect(fill = "white", colour="black"),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = pdf,plot=p,width=10,height=1)

data_FIGO <- melt(as.matrix(df_merge[,"FIGO",drop=FALSE]))
data_FIGO$Var1 <- factor(data_FIGO$Var1,levels=rownames(df_merge))

pdf <- "ovary_residual_cancer_test_prediction_rf_FIGO.pdf"

p <- ggplot(data_FIGO) +
geom_tile(aes(x=as.factor(Var1),y=1,fill=as.factor(value)))+
scale_fill_brewer(palette = "Set1") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = pdf,plot=p)


data_histtype <- melt(as.matrix(df_merge[,"hist_type",drop=FALSE]))
data_histtype$Var1 <- factor(data_histtype$Var1,levels=rownames(df_merge))
data_histtype$value <- factor(data_histtype$value,levels=c("serous","clearcell","endometrioid","mucinous","others"))

pdf <- "ovary_residual_cancer_test_prediction_rf_hist_type.pdf"

p <- ggplot(data_histtype) +
geom_tile(aes(x=as.factor(Var1),y=1,fill=as.factor(value)))+
scale_fill_manual(values = cols1) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = pdf,plot=p)


data_residual_cancer <- melt(as.matrix(df_merge[,"residual_cancer",drop=FALSE]))
data_residual_cancer$Var1 <- factor(data_residual_cancer$Var1,levels=rownames(df_merge))

pdf <- "ovary_residual_cancer_test_prediction_rf_residual_cancer.pdf"

p <- ggplot(data_residual_cancer) +
geom_tile(aes(x=as.factor(Var1),y=1,fill=as.factor(value)))+
scale_fill_manual(values = cols3) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = pdf,plot=p)


