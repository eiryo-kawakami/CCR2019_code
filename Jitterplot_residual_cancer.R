library(ggplot2)
library(reshape2)
library(RColorBrewer)

mycols <- c("#3D3939","#D5221E")

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

data_file <- "../data/ovary_cancer_update171205.txt"
data_all <- read.table(data_file,sep="\t",header=T,skip=1)
data_all <- data_all[data_all$stage!="benign",]
data_all <- data_all[data_all$FIGO!="I",]
data_all$CA125 <- log10(data_all$CA125)
data_all$CA19.9 <- log10(data_all$CA19.9)
data_all$CRP <- log10(data_all$CRP)
data_all$CEA <- log10(data_all$CEA)

top_num <- 8


df1 <- read.table('ovary_residual_cancer_prediction1_rf_importance.txt',header=TRUE,row.names=1)
df2 <- read.table('ovary_residual_cancer_prediction2_rf_importance.txt',header=TRUE,row.names=1)

df_merge <- data.frame(R0=df1[,1],R1=df2[,1])
rownames(df_merge) <- rownames(df1)

top_vals <- rownames(df_merge[head(rev(order(rowSums(df_merge))),top_num),])

data_all <- data.frame(data_all[,top_vals],residual_cancer=as.factor(data_all[,"residual_cancer"]))
#data_all$class <- factor(convert_class(data_all$class))
#data_all <- na.omit(data_all)

anova_pval <- c()

for (val in top_vals){
	df_ind <- na.omit(data.frame(residual_cancer=data_all$residual_cancer,value=data_all[,val]))
	ggplot(data=df_ind,aes(x=residual_cancer,y=value))+
	geom_boxplot(color="black",outlier.shape = NA)+
	geom_jitter(position=position_jitter(0.2),size=1,shape=20)+
	#scale_color_manual(values=mycols)+
	theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank()
	)
	ggsave(paste("ovary_cancer_residual_cancer",val,"jitterplot.pdf",sep="_"),width=2,height=2.5,useDingbats=FALSE)
}
