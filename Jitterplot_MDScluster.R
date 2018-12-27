library(ggplot2)
library(reshape2)
library(RColorBrewer)

cols <-brewer.pal(5, "Set1")

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
data_all <- read.table(data_file,sep="\t",header=T,skip=1,row.names=1)
#data_all <- data_all[data_all$FIGO=="III",]
data_all$CA125 <- log10(data_all$CA125)
data_all$CA19.9 <- log10(data_all$CA19.9)
data_all$CRP <- log10(data_all$CRP)
data_all$CEA <- log10(data_all$CEA)


exp_vars <- c("age","CA125","CA19.9","CEA","D.dimer","FDP","ALT","CRP","PT","APTT","Fbg","WBC","Hb","PLT","AST","LDH","CHE","Cr","TP","Alb","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso","residual_cancer")

use_exp_vars <- c()

for (exp_var in exp_vars){
 print(exp_var)
 if (length(na.omit(data_all[,exp_var])) >= 400){
   use_exp_vars <- c(use_exp_vars,exp_var)
 }
}

data <- na.omit(data_all[,c(use_exp_vars)])

MDSdata <- read.table('MDSdata_cluster.txt')

MDSdata_merged <- data.frame(MDSdata,data[rownames(MDSdata),])
MDSdata_merged$clust <- factor(MDSdata_merged$clust,levels=c(2,1))

anova_pval <- c()

for (val in use_exp_vars){
	df_ind <- na.omit(data.frame(cluster=MDSdata_merged[MDSdata_merged$stage=="early",]$clust,value=MDSdata_merged[MDSdata_merged$stage=="early",val]))
	ggplot(data=df_ind,aes(x=cluster,y=value))+
	geom_boxplot(color="black",outlier.shape = NA)+
	geom_jitter(aes(colour=cluster),position=position_jitter(0.2),size=1,shape=20)+
	scale_color_manual(values=cols[c(2,1)])+
	theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank()
	)
	ggsave(paste("ovary_cancer_cluster",val,"jitterplot.pdf",sep="_"),width=2.5,height=2.5,useDingbats=FALSE)

	if (val != "hist_type"){
		est <- lm(value~cluster,data=df_ind)
		anova_pval <- c(anova_pval,coef(summary(est))[2,4])
	} else {
		anova_pval <- c(anova_pval,NA)
	}
}


df_anova <- data.frame(clinical_items=use_exp_vars,anova_pval=anova_pval,anova_FDR=p.adjust(anova_pval))

write.table(df_anova,file="ovary_cancer_cluster_lm.txt",row.names=F,sep="\t")

