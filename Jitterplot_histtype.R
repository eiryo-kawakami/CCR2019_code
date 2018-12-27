library(ggplot2)
library(reshape2)
library(RColorBrewer)

mycols <-c(brewer.pal(5, "Set1")[c(1,2,3,5)],"#333333")

cutoff <- function(num){
	if(num > 100){
		return(num-240)
	}else return(num)
}

convert_class <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val == "benign"){
			tmp_vec <- c(tmp_vec,NA)
		} else {
			tmp_vec <- c(tmp_vec,val)
		}
	}

	return(tmp_vec)
}

data_file <- "../data/ovary_cancer_update170830_dummyvar.txt"
data_all <- read.table(data_file,sep="\t",header=T,skip=1)
data_all$CA125 <- log10(data_all$CA125)
data_all$CA19.9 <- log10(data_all$CA19.9)
data_all$CRP <- log10(data_all$CRP)
data_all$CEA <- log10(data_all$CEA)


top_num <- 12

df <- read.table('ovary_rf_importance_summary.txt',header=TRUE,row.names=1)

top_vals <- rownames(df[head(rev(order(rowSums(df[,c(3:6)]))),top_num),])

data_all <- na.omit(data.frame(data_all[,top_vals],class=factor(data_all[,"hist_type"],levels=c("HGSC","clearcell","endometrioid","mucinous")),stage=data_all[,"stage"]))
data_all$stage <- factor(convert_class(data_all$stage))
data_all <- na.omit(data_all)

for (val in top_vals){

	df_ind <- data.frame(class=data_all$class,value=data_all[,val])
	if (val == "APTT"){
		print(max(df_ind$value))
		df_ind$value <- sapply(df_ind$value,cutoff)
	}
	ggplot(data=df_ind,aes(x=class,y=value,color=class))+
	geom_boxplot(color="black",outlier.shape = NA)+
	geom_jitter(position=position_jitter(0.2),size=1,shape=20)+
	scale_color_manual(values=mycols)+
	theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1)
	)
	ggsave(paste("ovary_cancer_histtype",val,"jitterplot.pdf",sep="_"),width=2.7,height=2.5,useDingbats=FALSE)
}