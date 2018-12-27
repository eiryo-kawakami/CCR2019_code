library(ggplot2)
library(reshape2)
library(RColorBrewer)

mycols <-brewer.pal(4, "Spectral")

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

top_num <- 12

df <- read.table('ovary_rf_importance_summary.txt',header=TRUE,row.names=1)

top_vals <- rownames(df[head(rev(order(rowSums(df[,"stage",drop=F]))),top_num),])

data_all <- na.omit(data.frame(data_all[,top_vals],class=data_all[,"FIGO"]))
data_all$class <- factor(convert_class(data_all$class))
data_all <- na.omit(data_all)

for (val in top_vals){
	df_ind <- data.frame(class=data_all$class,value=data_all[,val])
	ggplot(data=df_ind,aes(x=class,y=value,color=class))+
	geom_boxplot(color="black",outlier.shape = NA)+
	geom_jitter(position=position_jitter(0.2),size=1,shape=20)+
	scale_color_manual(values=mycols)+
	theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank()
	)
	ggsave(paste("ovary_cancer_stage",val,"jitterplot.pdf",sep="_"),width=2.5,height=2.5,useDingbats=FALSE)
}