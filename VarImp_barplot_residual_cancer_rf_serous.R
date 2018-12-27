library(ggplot2)
library(reshape2)
library(RColorBrewer)


top_num <- 8

#mycols <- c("#3E3A39","#EB6100","#E60012","#E4007F","#601986","#036EB8","#2EA7E0","#009944")
mycols <- c("black","#777777")

df1 <- read.table('ovary_residual_cancer_prediction1_serouse_rf_importance.txt',header=TRUE,row.names=1)
df2 <- read.table('ovary_residual_cancer_prediction2_serouse_rf_importance.txt',header=TRUE,row.names=1)

df_merge <- data.frame(R0=df1[,1],R1=df2[,1])
rownames(df_merge) <- rownames(df1)

df_top <- df_merge[head(rev(order(rowSums(df_merge))),top_num),]

df_melt <- melt(as.matrix(df_top))

df_melt$Var1 <- factor(df_melt$Var1,levels=rownames(df_top))

ggplot(df_melt,aes(x=Var1,y=value,fill=Var2)) +
geom_bar(width=0.8,position = position_dodge(width = 0.8),stat="identity") +
scale_fill_manual(values=mycols)+
theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 1)
	)

ggsave("ovary_residual_cancer_prediction_importance_rf_serous.pdf",width=3,height=2)