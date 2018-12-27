library(ggplot2)
library(reshape2)
library(RColorBrewer)


top_num <- 8

#mycols <- c("#3E3A39","#EB6100","#E60012","#E4007F","#601986","#036EB8","#2EA7E0","#009944")
mycols <- c("#3E3A39",brewer.pal(7,"Paired"))

df <- read.table('ovary_rf_importance_summary.txt',header=TRUE,row.names=1)

df_top <- df[head(rev(order(rowSums(df[,"stage",drop=F]))),top_num),]

data <- data.frame(Var=factor(rownames(df_top),levels=rownames(df_top)),df_top)

ggplot(data,aes(x=Var,y=stage)) +
geom_bar(width=0.8,position = position_dodge(width = 0.8),stat="identity") +
scale_fill_manual(values=mycols)+
theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 1)
	)

ggsave("ovary_stage_importance_rf.pdf",width=3,height=2)