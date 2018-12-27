library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(matrixStats)


hist_list = c("serous","clear","endometrioid","mucinous")

df_median <- c()

for (hist in hist_list){
	df_ind <- read.table(paste('ovary',hist,'_rfimportance.txt',sep=""),header=TRUE,row.names=1)
	df_median <- rbind(df_median,rowMedians(as.matrix(df_ind)))
}

rownames(df_median) <- hist_list
colnames(df_median) <- use_exp_vars

top_num <- 33

df_top <- df_median[,head(rev(order(colSums(df_median))),top_num)]
df_top_m <- melt(as.matrix(df_top))

#mycols <- c("#3E3A39","#EB6100","#E60012","#E4007F","#601986","#036EB8","#2EA7E0","#009944")
mycols <-c(brewer.pal(5, "Set1")[c(1,2,3,5)],"#333333")

# df_all <- c()

# for (hist in hist_list){
# 	df_ind <- read.table(paste('ovary',hist,'_rfimportance.txt',sep=""),header=TRUE,row.names=1)
# 	df_ind <- df_ind[as.vector(colnames(df_top)),]
# 	df_ind_m <- melt(as.matrix(df_ind))
# 	df_ind_m <- data.frame(df_ind_m[,c(1,3)],hist=rep(hist,nrow(df_ind_m)))
# 	df_all <- rbind(df_all,df_ind_m)
# }

ggplot(df_top_m,aes(x=Var2,y=value,fill=Var1)) +
geom_bar(position="dodge",stat="identity") +
scale_fill_manual(values=mycols)+
theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		axis.text.x = element_text(size=10, angle = 90, hjust = 1),
		panel.grid.minor=element_blank()
	)

ggsave("ovary_histtype_importance_rep.pdf",width=10,height=3,useDingbats=FALSE)
