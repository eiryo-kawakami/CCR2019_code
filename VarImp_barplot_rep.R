library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(matrixStats)

data_train <- read.table("../data/ovary_cancer_update170830_dummyvar_train.txt",sep="\t",header=T)
data_test <- read.table("../data/ovary_cancer_update170830_dummyvar_test.txt",sep="\t",header=T)

#exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso","Blood_A","Blood_B","Blood_AB","Blood_O")
exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso")

use_exp_vars <- c()

for (exp_var in exp_vars){
	print(exp_var)
	if (length(na.omit(data_train[,exp_var])) >= 170){
		use_exp_vars <- c(use_exp_vars,exp_var)
	}
}

# method_list = c("glm","gbm","svmRadial","rf","cforest","nb","nnet","glmnet")
method_list = c("glm","rf")

df_median <- c()

for (method in method_list){
	df_ind <- read.table(paste('ovarycancer_age_171122_',method,'importance.txt',sep=""),header=TRUE,row.names=1)
	df_median <- rbind(df_median,rowMedians(as.matrix(df_ind)))
}

rownames(df_median) <- method_list
colnames(df_median) <- use_exp_vars

top_num <- 33

df_top <- df_median[,head(rev(order(colSums(df_median["rf",,drop=F]))),top_num)]
df_top_m <- melt(as.matrix(df_top))

#mycols <- c("#3E3A39","#EB6100","#E60012","#E4007F","#601986","#036EB8","#2EA7E0","#009944")
# mycols <- c("#3E3A39",brewer.pal(7,"Paired"))
mycols <- c("#3E3A39","#009944")

df_all <- c()

for (method in method_list){
	df_ind <- read.table(paste('ovarycancer_age_171122_',method,'importance.txt',sep=""),header=TRUE,row.names=1)
	df_ind <- df_ind[as.vector(colnames(df_top)),]
	df_ind_m <- melt(as.matrix(df_ind))
	df_ind_m <- data.frame(df_ind_m[,c(1,3)],method=rep(method,nrow(df_ind_m)))
	df_all <- rbind(df_all,df_ind_m)
}

ggplot(df_top_m,aes(x=Var2,y=value,fill=Var1)) +
geom_bar(position="dodge",stat="identity") +
scale_fill_manual(values=mycols)+
theme(
		panel.background = element_rect(fill = "white", colour="black"),
		panel.grid.major=element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 1),
		panel.grid.minor=element_blank()
	)

ggsave("ovary_cancer_age_171122_importance_rep_v2.pdf",width=10,height=2.5)
