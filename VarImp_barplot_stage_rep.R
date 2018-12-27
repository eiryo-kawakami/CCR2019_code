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

df_median <- c()

df_ind <- read.table('ovarystage_rfimportance.txt',header=TRUE,row.names=1)
df_median <- rowMedians(as.matrix(df_ind))

top_num <- 33

df_top <- df_ind[head(rev(order(df_median)),top_num),]
df_top_m <- melt(as.matrix(df_top))

#mycols <- c("#3E3A39","#EB6100","#E60012","#E4007F","#601986","#036EB8","#2EA7E0","#009944")
mycols <- c("#3E3A39",brewer.pal(7,"Paired"))

ggplot(df_top_m,aes(x=Var1,y=value,fill=Var1)) +
geom_boxplot() +
#scale_fill_manual(values=mycols)+
theme(
	panel.background = element_rect(fill = "white", colour="black"),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.text.x = element_text(size=10, angle = 90, hjust = 1),
	legend.position = "none"
)

ggsave("ovary_stage_importance_rep.pdf",width=10,height=3,useDingbats=FALSE)
