library(corrplot)
library(RColorBrewer)

data_file <- "../data/ovary_cancer_update170830_dummyvar.txt"
data_all <- read.table(data_file,sep="\t",header=T,skip=1)

top_num <- 8

df_imp <- read.table('ovary_cancer_age_171122_importance_summary.txt',header=TRUE,row.names=1)

top_vals <- rownames(df_imp[head(rev(order(rowSums(df_imp))),top_num),])

data <- na.omit(data_all[,c(top_vals,"stage")])
df <- data[,top_vals]

cor.res <- cor(df,method="spearman")

pdf("ovary_cancer_correlation.pdf",width=4,height=4,useDingbats=FALSE)
corrplot(cor.res,type="upper",order="hclust",col=rev(brewer.pal(n=8,name="RdYlBu")))
dev.off()