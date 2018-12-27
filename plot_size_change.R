library(ggplot2)
library(gmodels)

pd <- position_dodge(0.1) 

auc_data <- read.table('cancer_age_171122auc_2-fold_validation_size_change.txt',sep="\t")
acc_data <- read.table('cancer_age_171122accuracy_2-fold_validation_size_change.txt',sep="\t")

auc_ci <- c()
acc_ci <- c()

for (i in c(1:5)){
	auc_ci <- rbind(auc_ci,ci(as.vector(auc_data)[,i]))
	acc_ci <- rbind(acc_ci,ci(as.vector(acc_data)[,i]))
}

auc_ci <- data.frame(Fraction=c(20,40,60,80,100),auc_ci)
acc_ci <- data.frame(Fraction=c(20,40,60,80,100),acc_ci)

colnames(auc_ci) <- c("Fraction","mean","lower","upper","se")
colnames(acc_ci) <- c("Fraction","mean","lower","upper","se")

ggplot(auc_ci, aes(x=Fraction, y=mean),colour="black") + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd)+
   	theme_bw(base_size = 16)+
    theme(
          panel.background = element_rect(size=2,colour="black",fill="white"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.5,linetype = "dashed",colour="black")
    )
ggsave("cancer_age_171122auc_2-fold_validation_size_change.pdf", useDingbats=FALSE, width = 7, height = 4)

ggplot(acc_ci, aes(x=Fraction, y=mean),colour="black") + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.9, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd)+
   	theme_bw(base_size = 16)+
    theme(
          panel.background = element_rect(size=2,colour="black",fill="white"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.5,linetype = "dashed",colour="black")
    )
ggsave("cancer_age_171122accuracy_2-fold_validation_size_change.pdf", useDingbats=FALSE, width = 7, height = 4)
