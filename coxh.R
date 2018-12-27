library(survival)

convert_age_class <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (val >= 50){
			tmp_vec <- c(tmp_vec,"over50")
		} else {
			tmp_vec <- c(tmp_vec,"under50")
		}
	}
	return(tmp_vec)
}

convert_histotype_class <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (is.na(val)){
			tmp_vec <- c(tmp_vec,NA)
		} else if (val == "HGSC"){
			tmp_vec <- c(tmp_vec,"HGSC")
		} else {
			tmp_vec <- c(tmp_vec,"others")
		}
	}
	return(tmp_vec)
}

convert_residual_class1 <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (is.na(val)){
			tmp_vec <- c(tmp_vec,NA)
		} else if (val == 0){
			tmp_vec <- c(tmp_vec,"R0")
		} else {
			tmp_vec <- c(tmp_vec,"R1")
		}
	}
	return(tmp_vec)
}

convert_residual_class2 <- function(vec){
	tmp_vec <- c()
	for (val in vec){
		if (is.na(val)){
			tmp_vec <- c(tmp_vec,NA)
		} else if (val == 2){
			tmp_vec <- c(tmp_vec,"R2")
		} else {
			tmp_vec <- c(tmp_vec,"R0")
		}
	}
	return(tmp_vec)
}

MDSdata <- read.table('MDSdata_cluster.txt',sep="\t")
MDSdata$clust <- factor(MDSdata$clust,levels=c("2","1"))
data_file <- "../data/ovary_cancer_update171205_mod.txt"
data_all <- read.table(data_file,sep="\t",header=T,skip=1,row.names=1)

MDSdata_merged <- data.frame(MDSdata,age=data_all[as.vector(rownames(MDSdata)),"age"])
MDSdata_merged$age <- factor(convert_age_class(MDSdata_merged$age),levels=c("over50","under50"))
MDSdata_merged$stage <- factor(MDSdata_merged$stage,levels=c("early","late"))
data_all$stage <- factor(data_all$stage,levels=c("early","late"))
MDSdata_merged$hist_type <- factor(convert_histotype_class(MDSdata_merged$hist_type),levels=c("others","HGSC"))
data_all$hist_type <- factor(convert_histotype_class(data_all$hist_type),levels=c("others","HGSC"))

#MDSdata_merged <- na.omit(MDSdata_merged)

MDSdata_merged$resid1 <- factor(convert_residual_class1(MDSdata_merged$residual_cancer),levels=c("R0","R1"))
MDSdata_merged$resid2 <- factor(convert_residual_class2(MDSdata_merged$residual_cancer),levels=c("R0","R2"))

data_all$resid1 <- factor(convert_residual_class1(data_all$residual_cancer),levels=c("R0","R1"))
data_all$resid2 <- factor(convert_residual_class2(data_all$residual_cancer),levels=c("R0","R2"))

res <- c()

data_tmp <- na.omit(data_all[data_all$stage=="early" | data_all$stage=="late",c("time_recurrence","recurrence","age")])

est <- coxph(Surv(time_recurrence,recurrence)~age,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="early" | data_all$stage=="late",c("time_recurrence","recurrence","hist_type")])

est <- coxph(Surv(time_recurrence,recurrence)~hist_type,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="early" | data_all$stage=="late",c("time_recurrence","recurrence","stage")])

est <- coxph(Surv(time_recurrence,recurrence)~stage,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="early" | data_all$stage=="late",c("time_recurrence","recurrence","resid1")])

est <- coxph(Surv(time_recurrence,recurrence)~resid1,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="early" | data_all$stage=="late",c("time_recurrence","recurrence","resid2")])

est <- coxph(Surv(time_recurrence,recurrence)~resid2,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(MDSdata_merged[MDSdata_merged$stage=="early" | MDSdata_merged$stage=="late",c("time_recurrence","recurrence","clust")])

est <- coxph(Surv(time_recurrence,recurrence)~clust,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))
rownames(res) <- c("age","HGSC","stage","residual_1","residual_2","cluster")

write.table(res,file="ovary_cancer_coxph_all_stat.txt",sep="\t",quote=F)


res <- c()

data_tmp <- na.omit(data_all[data_all$stage=="early",c("time_recurrence","recurrence","age")])

est <- coxph(Surv(time_recurrence,recurrence)~age,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="early",c("time_recurrence","recurrence","hist_type")])

est <- coxph(Surv(time_recurrence,recurrence)~hist_type,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(MDSdata_merged[MDSdata_merged$stage=="early",c("time_recurrence","recurrence","clust")])

est <- coxph(Surv(time_recurrence,recurrence)~clust,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))
rownames(res) <- c("age","HGSC","cluster")

write.table(res,file="ovary_cancer_coxph_early_stat.txt",sep="\t",quote=F)


res <- c()

data_tmp <- na.omit(data_all[data_all$stage=="late",c("time_recurrence","recurrence","age")])

est <- coxph(Surv(time_recurrence,recurrence)~age,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="late",c("time_recurrence","recurrence","hist_type")])

est <- coxph(Surv(time_recurrence,recurrence)~hist_type,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="late",c("time_recurrence","recurrence","resid1")])
#data_tmp$FIGO <- factor(data_tmp$FIGO,levels=c("III","IV"))

est <- coxph(Surv(time_recurrence,recurrence)~resid1,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(data_all[data_all$stage=="late",c("time_recurrence","recurrence","resid2")])

est <- coxph(Surv(time_recurrence,recurrence)~resid2,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){ 
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))

data_tmp <- na.omit(MDSdata_merged[MDSdata_merged$stage=="late",c("time_recurrence","recurrence","clust")])

est <- coxph(Surv(time_recurrence,recurrence)~clust,method="exact",data=data_tmp)

calc_coxph_stat <- function(x){
	x <- summary(x)
	p.value<-signif(x$wald["pvalue"], digits=3)
	log.test<-signif(x$logtest["test"], digits=3)
	beta<-signif(x$coef[1], digits=2);#coeficient beta
	HR <-signif(x$coef[2], digits=2);#exp(beta)
	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
	HR <- paste0(HR, " (",
	HR.confint.lower, "-", HR.confint.upper, ")")
	res<-c(beta, HR, log.test, p.value)
	names(res)<-c("beta", "HR (95% CI for HR)", "likelifood-ratio.test", "p.value")
	return(res)
	#return(exp(cbind(coef(x),confint(x))))
}

res <- rbind(res,calc_coxph_stat(est))
rownames(res) <- c("age","HGSC","residual_1","residual_2","cluster")

write.table(res,file="ovary_cancer_coxph_late_stat.txt",sep="\t",quote=F)
