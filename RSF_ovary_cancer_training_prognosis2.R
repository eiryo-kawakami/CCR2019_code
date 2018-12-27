################## Load packages ##################
library("ggplot2")         # Graphics engine
library("RColorBrewer")    # Nice color palettes
#library("plot3D")          # for 3d surfaces. 
library("dplyr")           # Better data manipulations
library("parallel")        # mclapply for multicore processing

library("NMF")

# Analysis packages.
library("randomForestSRC") # random forest for survival, regression and 
                           # classification
library("ggRandomForests") # ggplot2 random forest figures (This!)

################ Default Settings ##################
theme_set(theme_bw())     # A ggplot2 theme with white background

## Set open circle for censored, and x for events 
event.marks <- c(1, 4)
event.labels <- c(FALSE, TRUE)

## We want red for death events, so reorder this set.
strCol <- brewer.pal(3, "Set1")[c(2,1,3)]

data_file <- "../data/ovary_cancer_update171205.txt"
#setwd("/Users/kosekikeita/Dropbox/Koseki-san/random_survival_forest/")
data_train <- read.table(data_file,sep="\t",header=T,skip=1)
val_type <- read.table(data_file,sep="\t")[1,]
data_train <- data_train[data_train$stage!="benign",]

#exp_vars <- c("age","CA125","CA19.9","CEA","PT","APTT","Fbg","D.dimer","FDP","WBC","Hb","PLT","AST","ALT","LDH","CHE","Cr","TP","Alb","CRP","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso","Blood_A","Blood_B","Blood_AB","Blood_O")
exp_vars <- c("age","CA125","CA19.9","CEA","D.dimer","FDP","ALT","CRP","PT","APTT","Fbg","WBC","Hb","PLT","AST","LDH","CHE","Cr","TP","Alb","T.Bil","D.Bil","ALP","LAP","gamma.GTP","UN","UA","Na","K","Cl","Hct","Neu","Lym","Mono","Eo","Baso","residual_cancer")
# log_vals <- c("CA125","CA19.9","CEA","D.dimer","FDP","ALT","CHE","CRP")

#data_all <- data_all[data_all$hist_type!="others",]
#data_all <- data_all[data_all$hist_type!="cyst",]
#data_all <- data_all[data_all$hist_type!="teratoma",]

use_exp_vars <- c()

for (exp_var in exp_vars){
 print(exp_var)
 if (length(na.omit(data_train[,exp_var])) >= 300){
   use_exp_vars <- c(use_exp_vars,exp_var)
 }
}

for (i in 1:length(val_type)){
	if (val_type[i] == "log"){
		data_train[,i] <- log(data_train[,i])
	}
}

#data1 <- log(data_all[,log_vals])
data <- na.omit(data_train[,c(use_exp_vars,"time","prog")])
#data <- cbind(data1,data2)

data$prog <- as.logical(data$prog)
data$residual_cancer <- as.factor(data$residual_cancer)
#data$Cr <- log(data$Cr)
#data$LDH <- log(data$LDH)

# df <- data[,use_exp_vars]

# for (i in 1:length(use_exp_vars)){
#     data[,i] <- rank(data[,i])
# }

# mat <- as.matrix(df)

# estim.r <- nmf(mat, 2:25, nrun = 5, seed = 123456)
# pdf("ovary_train_estim.pdf",useDingbats=FALSE)
# plot(estim.r)
# dev.off()

# set.seed(123456)
# res <- nmf(mat, 23)
# pdf("ovary_rank23.pdf",useDingbats=FALSE)
# #layout(cbind(1, 2))
# basismap(res)
# coefmap(res)
# dev.off()

#> head(data)
#  age    CA125   CA19.9        CEA  PT APTT Fbg   D.dimer  WBC   Hb  PLT AST      ALT LDH      CHE   Cr  TP Alb        CRP T.Bil D.Bil ALP gamma.GTP    UN  UA  Na   K  Cl  Hct  Neu  Lym Mono  Eo Baso time  prog
#1  52 3.828641 3.737670  1.9600948  92 34.0 581        NA 6600 13.8 27.9  21 2.197225 444 5.886104 0.84 8.0 4.0 -1.3862944   0.8   0.1 268        13 10.00 5.3 141 4.6 105 43.1 4200 1900  200 200    0  703 FALSE
#2  43 7.243513 6.837333  1.1631508  74 36.8 451        NA 6600 10.3 38.4  14 1.945910 240 5.455321 0.42 7.8 4.2 -0.2107210   0.4   0.0 132        14 10.00 2.6 139 4.2 104 33.1 4700 1300  400 100    0 1019 FALSE
#3  63 3.401197 4.290459  2.0541237  93 26.1 468        NA 9500 13.3 24.2  13 1.945910 153 5.924256 0.55 6.4 3.8 -0.4620355   0.4   0.0 183        19 17.00 6.1 143 4.4 110 38.7 6700 2100  400   0    0  534 FALSE
#4  46 3.737670 3.295837  0.5306283 100 31.1 506 0.2623643 6700 14.1 26.9  23 3.091042 217 5.493061 0.69 7.0 3.9 -1.3862944   0.6   0.0 219        23 10.00 3.3 140 3.9 105 41.4 5100 1100  300   0    0  939 FALSE
#5  68 3.688879 3.295837  0.8329091 100 30.6 616 0.2623643 9600 12.5 32.2  13 2.302585 190 5.826000 0.67 7.4 4.1 -0.2357223   0.5   0.1 281        13  0.67 4.3 141 4.4 105 37.8 6900  200  300 100  100  791 FALSE
#6  39 2.564949 3.496508 -0.3566749  88 33.7 226        NA 5100 12.4 28.1  16 2.833213 134 5.497168 0.57 6.9 4.1 -3.2188758   1.6   0.1 139        16 10.00 4.8 140 4.4 105   NA 3100 1200  100   0  100  759 FALSE

cls <- sapply(data, class) 
labels <- use_exp_vars

dta.labs <- data.frame(cbind(names = colnames(data), label = labels, type = cls))
# Put the "years" variable on top.
dta.labs <- rbind(dta.labs[nrow(dta.labs),], dta.labs[-nrow(dta.labs),])

st.labs <- as.character(dta.labs$label)
names(st.labs) <- rownames(dta.labs)

#The rfsrc function call grows the forest, determining the type of forest by the response supplied in the formula argument.
rfsrc_oc <- rfsrc(Surv(time, prog) ~ ., data = data,
                    ntree = 10000,
                    nsplit = 10, na.action = "na.impute",
                    tree.err = TRUE,importance = TRUE)
#The gg_error function operates on the random forest (rfsrc_pbc) object to extract the error estimates as a function of the number of trees in the forest.
#Figure 6: Random forest OOB prediction error estimates as a function of the number of trees in the forest.
p <- plot(gg_error(rfsrc_oc))
ggsave(file = "error_estimates_training_prognosis.pdf", plot = p)


#Figure 7: Random forest OOB predicted survival. Blue curves correspond to censored observations, red curves correspond to observations experiencing death events.
#ggRFsrc <- plot(gg_rfsrc(rfsrc_oc), alpha = 0.2) +
#  theme(legend.position = "none") +
#  labs(y = "Survival Probability", x = "time(days)") +
#  coord_cartesian(ylim = c(-0.01, 1.01))
#ggsave(file = "Random_forest_OOB_predicted_survival.pdf", plot = ggRFsrc)

#The gg_vimp function extracts VIMP measures for each of the variables used to grow the forest. The plot.gg_vimp function shows the variables, in VIMP rank order, labeled with the named vector in the lbls argument.
#Figure 10: Random forest Variable Importance (VIMP). Blue bars indicates positive VIMP, red indicates negative VIMP. Importance is relative to positive length of bars.
vimp_res <- gg_vimp(rfsrc_oc)
write.table(vimp_res,file="Variable_Importance_training_prognosis.txt",sep="\t")

p <- plot(gg_vimp(rfsrc_oc), lbls = st.labs) +
   theme(legend.position = c(0.8, 0.2)) +
   labs(fill = "VIMP > 0")
ggsave(file = "Variable_Importance_training_prognosis.pdf", plot = p)

#var.select function uses the minimal depth methodology for variable selection, returning an object with both minimal depth and vimp measures.
#gg_minimal_depth function is analogous to the gg_vimp function. Variables are ranked from most important at the top (minimal depth measure), to least at the bottom (maximal minimal depth).
varsel_oc <- var.select(rfsrc_oc)
gg_md <- gg_minimal_depth(varsel_oc, lbls = st.labs)
print(gg_md)
#Figure 11: Minimal Depth variable selection. Low minimal depth indicates important variables. The dashed line is the threshold of maximum value for variable selection.
p <- plot(gg_md, lbls = st.labs)
ggsave(file="Minimal_Depth_variable_selection_training_prognosis.pdf", plot=p)

#Figure 12: Comparing Minimal Depth and Vimp rankings. Points on the red dashed line are ranked equivalently, points above have higher VIMP ranking, those below have higher minimal depth ranking.
p <- plot(gg_minimal_vimp(gg_md), lbls = st.labs) +
   theme(legend.position=c(0.8, 0.2))
ggsave(file="Minimal_Depth_and_Vimp_rankings_training_prognosis.pdf", plot=p)

#Figure 13: Random forest predicted survival (Figure 7) with vertical dashed lines indicate the 1 and 3 year survival estimates.
#ggRFsrc + geom_vline(aes(xintercept = 300), linetype = "dashed") +
#geom_vline(aes(xintercept = 600), linetype = "dashed") +
#coord_cartesian(xlim = c(0, 1000))

#A variable dependence plot is generated from the predicted response value of each survival curve at the intersecting time line plotted against covariate value for that observation.

#gg_v <- gg_variable(rfsrc_oc, time = c(300, 600),
#                   time.labels = c("300 days", "600 days"))

#The gg_variable function extracts the training set variables and the predicted OOB response from rfsrc and predict objects.
#we store the gg_variable data object for later use (gg_v), as all remaining variable dependence plots can be constructed from this object.
#Figure 14: Variable dependence of survival at 1 and 3 years on bili variable. Individual cases are marked with blue circles (alive or censored) and red x (dead). Loess smooth curve with shaded 95% confidence band indicates decreasing survival with increasing bilirubin.
#p <- plot(gg_v, xvar = "TP", alpha = 0.4) + #, se=FALSE
#  labs(y = "Survival", x = st.labs["TP"]) +
#  theme(legend.position = "none") +
#  scale_color_manual(values = strCol, labels = event.labels) +
#  scale_shape_manual(values = event.marks, labels = event.labels) +
#  coord_cartesian(ylim = c(-0.01, 1.01))
#ggsave(file = "TP_dependence.pdf", plot = p)

#p <- plot(gg_v, xvar = "CRP", alpha = 0.4) + #, se=FALSE
#  labs(y = "Survival", x = st.labs["CRP"]) +
#  theme(legend.position = "none") +
#  scale_color_manual(values = strCol, labels = event.labels) +
#  scale_shape_manual(values = event.marks, labels = event.labels) +
#  coord_cartesian(ylim = c(-0.01, 1.01))
#ggsave(file = "CRP_dependence.pdf", plot = p)

#p <- plot(gg_v, xvar = "Alb", alpha = 0.4) + #, se=FALSE
#  labs(y = "Survival", x = st.labs["Alb"]) +
#  theme(legend.position = "none") +
#  scale_color_manual(values = strCol, labels = event.labels) +
#  scale_shape_manual(values = event.marks, labels = event.labels) +
#  coord_cartesian(ylim = c(-0.01, 1.01))
#ggsave(file = "Alb_dependence.pdf", plot = p)

#p <- plot(gg_v, xvar = "LDH", alpha = 0.4) + #, se=FALSE
#  labs(y = "Survival", x = st.labs["LDH"]) +
#  theme(legend.position = "none") +
#  scale_color_manual(values = strCol, labels = event.labels) +
#  scale_shape_manual(values = event.marks, labels = event.labels) +
#  coord_cartesian(ylim = c(-0.01, 1.01))
#ggsave(file = "LDH_dependence.pdf", plot = p)

#p <- plot(gg_v, xvar = "Cr", alpha = 0.4) + #, se=FALSE
#  labs(y = "Survival", x = st.labs["Cr"]) +
#  theme(legend.position = "none") +
#  scale_color_manual(values = strCol, labels = event.labels) +
#  scale_shape_manual(values = event.marks, labels = event.labels) +
#  coord_cartesian(ylim = c(-0.01, 1.01))
#ggsave(file = "Cr_dependence.pdf", plot = p)



#Figure 5: Kaplan–Meier survival estimates comparing different groups of Bilirubin measures (bili) for the pbc data set.
#data.TP <- data
#data.TP$TP_grp <- cut(data.TP$TP, breaks = c(0, 7,7.7, 10))
#p<-plot(gg_survival(interval = "time", censor = "prog", by = "TP_grp",
#               data = data.TP), error = "none") +
#labs(y = "Survival Probability", x = "Observation Time (days)",
#     color = "TP")
#ggsave(file = "TP_Kaplan_Meier.pdf", plot = p)


