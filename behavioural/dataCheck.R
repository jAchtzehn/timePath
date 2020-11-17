dsa_j = read.table(file = '/Users/jachtzehn/data/fMRI/timePath/behavioural/RT_data_all_subj.tsv', sep='\t', header=TRUE)
load("/Users/jachtzehn/Downloads/tp2.behav.stana.preproc.RData")
dsa_m = tp2.behav.stana.preproc

# delete all cross tasks from martin's data (already done in johannes' raw data)
dsa_m = subset(dsa_m, task!=4)

# delete all missed responses
dsa_m = subset(dsa_m, event==9)

# delete participant 15 (VP 22) from martin's data (already done in johannes' data)
dsa_m = subset(dsa_m, VP!=22)

# delete all RTs below 300ms from johannes' data
dsa_j = subset(dsa_j, RT >= 0.3)

# delete all NAN from Martin's data (these were all the trials with RT below 300ms)
dsa_m = subset(dsa_m, !is.na(RT))

# plot
#par(mfrow=c(2,1))
#plot(dsa_m$RT~dsa_m$VP)
#plot(dsa_j$RT~dsa_j$subject)


library(nlme)
model_j<-lme(RT~trial_type_nr*run,data=dsa_j,random=~1|subject,method="REML",na.action=na.exclude)
aov_j = anova(model_j)
print(summary(model_j))
model_m<-lme(RT~task*run,data=dsa_m,random=~1|VP,method="REML",na.action=na.exclude)
print(summary(model_m))
aov_m = anova(model_m)