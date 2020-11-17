dsa = read.table(file = '/Users/jachtzehn/data/fMRI/timePath/behavioural/vas_data_all_subj.tsv', sep='\t', header=TRUE)

# data clean up steps
# delete trials in which subject confused the button (because those were unrepresentatively easy)
dsa$vas[dsa$subject==4 & dsa$trial_type=='lumin' & dsa$session==1] = NA
dsa$vas[dsa$subject==5 & dsa$trial_type=='lumin' & dsa$session==1] = NA
dsa$vas[dsa$subject==6 & dsa$trial_type=='lumin' & dsa$session==1] = NA
dsa = subset(dsa, !is.na(vas))

#change values below 0 to 0 and above 1 to 1
dsa$vas[dsa$vas < 0] = 0
dsa$vas[dsa$vas > 1] = 1

# remove participant 15, as he missed more than 20% of responses
dsa = subset(dsa, subject!=15)

par(mfrow=c(1,1))
plot(dsa$vas~dsa$run)
#plot(dsa$vas~dsa$trial_type)

library(nlme)
lmemodel =lme(vas~trial_type*run,data=dsa,random=~1|subject,method="REML",na.action=na.exclude)
aov_lmemodel = anova(lmemodel)
#print(summary(lmemodel))
print(aov_lmemodel)

library(multcompView)
lmodel = lm('vas ~ trial_type', data=dsa)
aov_lmodel = aov(lmodel)

tukey_model = TukeyHSD(x=aov_lmodel)
print(tukey_model)