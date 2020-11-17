dsa = read.table(file = '/Users/jachtzehn/data/fMRI/timePath/behavioural/RT_data_all_subj.tsv', sep='\t', header=TRUE)

# data clean up steps
# delete trials in which subject confused the button (because those were unrepresentatively easy)
dsa$RT[dsa$subject==4 & dsa$trial_type=='lumin' & dsa$session==1] = NA
dsa$RT[dsa$subject==5 & dsa$trial_type=='lumin' & dsa$session==1] = NA
dsa$RT[dsa$subject==6 & dsa$trial_type=='lumin' & dsa$session==1] = NA
dsa = subset(dsa, !is.na(RT))

# remove participant 15, as he missed more than 20% of responses
dsa = subset(dsa, subject!=15)

# remove RT that are quicker than 300ms
dsa = subset(dsa, RT > 0.3)

### aggregation
a = aggregate(dsa$RT, by=list(dsa$trial_type,dsa$run,dsa$session,dsa$subject), mean, na.rm=T)
names(a)[names(a)=="Group.1"]<-"trial_type"
names(a)[names(a)=="Group.2"]<-"run"
names(a)[names(a)=="Group.3"]<-"session"
names(a)[names(a)=="Group.4"]<-"subject"
names(a)[names(a)=="x"]<-"RT"

### visual inspection of RTs
#par(mfrow=c(2,1))
#plot(a$RT~a$run)
#hist(a$RT)

# linear mixed effects model
library(nlme)
lmemodel =lme(RT~trial_type*run,data=dsa,random=~1|subject,method="REML",na.action=na.exclude)
print(summary(lmemodel))
aov_lmemodel = anova(lmemodel)
print(aov_lmemodel)
# anova
library(multcompView)
lmodel = lm('RT ~ trial_type', data=dsa)
aov_lmodel = aov(lmodel)

tukey_model = TukeyHSD(x=aov_lmodel)
print(tukey_model)