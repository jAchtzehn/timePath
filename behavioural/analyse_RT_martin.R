load("/Users/jachtzehn/Downloads/tp2.behav.stana.preproc.RData")
dsa_m <- tp2.behav.stana.preproc

### visual inspection of RTs
par(mfrow=c(2,1))
plot(dsa_m$RT~dsa_m$VP)

# delete trials in which subject confused the button (because those were unrepresentatively easy)
dsa_m$RT[dsa_m$VP==11 & dsa_m$task==3 & dsa_m$sess==1] <- NA
dsa_m$RT[dsa_m$VP==12 & dsa_m$task==3 & dsa_m$sess==1] <- NA
dsa_m$RT[dsa_m$VP==13 & dsa_m$task==3 & dsa_m$sess==1] <- NA

### aggregation
a_m <- aggregate(dsa_m$RT, by=list(dsa_m$task,dsa_m$run,dsa_m$session,dsa_m$VP), mean, na.rm=T)
names(a_m)[names(a_m)=="Group.1"]<-"task"
names(a_m)[names(a_m)=="Group.2"]<-"run"
names(a_m)[names(a_m)=="Group.3"]<-"sess"
names(a_m)[names(a_m)=="Group.4"]<-"VP"
names(a_m)[names(a_m)=="x"]<-"RT"
dsa_m <- subset(a_m, task!=4)

# remove vp22
dsa_m <- subset(dsa_m, VP!=22)

library(nlme)
model<-lme(RT~task*run,data=dsa_m,random=~1|VP,method="REML",na.action=na.exclude)
aov_m<-anova(model)
