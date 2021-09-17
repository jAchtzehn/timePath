###################################################################################
### define directory
###################################################################################

#setwd("...")   # insert directory
setwd("/Volumes/Elements/Martin/projects/timepath/tp2/paper/3_Neuroimage/scripts")   # insert directory






#####################################################################################
### Analysis of perceived difficulty

load("difficulty_ratings.RData")
dsa <- difficulty_ratings

# definition of column names:
# 
# subj = subject number
# session = session number (first or second scanning session)
# run = run number (runs 1-4 belong to session 1, run 5-8 belong to session 2)
# task = task type (0: time; 1: space; 2: numerosity; 3: luminance)
# RT = reaction time
# resp = response (from 0 ["easy"] to 1 ["difficult"])

dsa$task <- as.factor(dsa$task)

a1 <- aov(dsa$resp ~ dsa$task * dsa$run)
anova(a1)
summary(a1)
# calculate partial_eta_square (according to Lakens (2013) Front Psychol)
set = anova(a1)
(set$`F value`*set$Df)/(set$`F value`*set$Df+748)
posthoc <- TukeyHSD(x=a1, 'dsa$task', conf.level=0.95)
posthoc

# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
# for time vs luminance:
mdiff = mean(dsa$resp[dsa$task==0]-dsa$resp[dsa$task==3], na.rm=T)
sd1 = sd(dsa$resp[dsa$task==0], na.rm=T)
sd2 = sd(dsa$resp[dsa$task==3], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)
# for space vs luminance:
mdiff = mean(dsa$resp[dsa$task==1]-dsa$resp[dsa$task==3], na.rm=T)
sd1 = sd(dsa$resp[dsa$task==1], na.rm=T)
sd2 = sd(dsa$resp[dsa$task==3], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)
# for numerosity vs luminance:
mdiff = mean(dsa$resp[dsa$task==2]-dsa$resp[dsa$task==3], na.rm=T)
sd1 = sd(dsa$resp[dsa$task==2], na.rm=T)
sd2 = sd(dsa$resp[dsa$task==3], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)
# for time vs space:
mdiff = mean(dsa$resp[dsa$task==0]-dsa$resp[dsa$task==1], na.rm=T)
sd1 = sd(dsa$resp[dsa$task==0], na.rm=T)
sd2 = sd(dsa$resp[dsa$task==1], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)
# for time vs numerosity:
mdiff = mean(dsa$resp[dsa$task==0]-dsa$resp[dsa$task==2], na.rm=T)
sd1 = sd(dsa$resp[dsa$task==0], na.rm=T)
sd2 = sd(dsa$resp[dsa$task==2], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)
# for space vs numerosity:
mdiff = mean(dsa$resp[dsa$task==1]-dsa$resp[dsa$task==2], na.rm=T)
sd1 = sd(dsa$resp[dsa$task==1], na.rm=T)
sd2 = sd(dsa$resp[dsa$task==2], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)











#####################################################################################
### Analysis of reaction times

load("reaction_times.RData")
dsa <- reaction_times

# definition of column names:
# 
# subj = subject number
# session = session number (first or second scanning session)
# run = run number (runs 1-4 belong to session 1, run 5-8 belong to session 2)
# trial = trial number
# task = task type (0: time; 1: space; 2: numerosity; 3: luminance)
# timeStand = travel time (either 2.8s or 4.8s)
# distStand = traveled distance (either 11.5m or 19.7m)
# dotsStand = number of dots (either 45 or 77)
# luminStand = white content (either 0.16 or 0.28)
# RT = reaction time
# resp = response (either 1 [left key] or 2 [right key])

# delete trials in which subject confused the button (because those were unrepresentatively easy)
dsa$RT[dsa$subj==4 & dsa$task==3 & dsa$session==1] <- NA
dsa$RT[dsa$subj==5 & dsa$task==3 & dsa$session==1] <- NA
dsa$RT[dsa$subj==6 & dsa$task==3 & dsa$session==1] <- NA


### aggregation
a <- aggregate(dsa$RT, by=list(dsa$task,dsa$run,dsa$subj), mean, na.rm=T)
names(a)[names(a)=="Group.1"]<-"task"
names(a)[names(a)=="Group.2"]<-"run"
names(a)[names(a)=="Group.3"]<-"subj"
names(a)[names(a)=="x"]<-"RT"
dsa <- subset(a, task!=4)

# remove subject 15
dsa <- subset(dsa, subj!=15)

dsa$task <- as.factor(dsa$task)

a1 <- aov(dsa$RT ~ dsa$task * dsa$run)
anova(a1)
summary(a1)
# calculate partial_eta_square (according to Lakens (2013) Front Psychol)
set = anova(a1)
(set$`F value`*set$Df)/(set$`F value`*set$Df+748)










###################################################################################
### Analysis of psychometric functions

load("pse_data_all_common.RData")
dsa <- pse_data_all_common
library(coin)



###################################
### Participants were able to differentiate between the small and the large value of each dimension:

# Time:
shortTime <- subset(dsa, task==0 & stand==1)
longTime <- subset(dsa, task==0 & stand==2)
t.test(longTime$thr50, shortTime$thr50, alternative="greater", paired=T, var.equal=F, con.level=0.95)
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(shortTime$thr50-longTime$thr50, na.rm=T)
abs(mdiff)/((sd(shortTime$thr50, na.rm=T)+sd(longTime$thr50, na.rm=T))/2)

# Space:
shortDistance <- subset(dsa, task==1 & stand==1)
longDistance <- subset(dsa, task==1 & stand==2)
t.test(longDistance$thr50, shortDistance$thr50, alternative="greater", paired=T, var.equal=F, con.level=0.95)
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(shortDistance$thr50-longDistance$thr50, na.rm=T)
abs(mdiff)/((sd(shortDistance$thr50, na.rm=T)+sd(longDistance$thr50, na.rm=T))/2)

# Numerosity:
smallAmount <- subset(dsa, task==2 & stand==1)
largeAmount <- subset(dsa, task==2 & stand==2)
t.test(largeAmount$thr50, smallAmount$thr50, alternative="greater", paired=T, var.equal=F, con.level=0.95)
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(smallAmount$thr50-largeAmount$thr50, na.rm=T)
abs(mdiff)/((sd(smallAmount$thr50, na.rm=T)+sd(largeAmount$thr50, na.rm=T))/2)

# Luminance:
lowValue <- subset(dsa, task==3 & stand==1)
highValue <- subset(dsa, task==3 & stand==2)
t.test(highValue$thr50, lowValue$thr50, alternative="greater", paired=T, var.equal=F, con.level=0.95)
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(lowValue$thr50-highValue$thr50, na.rm=T)
abs(mdiff)/((sd(lowValue$thr50, na.rm=T)+sd(highValue$thr50, na.rm=T))/2)




###################################
### Precision decreased with increasing magnitude (for time, space and numerosity, but not for luminance)

# Time:
t.test(longTime$width, shortTime$width, alternative="greater", paired=T, var.equal=F, con.level=0.95)
abs(mean(shortTime$width-longTime$width, na.rm=T))/((sd(shortTime$width, na.rm=T)+sd(longTime$width, na.rm=T))/2)

# Space:
t.test(longDistance$width, shortDistance$width, alternative="greater", paired=T, var.equal=F, con.level=0.95)
abs(mean(shortDistance$width-longDistance$width, na.rm=T))/((sd(shortDistance$width, na.rm=T)+sd(longDistance$width, na.rm=T))/2)

# Numerosity:
t.test(largeAmount$width, smallAmount$width, alternative="greater", paired=T, var.equal=F, con.level=0.95)
abs(mean(smallAmount$width-largeAmount$width, na.rm=T))/((sd(smallAmount$width, na.rm=T)+sd(largeAmount$width, na.rm=T))/2)

# Luminance:
t.test(highValue$width, lowValue$width, alternative="greater", paired=T, var.equal=F, con.level=0.95)
abs(mean(lowValue$width-highValue$width, na.rm=T))/((sd(lowValue$width, na.rm=T)+sd(highValue$width, na.rm=T))/2)




###################################################################################
### Correlation of precision between tasks (for each task, small and large standards are standardised and analysed together)

load("pse_data_all_sep_std_value.RData")
dsa <- pse_data_all_sep_std_value


# between time and space
cor.test(dsa$width[dsa$task==0], dsa$width[dsa$task==1], alternative="greater", method="pearson", conf.level = 0.95)

# between time and numerosity
cor.test(dsa$width[dsa$task==0], dsa$width[dsa$task==2], alternative="greater", method="pearson", conf.level = 0.95)

# between space and numerosity
cor.test(dsa$width[dsa$task==1], dsa$width[dsa$task==2], alternative="greater", method="pearson", conf.level = 0.95)


### test correlation coefficients against each other (https://www.personality-project.org/r/html/paired.r.html)
library('psych')

# time-space-correlation vs space-numerosity-correlation:
paired.r(0.389308, -0.02580663, 0.2430187, n=24, twotailed=F)

# time-space-correlation vs time-numerosity-correlation:
paired.r(0.389308, 0.2430187, -0.02580663, n=24, twotailed=F)







###################################################################################
### cross-dimensional interferences

load("pse_data_cross_dim_individual_norm.RData")
dsa <- pse_data_cross_dim_individual_norm


### delete outliers:
# for space-on-time interference
depDim = 0
indDim = 1
m = mean(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
s = sd(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
nrow(dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),])

# for numerosity-on-time interference
depDim = 0
indDim = 2
m = mean(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
s = sd(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
nrow(dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),])

# for time-on-space interference
depDim = 1
indDim = 0
m = mean(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
s = sd(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
nrow(dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),])
dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),]
# exclude this outlier: interference effect from time on space for subject 04
dsa$thr50_diff[dsa$subj==4 & dsa$task==1 & dsa$condTask==0] <- NA

# for numerosity-on-space interference
depDim = 1
indDim = 2
m = mean(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
s = sd(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
nrow(dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),])

# for time-on-numerosity interference
depDim = 2
indDim = 0
m = mean(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
s = sd(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
nrow(dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),])

# for space-on-numerosity interference
depDim = 2
indDim = 1
m = mean(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
s = sd(dsa$thr50_diff[dsa$task==depDim & dsa$condTask==indDim])
nrow(dsa[dsa$task==depDim & dsa$condTask==indDim & (dsa$thr50_diff>m+(3*s) | dsa$thr50_diff<m-(3*s)),])



### significance of interference effects:

# effect of space on time:
t.test(dsa$thr50_diff[dsa$condTask==1 & dsa$task==0], NULL, alternative="two.sided", paired=F, var.equal=FALSE, con.level=0.95)  
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(dsa$thr50_diff[dsa$condTask==1 & dsa$task==0], na.rm=T)
sd1 = sd(dsa$thr50_large[dsa$condTask==1 & dsa$task==0], na.rm=T)
sd2 = sd(dsa$thr50_small[dsa$condTask==1 & dsa$task==0], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)

# effect of numerosity on time:
t.test(dsa$thr50_diff[dsa$condTask==2 & dsa$task==0], NULL, alternative="two.sided", paired=F, var.equal=FALSE, con.level=0.95)  
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(dsa$thr50_diff[dsa$condTask==2 & dsa$task==0], na.rm=T)
sd1 = sd(dsa$thr50_large[dsa$condTask==2 & dsa$task==0], na.rm=T)
sd2 = sd(dsa$thr50_small[dsa$condTask==2 & dsa$task==0], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)

# effect of time on space:
t.test(dsa$thr50_diff[dsa$condTask==0 & dsa$task==1], NULL, alternative="two.sided", paired=F, var.equal=FALSE, con.level=0.95)  
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(dsa$thr50_diff[dsa$condTask==0 & dsa$task==1], na.rm=T)
sd1 = sd(dsa$thr50_large[dsa$condTask==0 & dsa$task==1], na.rm=T)
sd2 = sd(dsa$thr50_small[dsa$condTask==0 & dsa$task==1], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)

# effect of numerosity on space:
t.test(dsa$thr50_diff[dsa$condTask==2 & dsa$task==1], NULL, alternative="two.sided", paired=F, var.equal=FALSE, con.level=0.95)  
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(dsa$thr50_diff[dsa$condTask==2 & dsa$task==1], na.rm=T)
sd1 = sd(dsa$thr50_large[dsa$condTask==2 & dsa$task==1], na.rm=T)
sd2 = sd(dsa$thr50_small[dsa$condTask==2 & dsa$task==1], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)

# effect of time on numerosity:
t.test(dsa$thr50_diff[dsa$condTask==0 & dsa$task==2], NULL, alternative="two.sided", paired=F, var.equal=FALSE, con.level=0.95)  
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(dsa$thr50_diff[dsa$condTask==0 & dsa$task==2], na.rm=T)
sd1 = sd(dsa$thr50_large[dsa$condTask==0 & dsa$task==2], na.rm=T)
sd2 = sd(dsa$thr50_small[dsa$condTask==0 & dsa$task==2], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)

# effect of space on numerosity:
t.test(dsa$thr50_diff[dsa$condTask==1 & dsa$task==2], NULL, alternative="two.sided", paired=F, var.equal=FALSE, con.level=0.95)  
# calculate Cohen's d(av) (according to Lakens (2013) Front Psychol: Equation 10 on page 5)
mdiff = mean(dsa$thr50_diff[dsa$condTask==1 & dsa$task==2], na.rm=T)
sd1 = sd(dsa$thr50_large[dsa$condTask==1 & dsa$task==2], na.rm=T)
sd2 = sd(dsa$thr50_small[dsa$condTask==1 & dsa$task==2], na.rm=T)
abs(mdiff)/((sd1+sd2)/2)



### direct comparison between interference effects

# effect "space on time" is larger than effect "numerosity on time":
t.test(dsa$thr50_diff[dsa$condTask==1 & dsa$task==0], dsa$thr50_diff[dsa$condTask==2 & dsa$task==0], alternative="two.sided", paired=T, var.equal=FALSE, con.level=0.95)  
# effect "time on space" is larger than effect "numerosity on space":
t.test(dsa$thr50_diff[dsa$condTask==0 & dsa$task==1], dsa$thr50_diff[dsa$condTask==2 & dsa$task==1], alternative="two.sided", paired=T, var.equal=FALSE, con.level=0.95)  
# effect "time on numerosity" is not different from effect "space on numerosity":
t.test(dsa$thr50_diff[dsa$condTask==0 & dsa$task==2], dsa$thr50_diff[dsa$condTask==1 & dsa$task==2], alternative="two.sided", paired=T, var.equal=FALSE, con.level=0.95)  
# effect "space on time" is larger than effect "time on space":
t.test(dsa$thr50_diff[dsa$condTask==1 & dsa$task==0], dsa$thr50_diff[dsa$condTask==0 & dsa$task==1], alternative="two.sided", paired=T, var.equal=FALSE, con.level=0.95)  




