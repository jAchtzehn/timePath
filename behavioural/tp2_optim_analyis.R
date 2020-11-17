load("/Users/jachtzehn/Documents/Development/dzne/branches/timePath/behavioral/r/tp2.behav.stana.basepointdata.RData")
dsa <- tp2.behav.stana.basepointdata



### calculate relative comparison value (relative to a standard of 1)
dsa$relComp <- NA
dsa$relComp <- ifelse(dsa$task==0 & dsa$stand==1, dsa$comp/2.8, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==0 & dsa$stand==2, dsa$comp/4.8, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==1 & dsa$stand==1, dsa$comp/11.5, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==1 & dsa$stand==2, dsa$comp/19.7, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==2 & dsa$stand==1, dsa$comp/45, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==2 & dsa$stand==2, dsa$comp/77, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==3 & dsa$stand==1, dsa$comp/0.16, dsa$relComp)
dsa$relComp <- ifelse(dsa$task==3 & dsa$stand==2, dsa$comp/0.28, dsa$relComp)



### draw logistic function using optim
drawLogFun <- function(v,t,s,ct,cs,color,linetype,pointtype,slStartVal,thrStartVal,minx,maxx){
  sub <- subset(dsa,VP==v & task==t & stand==s & condTask==ct & condStand==cs)
  x <- sub$comp
  r <- sub$numR1
  m <- sub$abs_sum_resp
  fn <- function(par){
    rhat <- 1/(1+exp(1)^(-par[1]*(x-par[2])))
    sum((r/m - rhat)^2)
  }
  a<-optim(c(slStartVal,thrStartVal), fn, method="BFGS", control=list(reltol=1e-9))
  points(r/m~x, col=color, pch=pointtype)
  z <- seq(minx,maxx,0.01)
  y <- 1/(1+exp(1)^(-a$par[1]*(z-a$par[2])))
  lines(y~z, col=color, lty=linetype)
  return(plot)
}



### get parameters for logistic function using optim
getLogFunParam <- function(v,t,s,ct,cs,slStartVal,thrStartVal){
  sub <- subset(dsa,VP==v & task==t & stand==s & condTask==ct & condStand==cs)
  x <- sub$comp
  r <- sub$numR1
  m <- sub$abs_sum_resp
  fn <- function(par){
    rhat <- 1/(1+exp(1)^(-par[1]*(x-par[2])))
    sum((r/m - rhat)^2)
  }
  a<-optim(c(slStartVal,thrStartVal), fn, method="BFGS", control=list(reltol=1e-9))
  thr25 <- ((log((1/0.25)-1))/-a$par[1])+a$par[2]
  thr50 <- a$par[2]
  thr75 <- ((log((1/0.75)-1))/-a$par[1])+a$par[2]
  sl50 <- a$par[1]
  dat <- data.frame(v,t,s,ct,cs,thr25,thr50,thr75,sl50)
  return(dat)
}



### draw functions for a specific subject (remember: sub-01 == vp08)
i <- 15
par(mfrow=c(4,2))
for (t in 0:3){
  thrStartVal <- mean(dsa$comp[dsa$task==t & dsa$stand==1 & dsa$VP==i])
  minx <- min(dsa$comp[dsa$task==t & dsa$stand==1])# & dsa$VP==i])
  maxx <- max(dsa$comp[dsa$task==t & dsa$stand==1])# & dsa$VP==i])
  plot(-4,1,col="white",xlim=c(minx,maxx),ylim=c(0,1), main=t, ylab=i)
  # drawLogFun(VP,task,stand,condTask,condStand,color,linetype,pointtype,slStartVal,thrStartVal,minx,maxx)
  drawLogFun(i,t,1,t,1,"black",1,1,1,thrStartVal,minx,maxx)
  
  thrStartVal <- mean(dsa$comp[dsa$task==t & dsa$stand==2 & dsa$VP==i])
  minx <- min(dsa$comp[dsa$task==t & dsa$stand==2])# & dsa$VP==i])
  maxx <- max(dsa$comp[dsa$task==t & dsa$stand==2])# & dsa$VP==i])
  plot(-4,1,col="white",xlim=c(minx,maxx),ylim=c(0,1), main=t, ylab=i)
  # drawLogFun(VP,task,stand,condTask,condStand,color,linetype,pointtype,slStartVal,thrStartVal,minx,maxx)
  drawLogFun(i,t,2,t,2,"black",1,1,1,thrStartVal,minx,maxx)
}



# get new dataset with parameters
newdsa <- c(1:9)
for (t in 0:3){
  for (i in unique(dsa$VP)){
    for (j in 1:2){
      thrStartVal <- mean(dsa$comp[dsa$task==t & dsa$stand==j & dsa$VP==i])
      # getLogFunParam(VP,task,stand,condTask,condStand,slStartVal,thrStartVal)
      b <- getLogFunParam(i,t,j,t,j,1,thrStartVal)
      newdsa <- rbind(newdsa,b)
    }
  }
}
dsa2 <- newdsa
dsa2 <- dsa2[2:(nrow(dsa2)),]   # delete dummy row   

names(dsa2)[names(dsa2)=="v"] <- "VP"
names(dsa2)[names(dsa2)=="t"] <- "task"
names(dsa2)[names(dsa2)=="s"] <- "stand"
names(dsa2)[names(dsa2)=="ct"] <- "condTask"
names(dsa2)[names(dsa2)=="cs"] <- "condStand"

dsa2$dl <- (dsa2$thr75 - dsa2$thr25)/2
dsa2$dl.log <- log(dsa2$dl)





