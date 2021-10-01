reentry_data<-read.spss("transformed data.sav", to.data.frame = T)
saveRDS(reentry_data,file="reentry_data.RDS")
reentry_data<-readRDS("reentry_data.RDS")

require(lme4)
require(lcmm)
require(dplyr)
require(ggplot2)
require(tidyr)

#modeling time
linear_model<-lmer(dynamic~time_assess+(1+time_assess|ID),
                    data=reentry_data)
summary(linear_model)
quad_model<-lmer(dynamic~poly(time_assess,2)+(1+time_assess|ID),
                  data=reentry_data)
summary(quad_model)
cubic_model<-lmer(dynamic~poly(time_assess,3)+(1+time_assess|ID),	data=reentry_data)
summary(cubic_model)
anova(linear_model,quad_model,cubic_model)

#adding an optimizer
Quad_model<-lmer(dynamic~poly(time_assess,2)+(1+time_assess|ID),
                 data=reentry_data,
                 control = lmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1e5)))
summary(Quad_model)

#test slope and intercept variability
m1<-lmer(dynamic~(1|ID),data= reentry_data)
summary(m1)
m2<-lmer(dynamic~time_assess+(1|ID),data= reentry_data)
summary(m2)
m3<-lmer(dynamic~time_assess+(1+time_assess|ID),data= reentry_data)
summary(m3)
anova(m1,m2,m3)

#initial model
mod_init<-Jointlcmm(dynamic~time_assess+static,random=~time_assess,
                      survival=Surv(time,status)~static,
                      hazard="Weibull",subject="ID",data=reentry_data,
                    ng=1,link="splines")
summary(mod_init)

#selection models

#2 groups
mod_2g<-gridsearch(rep=3,maxiter=50,minit=mod_init,
                 Jointlcmm(dynamic~time_assess+static,
                           mixture=~time_assess,
                           random=~time_assess,
                           survival=Surv(time,status)~static,
                           classmb=~static,
                           hazard="Weibull",subject="ID",
                           data=reentry_data,
                           ng=2,link="splines",verbose=FALSE))
summary(mod_2g)
postprob(mod_2g)

#2 groups
mod_2g<-gridsearch(rep=3,maxiter=50,minit=mod_init,
                   Jointlcmm(dynamic~time_assess+static,
                             mixture=~time_assess,
                             random=~time_assess,
                             survival=Surv(time,status)~static,
                             classmb=~static,
                             hazard="Weibull",subject="ID",
                             data=reentry_data,
                             ng=2,link="splines",verbose=FALSE))
summary(mod_2g)
postprob(mod_2g)

entropy <- function(mod_2g) {
  pp <- mod_2g$pprob[3:ncol(mod_2g$pprob)]
  N <- nrow(pp)                   # number of individuals
  J <- mod_2g$ng                     # number of latent classes
  ent <- -sum(pp * log(pp))       # entropy
  res <- 1 - (ent / (N * log(J))) # relative entropy
  res
}
entropy(mod_2g)

#plotting
newdata<-data.frame(time_assess=seq(0,25,length=50),
                    static=rep(0.46,50))
mod2_predtraj<-predictY(mod_2g,newdata,	var.time="time_assess")
plot(mod2_predtraj)

#3 groups
mod_3g<-gridsearch(rep=5,maxiter=100,minit=mod_init,
                   Jointlcmm(dynamic~time_assess+static,
                             mixture=~time_assess,
                             random=~time_assess,
                             survival=Surv(time,status)~static,
                             classmb=~static,
                             hazard="Weibull",subject="ID",
                             data=reentry_data,
                             ng=3,link="splines",verbose=FALSE))
summary(mod_3g)
postprob(mod_3g)

entropy <- function(mod_3g) {
  pp <- mod_3g$pprob[3:ncol(mod_3g$pprob)]
  N <- nrow(pp)                   # number of individuals
  J <- mod_3g$ng                     # number of latent classes
  ent <- -sum(pp * log(pp))       # entropy
  res <- 1 - (ent / (N * log(J))) # relative entropy
  res
}
entropy(mod_3g)

#4 groups
mod_4g<-gridsearch(rep=6,maxiter=100,minit=mod_init,
                   Jointlcmm(dynamic~time_assess+static,
                             mixture=~time_assess,
                             random=~time_assess,
                             survival=Surv(time,status)~static,
                             classmb=~static,
                             hazard="Weibull",subject="ID",
                             data=reentry_data,
                             ng=4,link="splines",verbose=FALSE))
summary(mod_4g)
postprob(mod_4g)

entropy <- function(mod_4g) {
  pp <- mod_4g$pprob[3:ncol(mod_4g$pprob)]
  N <- nrow(pp)                   # number of individuals
  J <- mod_4g$ng                     # number of latent classes
  ent <- -sum(pp * log(pp))       # entropy
  res <- 1 - (ent / (N * log(J))) # relative entropy
  res
}
entropy(mod_4g)

#5 groups
mod_5g<-gridsearch(rep=5,maxiter=100,minit=mod_init,
                   Jointlcmm(dynamic~time_assess+static,
                             mixture=~time_assess,
                             random=~time_assess,
                             survival=Surv(time,status)~static,
                             classmb=~static,
                             hazard="Weibull",subject="ID",
                             data=reentry_data,
                             ng=5,link="splines",verbose=FALSE))
summary(mod_5g)
postprob(mod_5g)

entropy <- function(mod_5g) {
  pp <- mod_5g$pprob[3:ncol(mod_5g$pprob)]
  N <- nrow(pp)                   # number of individuals
  J <- mod_5g$ng                     # number of latent classes
  ent <- -sum(pp * log(pp))       # entropy
  res <- 1 - (ent / (N * log(J))) # relative entropy
  res
}
entropy(mod_5g)

#Plot predicted mean trajectories, hazard rates, and survival curves
newdata<-data.frame(time_assess=seq(0,25,length=50),
                    static=rep(0.46,50))
predtraj<-predictY(mod_2g,newdata,	var.time="time_assess")
plot(predtraj)
#with confidence intervals
predtraj_MC<-predictY(mod_2g,newdata,draws=TRUE,var.time="time_assess")
plot(predtraj, bty = "l",
     ylim = c(0, 10), xlim=c(0,25),
     ylab = "Dynamic Risk",
     xlab = "Week", lwd = 1,col=1)
#hazard plot
plot(mod_2g, which = "hazard",
     lwd = 2,bty = "l",
     xlab = "Weeks",
     ylab = "Hazard for any Recidivism",
     xlim=c(0,25),col=1)
#survival plot
plot(mod_2g, which = "survival",
     lwd = 2,bty = "l",
     xlab = "Weeks",
     ylab = "Recidivism-Free probability",
     xlim=c(0,25),col=1)

#illustrate within group noise
m2_pprob<-as.data.frame(mod_2g$pprob[,1:4])
reentry_data<-left_join(reentry_data, m2_pprob)
#extract group 1 sample trajectories
m2_group1<-reentry_data[which(reentry_data$class==1),]
ID<-as.character(m2_group1$ID)
ID<-as.data.frame(ID,ID=ID)
ID<-ID[!duplicated(ID$ID),]
ID<-as.character(ID)
ID<-as.data.frame(ID, ID=ID)
ID$ID<-as.character(ID$ID)
m2_group1_sample<-sample(x=ID$ID, size = 50, replace = F, 
                                set.seed(1))
m2_group1_sample<-data.frame(m2_group1_sample, ID=m2_group1_sample)
m2_group1_sample<-m2_group1[m2_group1$ID%in%m2_group1_sample$ID,]
#extract mean trajectory
predtraj_df<-as.data.frame(predtraj$pred[,1:2],predtraj$times[,1])
#plot
p1<-ggplot()+geom_line(data=m2_group1_sample,
                       aes(x=time_assess,y=dynamic,group=ID),position="jitter")+
  geom_line(data=predtraj_df,aes(x=as.numeric(row.names(predtraj_df)),y=Ypred_class1,size=12))+
  ylab("Dynamic Score")+
  xlab("Weeks in Community")+
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+
  scale_x_continuous(breaks=seq(0,25,2))+
  scale_y_continuous(breaks=seq(0,15,1))+
  theme_classic()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
p1

#prepare for AUC and Brier Score testing 
dpred <-dynpred(mod_2g, landmark = c(4,12,24),
                newdata = reentry_data  ,horizon=c(10),
                var.time = "time_assess", event = 1)

pred_df<-dpred$pred%>%dplyr::select(ID, landmark, pred)
pred_df<-pivot_wider(pred_df, id_cols = ID,
                     names_from = landmark, values_from = pred)

timestatus<-reentry_data[c("ID","time","status")]
timestatus<-timestatus[!duplicated(timestatus),]
pred_df<-left_join(pred_df,timestatus )

names(pred_df)<-c("ID","pred.dynamic.4","pred.dynamic.12","pred.dynamic.24","time","status")
write.csv(pred_df,"pred_df.csv")