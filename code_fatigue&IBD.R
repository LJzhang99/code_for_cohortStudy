#~加载R包####
library(tidyr)
library(dplyr)
library(plyr)
library(tibble)
library(lubridate)
library(readr)
library(data.table)
library(stringr)
library(Hmisc)
library(broom)

setwd("~/fatigue_IBD")
#fatigue#####
fatigue<-read_csv("./data/2080_participant.csv")
fatigue=fatigue[,1:2]
#二分类 ###
#一半以上为疲劳

fatigueY<-filter(fatigue,`2080-0.0`== "More than half the days" | `2080-0.0`=="Nearly every day")
fatigueY$fatigue<-as.factor(1:1)
fatigueN<-filter(fatigue,`2080-0.0`== "Several days" | `2080-0.0`=="Not at all")
fatigueN$fatigue<-as.factor(0:0)
fatigue<-full_join(fatigueY,fatigueN)
fatigue <- fatigue %>%
  mutate_at(vars(2), ~case_when(
    . == "Not at all" ~ 0,
    . == "Several days" ~ 1,
    . == "More than half the days" ~ 2,
    . == "Nearly every day" ~ 3,
    TRUE ~ as.numeric(.)
  ))
#起始时间###
baseline<-read_csv("/home/data/PubData/epidemiological/ukb_original/53.csv")
baseline=baseline[,c(1,2)]
names(baseline)<-c("eid","beginDate")
fatigue<-inner_join(fatigue,baseline)
colnames(fatigue)[2]<-c("tiredness")

#IBD
#ICD 10
digestive<-read_csv("./data/digestiveDisease_firstOccur.csv") #first occur 
CD<-na.omit(digestive[,c(1,24)])

UC<-na.omit(digestive[,c(1,25)])

IBD<-full_join(CD,UC)
colnames(IBD)<-c("eid","CD","UC")

IBD<-data.frame(apply(IBD,2,function(x){
  x[is.na(x)] = "2222-02-09";x}))
IBD$CD<-as.Date(IBD$CD)
IBD$UC<-as.Date(IBD$UC)

getMinTime <- function (t) {
  return(min(t))
}

minTime<-apply(IBD[,c(2,3)],1,getMinTime)
IBD[c("endTime")] = as.Date(minTime) #获得诊断IBD的最早时间
IBD$status<-as.numeric(1:1)
IBD$eid<-as.double(IBD$eid)
eid<-read_csv("~/dataCommon/eid.csv")
IBDno<-anti_join(eid,IBD)
IBDno$status<-as.numeric(0:0)
IBDno$statusUC<-as.numeric(0:0)
IBDno$statusCD<-as.numeric(0:0)
IBDno$endTime<-as.Date("2022-10-30")
IBDall<-full_join(IBD,IBDno)
IBDall<-IBDall[,c(1,4,5)]

pheo<-IBDall[,c(1,5)]
 pheo$IID<-pheo$eid
 pheo<-pheo[,c(1,3,2)]
 names(pheo)<-c("FID","IID","IBD")
 fwrite(pheo,"Pheno.IBD",sep="\t")
TableFatigueIBD<-inner_join(fatigue,IBDall)
TableFatigueIBD<-filter(TableFatigueIBD,beginDate < endTime) #排除随访前诊断IBD

#死亡和失访#
#死亡
death<-read_csv("/mnt/UKBdataNew/40000.csv")
death=death[,1:2]
death<-na.omit(death)
death<-unique(death)
names(death)<-c("eid","deathDate")
TableFatigueIBD<-left_join(TableFatigueIBD,death)
TableFatigueIBD$deathDate = impute(TableFatigueIBD$deathDate, "3000-01-01") 
#死亡提前结束随访
TableFatigueIBD$endTime <- ifelse(TableFatigueIBD$deathDate <= TableFatigueIBD$endTime , TableFatigueIBD$deathDate, TableFatigueIBD$endTime)
TableFatigueIBD$endTime <-as.Date(TableFatigueIBD$endTime,format = "%Y-%m-%d",origin = "1970-01-01")

#排除失访
exclude<-read_csv("~/longCovid/data/X190.csv")
exclude<-na.omit(exclude)
e<-inner_join()
TableFatigueIBD<-anti_join(TableFatigueIBD,exclude)

#survival time
TableFatigueIBD[c("time")]<- as.numeric(difftime(TableFatigueIBD$endTime,TableFatigueIBD$beginDate,
                                         units = c("days"))/365.25)
#IBD来源#####
#排除self
source<-read_csv("data/source_of_IBD_participant.csv")
names(source)<-c("eid","CD","UC")
Selfreport<-filter(source,CD == 50 | UC == 50)

#协变量####
#失眠#
sleep<-read_csv("~/dataCommon/sleep_participant.csv")
sleep<-sleep[,c(1,6)]
sleep <- sleep %>%
  mutate_at(vars(2), ~case_when(
    . == "Never/rarely" ~ 0,
    . == "Sometimes" ~ 1,
    . == "Usually" ~ 2,
    TRUE ~ as.numeric(.)
  ))
colnames(sleep)[2]<-c("sleepness")

#education college==1
education<-read_csv("/mnt/UKBdata/6138.csv")
education=education[,c(1,2)]
education=na.omit(education)

education1<-education[apply(education[,c("6138-0.0")],1,function(x){any(grepl("1",x))}),]
education1$education<-as.factor(1:1)
education0<-anti_join(education,education1)
education0$education<-as.factor(0:0)
education<-full_join(education0,education1)
education=education[,c(1,3)]

#其他协变量##
covariate<-read_csv("/mnt3/mayuying/fatigue_IBD/cova_Before_impute.csv")

covariate<-inner_join(covariate,sleep)
covariate<-left_join(covariate,education)
varsToFactor <- c("medication","vitamin","Ethnicity","house","alcohol","smoking","sex"
                ,"sleepness")
covariate[varsToFactor]<-lapply(covariate[varsToFactor],factor)


#插补
library(mice)

cova_mice<-mice(covariate,m=5,maxit=50,meth="pmm",seed=500)

summary(cova_mice)
resultCova<-complete(cova_mice,action=3) 
sum(is.na(resultCova))

write_csv(resultCova,file="~/dataCommon/covariateImpute.csv")
covariate<-read_csv("~/dataCommon/covariateImpute.csv")
#补充协变量
addCova<-read_csv("/mnt3/mayuying/Previous/MS_GERD/covBaseline.csv")
addCova<-addCova[,c(1,15:23)]
varsFactor <- c("hypertension","heartFailure","renalFailure","asthma","dementia","MI","stroke",
                  "copd","diabetes")
addCova[varsFactor]<-lapply(addCova[varsFactor],factor)
#焦虑抑郁
Mental<-read_csv("~/fatigue_IBD/data/depression_anxiety_participant.csv")
Depression<-Mental[,c(1:3)]
Anxiety<-Mental[,c(1,4:5)]

Depression<-inner_join(Depression,baseline)
Anxiety<-inner_join(Anxiety,baseline)

X130894 <- Mental[,c(1,2)]
X130896 <- Mental[,c(1,3)]
X130904 <- Mental[,c(1,4)]
X130906 <- Mental[,c(1,5)]
write_csv(X130894,"130894.csv")
write_csv(X130896,"130896.csv")
write_csv(X130904,"130904.csv")
write_csv(X130906,"130906.csv")

names(Depression)<-c("eid","Episode","Recurrent","baseline")
names(Anxiety)<-c("eid","Phobic","other","baseline")

##新药物####
medication<-read_csv("~/fatigue_IBD/data/medication_NA.csv")

eid<-read_csv("~/dataCommon/eid.csv")
Depression<-filter(Depression,Episode <= baseline | Recurrent <= baseline)
Depression$depression<-as.factor(1:1)
DepressionN<-anti_join(eid,Depression)
DepressionN$depression<-as.factor(0:0)
DepressionAll<-full_join(Depression,DepressionN)
DepressionAll=DepressionAll[,c(1,5)]

Anxiety<-filter(Anxiety, Phobic <= baseline | other <= baseline)
Anxiety$anxiety<-as.factor(1:1)
AnxietyN<-anti_join(eid,Anxiety)
AnxietyN$anxiety<-as.factor(0:0)
AnxietyAll<-full_join(Anxiety,AnxietyN)
AnxietyAll=AnxietyAll[,c(1,5)]

Mental<-full_join(AnxietyAll,DepressionAll)
Mental<-Mental[,c(1,2,6)]
write_csv(Mental,"data/MentaBaseline.csv")
Mental<-read_csv("data/MentaBaseline.csv")
  
#合并协变量
TableFatigueIBD<-inner_join(TableFatigueIBD,resultCova)
TableFatigueIBD<-inner_join(TableFatigueIBD,Mental)
TableFatigueIBD<-inner_join(TableFatigueIBD,addCova)
TableFatigueIBD<-TableFatigueIBD[,-c(9)]
TableFatigueIBD<-inner_join(TableFatigueIBD,PRS)
str(TableFatigueIBD)
summary(TableFatigueIBD)

#排除self-reported
TableFatigueIBD<-anti_join(TableFatigueIBD,Selfreport) #排除27人
TableFatigueIBD<-inner_join(TableFatigueIBD,medication)
#IBD分析######
varsToFactor <- c("hypertension","heartFailure","renalFailure","asthma","dementia","MI","stroke",
                "copd","diabetes","medication","vitamin","Ethnicity","house","alcohol","smoking","sex",
                "education","sleepness","anxiety","depression")
TableFatigueIBD[varsToFactor]<-lapply(TableFatigueIBD[varsToFactor],factor)
TableFatigueIBD$fatigue<-factor(TableFatigueIBD$fatigue,levels = c("0","1"))

coxFactor<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                   vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+dementia+MI+
                   stroke+copd+diabetes+prsCD+prsUC,
                data=TableFatigueIBD)
summary(coxFactor)


coxFactor<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI,
                 data=TableFatigueIBD)

reduced_model <- coxph(Surv(time, status) ~ 1, data = TableFatigueIBD)

# Calculate predicted survival probabilities
predicted_survival <- predict(coxFactor, type = "expected")

# Calculate baseline survival probabilities
baseline_survival <- survfit(Surv(time, status) ~ 1, data = TableFatigueIBD)

# Extract survival probabilities at each event time
baseline_survival_probs <- baseline_survival$surv

# Calculate Cox-Snell residuals
cox_snell_residuals <- qnorm(predicted_survival / baseline_survival_probs)

# Calculate the proportion of variance explained
variance_elucidation <- 1 - exp(-sum(cox_snell_residuals^2) / (2 * length(TableFatigueIBD$status)))
variance_elucidation 

testF

ggsurvplot(survfit(Surv(time,status)~ fatigue, data = TableFatigueIBD),
           pval = TRUE, conf.int = TRUE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B"),
           ylim=c(0.98,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

TableFatigueIBD$tiredness<-factor(TableFatigueIBD$tiredness,levels = c("1","2","3","4"))
coxNum<-coxph(Surv(time,status)~ tiredness + sex + Ethnicity +smoking + alcohol + age + BMI +
                   vegetable + medication+ + vitamin + house + fruit + MET + TDindex  +sleepness,
                 data=TableFatigueIBD)
summary(coxNum)
testN<-cox.zph(coxNum)
testN

ggsurvplot(survfit(Surv(time,status)~ tiredness, data = TableFatigueIBD),
           ylim=c(0.98,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

#剂量反应曲线#
install.packages("drc")
library(drc)


##UC#####
#排除UC前CD# 
UC<-filter(IBD,UC<=CD)
UC$statusUC<-as.numeric(1:1)
UC=UC[,c(1,3,5)]
names(UC)<-c("eid","endTime","statusUC")
IBDno$statusUC<-as.numeric(0:0)
TableUC<-full_join(UC,IBDno)
TableUC<-TableUC[,-4]
TableUC<-inner_join(TableUC,baseline)
TableUC<-filter(TableUC,beginDate < endTime)
TableUC<-left_join(TableUC,death)
TableUC$deathDate = impute(TableUC$deathDate, "3000-01-01") 
#死亡提前结束随访
TableUC$endTime <- ifelse(TableUC$deathDate <= TableUC$endTime , TableUC$deathDate, TableUC$endTime)
TableUC$endTime <-as.Date(TableUC$endTime,format = "%Y-%m-%d",origin = "1970-01-01")
TableUC<-anti_join(TableUC,exclude)
TableUC<-inner_join(TableUC,resultCova)
TableUC[c("time")]<- as.numeric(difftime(TableUC$endTime,TableUC$beginDate,
                                                 units = c("days"))/365.25)
str(TableUC)
summary(TableUC)
TableUC<-inner_join(TableUC,fatigue)
#UC分析######
TableUC<-inner_join(TableUC,addCova)
TableUC<-inner_join(TableUC,Mental)
TableUC<-inner_join(TableUC,PRS)
TableUC<-anti_join(TableUC,Selfreport)
TableUC<-inner_join(TableUC,medication)

write_csv(TableUC,"TableFatigueUC.csv")
TableUC<-read_csv("TableUCnew.csv")
TableUC[varsToFactor]<-lapply(TableUC[varsToFactor],factor)
TableUC$fatigue<-factor(TableUC$fatigue,levels = c("0","1"))

TableUC$PRSGroupUC<-as.factor(TableUC$PRSGroupUC)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                   vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+
                   copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)


coxFactor<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI ,
                 data=TableUC)
reduced_model <- coxph(Surv(time, statusUC) ~ 1, data = TableUC)

# Calculate predicted survival probabilities
predicted_survival <- predict(coxFactor, type = "expected")

# Calculate baseline survival probabilities
baseline_survival <- survfit(Surv(time, statusUC) ~ 1, data = TableUC)

# Extract survival probabilities at each event time
baseline_survival_probs <- baseline_survival$surv

# Calculate Cox-Snell residuals
cox_snell_residuals <- qnorm(predicted_survival / baseline_survival_probs)

# Calculate the proportion of variance explained
variance_elucidation <- 1 - exp(-sum(cox_snell_residuals^2) / (2 * length(TableUC$statusUC)))
variance_elucidation 

testF<-cox.zph(coxFactor)
testF
TableUC$tiredness<-factor(TableUC$tiredness,levels = c("0","1","2","3"))

coxNum<-coxph(Surv(time,statusUC)~ tiredness + sex + Ethnicity +smoking + alcohol + age + BMI +
                vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
              data=TableUC)
summary(coxNum)
testN<-cox.zph(coxNum)
testN
##CD#####
#排除CD前CD# 
CD<-filter(IBD,CD<=UC)
CD$statusCD<-as.numeric(1:1)
CD=CD[,c(1,2,6)]
names(CD)<-c("eid","endTime","statusCD")
IBDno$statusCD<-as.numeric(0:0)
TableCD<-full_join(CD,IBDno)
TableCD<-TableCD[,-c(4:5)]
TableCD<-inner_join(TableCD,baseline)
TableCD<-filter(TableCD,beginDate < endTime)
TableCD<-left_join(TableCD,death)
TableCD$deathDate = impute(TableCD$deathDate, "3000-01-01") 
#死亡提前结束随访
TableCD$endTime <- ifelse(TableCD$deathDate <= TableCD$endTime , TableCD$deathDate, TableCD$endTime)
TableCD$endTime <-as.Date(TableCD$endTime,format = "%Y-%m-%d",origin = "1970-01-01")
TableCD<-anti_join(TableCD,exclude)
TableCD<-inner_join(TableCD,resultCova)
TableCD[c("time")]<- as.numeric(difftime(TableCD$endTime,TableCD$beginDate,
                                         units = c("days"))/365.25)
str(TableCD)
summary(TableCD)
TableCD<-inner_join(TableCD,fatigue)
#UC分析######
TableCD<-inner_join(TableCD,addCova)
TableCD<-inner_join(TableCD,Mental)

TableCD=TableCD[,-c(6)]
TableCD<-inner_join(TableCD,medication)
write_csv(TableCD,"TableFatigueCD.csv")

TableCD<-na.omit(read_csv("TableCDnew.csv"))
TableCD[varsToFactor]<-lapply(TableCD[varsToFactor],factor)
TableCD$fatigue<-factor(TableCD$fatigue,levels = c("0","1"))
coxFactor<-coxph(Surv(time,statusCD)~ fatigue ,
                 data=TableCD)
summary(coxFactor)

reduced_model <- coxph(Surv(time, statusCD) ~ 1, data = TableCD)

# Calculate predicted survival probabilities
predicted_survival <- predict(coxFactor, type = "expected")

# Calculate baseline survival probabilities
baseline_survival <- survfit(Surv(time, statusCD) ~ 1, data = TableCD)

# Extract survival probabilities at each event time
baseline_survival_probs <- baseline_survival$surv

# Calculate Cox-Snell residuals
cox_snell_residuals <- qnorm(predicted_survival / baseline_survival_probs)

# Calculate the proportion of variance explained
variance_elucidation <- 1 - exp(-sum(cox_snell_residuals^2) / (2 * length(TableCD$statusCD)))
variance_elucidation 

testF<-cox.zph(coxFactor)
testF
coxNum<-coxph(Surv(time,statusCD)~ tiredness + sex + Ethnicity +smoking + alcohol + age + BMI +
                vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
              data=TableCD)
summary(coxNum)
testN<-cox.zph(coxNum)
testN


#PRS####
PRS<-read_csv("data/UC_CD_PRS_participant.csv")
names(PRS)<-c("eid","prsCD","prsUC")
PRS<-inner_join(PRS,eid)
PRS<-inner_join(covariate,PRS)
PRS<-inner_join(PRS,addCova)
PRS<-inner_join(PRS,Mental)

#prs插补####
PRS[varsToFactor]<-lapply(PRS[varsToFactor],factor)
PRS=PRS[,c(1,8:13,17:18)]
Prs_mice<-mice(PRS,m=5,maxit=50,meth="pmm",seed=500)

summary(Prs_mice)
resultPRS<-complete(Prs_mice,action=3) 
sum(is.na(resultPRS))
resultPRS=resultPRS[,c(1,8:9)]
#不插补去掉没有PRS
PRS=na.omit(PRS)
UCprs<-na.omit(PRS[,c(1,3)])
CDprs<-na.omit(PRS[,c(1,2)])

TableUC<-inner_join(TableUC,UCprs)
TableUC$PRSGroup<-cut(TableUC$prsUC,breaks = c(-5,-0.3004540,0.5481677,10) ,labels=FALSE)
TableUC$PRSGroupUC<-cut(TableUC$prsUC,quantile(TableUC$prsUC,seq(0,1,1/3)) ,labels=FALSE)
TableUC$PRSGroupCD<-cut(TableUC$prsUC,quantile(TableUC$prsCD,seq(0,1,1/3)) ,labels=FALSE)

#用插补后的PRS####
resultPRS<-read_csv("data/PRSmultiple.csv")
TableFatigueIBD=TableFatigueIBD[,-c(35:36)]
TableFatigueIBD<-inner_join(TableFatigueIBD,resultPRS)
TableFatigueIBD$PRSGroupCD<-cut(TableFatigueIBD$prsCD,quantile(TableFatigueIBD$prsCD,seq(0,1,1/3)) ,labels=FALSE)
TableFatigueIBD$PRSGroupUC<-cut(TableFatigueIBD$prsUC,quantile(TableFatigueIBD$prsUC,seq(0,1,1/3)) ,labels=FALSE)


#IBD+prs######
TableFatigueIBD<-read_csv("TableIBDnew.csv")
#插补IBD的PRS###
prsIBD<-fread("~/fatigue_IBD/data/ibd_prs.0.001.profile")
prsIBD<-prsIBD[,c(1,6)]
names(prsIBD)<-c("eid","prs")
prsIBD<-left_join(eid,prsIBD)
prsIBD<-inner_join(covariate,prsIBD)
prsIBD<-inner_join(prsIBD,addCova)
prsIBD<-inner_join(prsIBD,Mental)
varsToFactor <- c("hypertension","heartFailure","renalFailure","asthma",
                  "copd","diabetes","medication","vitamin","Ethnicity","house","alcohol","smoking","sex",
                  "education","sleepness","anxiety","depression")
prsIBD[varsToFactor]<-lapply(prsIBD[varsToFactor],factor)
library(mice)
prsIBD_mice<-mice(prsIBD,m=5,maxit=50,meth="pmm",seed=500)

summary(prsIBD_mice)
resultprsIBD<-complete(prsIBD_mice,action=3) 
sum(is.na(resultprsIBD))
resultprsIBD<-resultprsIBD[,c(1,17)]
write_csv(resultprsIBD,"~/fatigue_IBD/data/ImputeIBDprs.csv")
prsIBD<-resultprsIBD

TableFatigueIBD<-inner_join(TableFatigueIBD,prsIBD)

TableFatigueIBD$PRSGroup<-cut(TableFatigueIBD$prs,quantile(TableFatigueIBD$prs,seq(0,1,1/3)) ,labels=FALSE, righe=TRUE,include.lowest = TRUE)
table(TableFatigueIBD$PRSGroup)
varsToFactor <- c("hypertension","heartFailure","renalFailure","asthma",
                  "copd","diabetes","medication","vitamin","Ethnicity","house","alcohol","smoking","sex",
                  "education","sleepness","anxiety","depression","PRSGroup")
TableFatigueIBD[varsToFactor]<-lapply(TableFatigueIBD[varsToFactor],factor)
TableFatigueIBD$fatigue<-factor(TableFatigueIBD$fatigue,levels = c("0","1"))

coxFactor<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                   vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)

coxNum<-coxph(Surv(time,status)~ tiredness+ sex + Ethnicity +smoking + alcohol + age + BMI +
                vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+dementia+MI+
                stroke+copd+diabetes+prsCD+prsUC,
              data=TableFatigueIBD)
summary(coxNum)
testN<-cox.zph(coxNum)
testN

TablePrs1IBD<-filter(TableFatigueIBD,PRSGroup==1)
TablePrs2IBD<-filter(TableFatigueIBD,PRSGroup==2)
TablePrs3IBD<-filter(TableFatigueIBD,PRSGroup==3)
coxPRS1<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma++copd+diabetes,
               data=TablePrs1IBD)
summary(coxPRS1)
testP1<-cox.zph(coxPRS1)
testP1

coxPRS2<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2IBD)
summary(coxPRS2)
testP2<-cox.zph(coxPRS2)
testP2

coxPRS3<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3IBD)
summary(coxPRS3)
testP3<-cox.zph(coxPRS3)
testP3

#UC######
TableUC=TableUC[,-c(24,25)]
TableUC<-inner_join(TableUC,resultPRS)

TableCD=TableCD[,-c(23,24)]
TableCD<-inner_join(TableCD,resultPRS)


TablePrs1UC<-filter(TableUC,PRSGroup==1)
TablePrs2UC<-filter(TableUC,PRSGroup==2)
TablePrs3UC<-filter(TableUC,PRSGroup==3)
coxPRS1<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                                                           vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
                                                         data=TablePrs1UC)
summary(coxPRS1)
testP1<-cox.zph(coxPRS1)
testP1

coxPRS2<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
               data=TablePrs2UC)
summary(coxPRS2)
testP2<-cox.zph(coxPRS2)
testP2

coxPRS3<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
               data=TablePrs3UC)
summary(coxPRS3)
testP3<-cox.zph(coxPRS3)
testP3


TableCD<-inner_join(TableCD,CDprs)

TableCD$PRSGroupCD<-cut(TableCD$prsCD,quantile(TableCD$prsCD,seq(0,1,1/3)) ,labels=FALSE)
TableCD$PRSGroupUC<-cut(TableCD$prsCD,quantile(TableCD$prsUC,seq(0,1,1/3)) ,labels=FALSE)

TablePrs1CD<-filter(TableCD,PRSGroup==1)
TablePrs2CD<-filter(TableCD,PRSGroup==2)
TablePrs3CD<-filter(TableCD,PRSGroup==3)
coxPRS1<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
               data=TablePrs1CD)
summary(coxPRS1)
testP1<-cox.zph(coxPRS1)
testP1

coxPRS2<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
               data=TablePrs2CD)
summary(coxPRS2)
testP2<-cox.zph(coxPRS2)
testP2

coxPRS3<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI +
                 vegetable + medication + vitamin + house + fruit + MET + TDindex  +sleepness,
               data=TablePrs3CD)
summary(coxPRS3)
testP3<-cox.zph(coxPRS3)
testP3
 
TableFatigueIBD<-filter(TableFatigueIBD,endTime >= "2020-12-01")
ee<-filter(TableFatigueIBD,status==1)

EasyFatigue<-read_csv("data/Fatigue_participant.csv")
EasyFatigue=EasyFatigue[,c(1,3)]
EasyFatigue=na.omit(EasyFatigue)

names(EasyFatigue)<-c("eid","sufferFatigue")

EasyFatigue$sufferFatigue<-as.factor(EasyFatigue$sufferFatigue)
EasyFatigue<-filter(EasyFatigue,sufferFatigue==1|sufferFatigue==0)
EasyFatigue=na.omit(EasyFatigue)
TableNewIBD<-inner_join(TableFatigueIBD,EasyFatigue)
TableNewIBD$sufferFatigue<-factor(TableNewIBD$sufferFatigue,levels = c("0","1"))
TableNewIBD[c("timeNew")]<- as.numeric(difftime(TableNewIBD$endTime,"2020-12-01",
                                         units = c("days"))/365.25)
TableNewIBD<-filter(TableNewIBD,endTime >= "2020-12-01" & deathDate > "2020-12-01")
coxFactor<-coxph(Surv(timeNew,status)~ sufferFatigue+ sex + Ethnicity +smoking + alcohol + age + BMI +
                   vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness,
                 data=TableNewIBD)
summary(coxFactor)
testF<-cox.zph(coxFactor)
testF


#画图####
library(survminer)
TableFatigueIBD<-na.omit(read_csv("TableIBDnew.csv"))

#PSM
PsmEID<-read_csv("/mnt3/mayuying/fatigue_IBD/match_IBD_eid.csv")
TableFatigueIBD<-inner_join(TableFatigueIBD,PsmEID,by="eid")
TableUC<-inner_join(TableUC,PsmEID,by="eid")
TableCD<-inner_join(TableCD,PsmEID,by="eid")
ggsurvplot(survfit(Surv(time,statusUC)~ fatigue, data = TableUC),
           pval = TRUE, conf.int = FALSE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B"),
           ylim=c(0.985,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

ggsurvplot(survfit(Surv(time,statusCD)~ fatigue, data = TableCD),
           pval = TRUE, conf.int = FALSE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B"),
           ylim=c(0.994,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)
ggsurvplot(survfit(Surv(time,status)~ fatigue, data = TableFatigueIBD),
           pval = TRUE, conf.int = FALSE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B"),
           ylim=c(0.985,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

ggsurvplot(survfit(Surv(time,status)~ tiredness, data = TableFatigueIBD),
           pval = TRUE, conf.int = FALSE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B","#42B540","#0099B3"),
           ylim=c(0.985,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

ggsurvplot(survfit(Surv(time,statusCD)~ tiredness, data = TableCD),
           pval = TRUE, conf.int = FALSE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B","#42B540","#0099B3"),
           ylim=c(0.994,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

ggsurvplot(survfit(Surv(time,statusUC)~ tiredness, data = TableUC),
           pval = TRUE, conf.int = FALSE,# 增加置信区间
           pval.coord = c(9,0.8),
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           risk.table = TRUE, # 添加风险表
           risk.table.height = 0.3, # 风险表高度
           risk.table.fontsize = 2.8, # 风险表字体大小
           risk.table.col = "strata", # 根据分层更改风险表颜色
           tables.y.text.col = TRUE,
           tables.y.text = FALSE,
           xlab = "Years",palette =c("#08519C","#B2182B","#42B540","#0099B3"),
           ylim=c(0.985,1), xlim=c(0,16),
           break.y.by=0.005,break.x.by=2,
           censor.size=1.0)

#IBD#####
varsToFactor <- c("sex","Ethnicity" ,"smoking", "alcohol", "medication" , "vitamin" ,
                  "house" , "sleepness","anxiety","depression","hypertension","heartFailure",
                  "renalFailure","asthma","copd","diabetes","PRSGroup")
TableFatigueIBD[varsToFactor]<-lapply(TableFatigueIBD[varsToFactor],factor)
TableFatigueIBD$fatigue<-factor(TableFatigueIBD$fatigue,levels = c("0","1"))

coxFactor<-coxph(Surv(time,status)~ fatigue ,
                 data=TableFatigueIBD)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI, 
                 data=TableFatigueIBD)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue+ sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup, 
                 data=TableFatigueIBD)
summary(coxFactor)
testF<-cox.zph(coxFactor)

testF

TableFatigueIBD$tiredness<-factor(TableFatigueIBD$tiredness,levels = c("1","2","3","4"))

coxNum<-coxph(Surv(time,status)~ tiredness ,
              data=TableFatigueIBD)
summary(coxNum)
coxNum<-coxph(Surv(time,status)~ tiredness + sex + Ethnicity +smoking + alcohol + age + BMI ,
              data=TableFatigueIBD)
summary(coxNum)
coxNum<-coxph(Surv(time,status)~ tiredness+ sex + Ethnicity +smoking + alcohol + age +
                BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
              data=TableFatigueIBD)
summary(coxNum)
testN<-cox.zph(coxNum)
testN

#UC####
varsToFactor <- c("sex","Ethnicity" ,"smoking", "alcohol", "medication" , "vitamin" ,
                  "house" , "sleepness","anxiety","depression","hypertension","heartFailure",
                  "renalFailure","asthma","copd","diabetes","PRSGroupUC")
TableUC[varsToFactor]<-lapply(TableUC[varsToFactor],factor)
TableUC$fatigue<-factor(TableUC$fatigue,levels = c("0","1"))
TableUC$tiredness<-factor(TableUC$tiredness,levels = c("0","1","2","3"))
coxFactor<-coxph(Surv(time,statusUC)~ fatigue ,
                 data=TableUC)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI ,
                                   data=TableUC)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)
testF<-cox.zph(coxFactor)
testF

coxNum<-coxph(Surv(time,statusUC)~ tiredness ,
              data=TableUC)
summary(coxNum)
coxNum<-coxph(Surv(time,statusUC)~ tiredness + sex + Ethnicity +smoking + alcohol + age + BMI ,
              data=TableUC)
summary(coxNum)
coxNum<-coxph(Surv(time,statusUC)~ tiredness+ sex + Ethnicity +smoking + alcohol + age +
                BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
              data=TableUC)
summary(coxNum)

#CD####
varsToFactor <- c("sex","Ethnicity" ,"smoking", "alcohol", "medication" , "vitamin" ,
                  "house" , "sleepness","anxiety","depression","hypertension","heartFailure",
                  "renalFailure","asthma","copd","diabetes","PRSGroupCD")
TableCD[varsToFactor]<-lapply(TableCD[varsToFactor],factor)
TableCD$fatigue<-factor(TableCD$fatigue,levels = c("0","1"))
TableCD$tiredness<-factor(TableCD$tiredness,levels = c("0","1","2","3"))

coxFactor<-coxph(Surv(time,statusCD)~ fatigue ,
                 data=TableCD)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age + BMI ,
                 data=TableCD)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)
testF<-cox.zph(coxFactor)
testF

coxNum<-coxph(Surv(time,statusCD)~ tiredness ,
              data=TableCD)
summary(coxNum)
coxNum<-coxph(Surv(time,statusCD)~ tiredness + sex + Ethnicity +smoking + alcohol + age + BMI ,
              data=TableCD)
summary(coxNum)
coxNum<-coxph(Surv(time,statusCD)~ tiredness + sex + Ethnicity +smoking + alcohol + age +
                BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
              data=TableCD)
summary(coxNum)
#tableone####
library(tableone)
CreateTableOne(data=TableFatigueIBD)

table(TableFatigueIBD$PRSGroupCD)

prop.table(table(TableFatigueIBD$hypertension))
prop.table(table(TableFatigueIBD$diabetes))
prop.table(table(TableFatigueIBD$heartFailure))
prop.table(table(TableFatigueIBD$renalFailure))
prop.table(table(TableFatigueIBD$depression))
prop.table(table(TableFatigueIBD$anxiety))
prop.table(table(TableFatigueIBD$copd))
prop.table(table(TableFatigueIBD$asthma))
prop.table(table(TableFatigueIBD$dementia))
prop.table(table(TableFatigueIBD$stroke))
prop.table(table(TableFatigueIBD$medication))
prop.table(table(TableFatigueIBD$PRSGroup))
prop.table(table(TableFatigueIBD$PRSGroupCD))
prop.table(table(TableFatigueIBD$PRSGroupUC))

TableNo<-filter(TableFatigueIBD,fatigue == 0)
TableYES<-filter(TableFatigueIBD,fatigue==1)
CreateTableOne(data=TableNo)
TableNo[varsToFactor]<-lapply(TableNo[varsToFactor],factor)
summary(TableNo)
table(TableNo$PRSGroupCD)
table(TableNo$PRSGroupUC)
prop.table(table(TableNo$sex))
prop.table(table(TableNo$Ethnicity))
prop.table(table(TableNo$house))
prop.table(table(TableNo$vitamin))
prop.table(table(TableNo$sleepness))
prop.table(table(TableNo$alcohol))
prop.table(table(TableNo$smoking))
prop.table(table(TableNo$hypertension))
prop.table(table(TableNo$diabetes))
prop.table(table(TableNo$heartFailure))
prop.table(table(TableNo$renalFailure))
prop.table(table(TableNo$depression))
prop.table(table(TableNo$anxiety))
prop.table(table(TableNo$copd))
prop.table(table(TableNo$asthma))

prop.table(table(TableNo$medication))
prop.table(table(TableNo$PRSGroup))
prop.table(table(TableNo$PRSGroupUC))
prop.table(table(TableNo$PRSGroupCD))
CreateTableOne(data=TableYES)
TableYES[varsToFactor]<-lapply(TableYES[varsToFactor],factor)
table(TableYES$PRSGroupCD)
table(TableYES$PRSGroupUC)
prop.table(table(TableYES$sex))
prop.table(table(TableYES$Ethnicity))
prop.table(table(TableYES$house))
prop.table(table(TableYES$sleepness))
prop.table(table(TableYES$alcohol))
prop.table(table(TableYES$smoking))
prop.table(table(TableYES$hypertension))
prop.table(table(TableYES$diabetes))
prop.table(table(TableYES$heartFailure))
prop.table(table(TableYES$renalFailure))
prop.table(table(TableYES$depression))
prop.table(table(TableYES$anxiety))
prop.table(table(TableYES$vitamin))
prop.table(table(TableYES$copd))
prop.table(table(TableYES$asthma))
prop.table(table(TableYES$dementia))
prop.table(table(TableYES$stroke))
prop.table(table(TableYES$medication))
prop.table(table(TableYES$PRSGroup))
prop.table(table(TableYES$PRSGroupUC))
prop.table(table(TableYES$PRSGroupCD))


#PRS

TablePrs1UC<-filter(TableUC,PRSGroupUC==1)
TablePrs2UC<-filter(TableUC,PRSGroupUC==2)
TablePrs3UC<-filter(TableUC,PRSGroupUC==3)
coxPRS1<-coxph(Surv(time,statusUC)~ fatigue+sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs1UC)
summary(coxPRS1)
testP1<-cox.zph(coxPRS1)
testP1

coxPRS2<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2UC)
summary(coxPRS2)
testP2<-cox.zph(coxPRS2)
testP2

coxPRS3<-coxph(Surv(time,statusUC)~ fatigue + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3UC)
summary(coxPRS3)
testP3<-cox.zph(coxPRS3)
testP3

#NUM
coxPRS1<-coxph(Surv(time,statusUC)~ tiredness+sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs1UC)
summary(coxPRS1)


coxPRS2<-coxph(Surv(time,statusUC)~ tiredness + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2UC)
summary(coxPRS2)


coxPRS3<-coxph(Surv(time,statusUC)~ tiredness + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3UC)
summary(coxPRS3)
#CD####
TableCD<-inner_join(TableCD,CDprs)

TableCD$PRSGroupCD<-cut(TableCD$prsCD,quantile(TableCD$prsCD,seq(0,1,1/3)) ,labels=FALSE)
TableCD$PRSGroupUC<-cut(TableCD$prsCD,quantile(TableCD$prsUC,seq(0,1,1/3)) ,labels=FALSE)

TablePrs1CD<-filter(TableCD,PRSGroupCD==1)
TablePrs2CD<-filter(TableCD,PRSGroupCD==2)
TablePrs3CD<-filter(TableCD,PRSGroupCD==3)
coxPRS1<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs1CD)
summary(coxPRS1)

coxPRS2<-coxph(Surv(time,statusCD)~ fatigue + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2CD)
summary(coxPRS2)


coxPRS3<-coxph(Surv(time,statusCD)~ fatigue  + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3CD)
summary(coxPRS3)

coxPRS1<-coxph(Surv(time,statusCD)~ tiredness + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs1CD)
summary(coxPRS1)

coxPRS2<-coxph(Surv(time,statusCD)~ tiredness+ sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2CD)
summary(coxPRS2)


coxPRS3<-coxph(Surv(time,statusCD)~ tiredness  + sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3CD)
summary(coxPRS3)

TablePrs1IBD<-filter(TableFatigueIBD,PRSGroup==1)
TablePrs2IBD<-filter(TableFatigueIBD,PRSGroup==2)
TablePrs3IBD<-filter(TableFatigueIBD,PRSGroup==3)

coxPRS1<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs1IBD)
summary(coxPRS1)


coxPRS2<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2IBD)
summary(coxPRS2)


coxPRS3<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3IBD)
summary(coxPRS3)

coxPRS1<-coxph(Surv(time,status)~ tiredness +sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs1IBD)
summary(coxPRS1)


coxPRS2<-coxph(Surv(time,status)~ tiredness +sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs2IBD)
summary(coxPRS2)


coxPRS3<-coxph(Surv(time,status)~ tiredness +sex + Ethnicity +smoking + alcohol + age +
                 BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes,
               data=TablePrs3IBD)
summary(coxPRS3)

#三、亚组####

#IBD####

#年龄#
tableYoung<-filter(TableFatigueIBD,age<=60)
tableOld<-filter(TableFatigueIBD,age>60)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + 
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableYoung)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol  +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableOld)
summary(coxFactor)
#性别#
tableFemale<-filter(TableFatigueIBD,sex==0)
tableMale<-filter(TableFatigueIBD,sex==1)
coxFactor<-coxph(Surv(time,status)~ fatigue  + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableFemale)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMale)
summary(coxFactor)

#TD
mean(TableFatigueIBD$TDindex)
tableMin<-filter(TableFatigueIBD,`TDindex`<=-1.325645)
tableMax<-filter(TableFatigueIBD,`TDindex`>-1.325645)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#BMI
tableMin<-filter(TableFatigueIBD, BMI <= 30)
tableMax<-filter(TableFatigueIBD, BMI > 30)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                 + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                  + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#alcohol
tableMin<-filter(TableFatigueIBD, alcohol == 1 | alcohol == 2)
tableMax<-filter(TableFatigueIBD, alcohol == 3 | alcohol == 4| alcohol == 5 |alcohol == 6)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#smoking
tableMin<-filter(TableFatigueIBD, smoking == 0 )
tableMedian<-filter(TableFatigueIBD, smoking == 1)
tableMax<-filter(TableFatigueIBD, smoking == 2)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity  + alcohol + age + BMI +
                   vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+dementia+MI+
                   stroke+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity  + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMedian)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity+ alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)


#MET
mean(TableFatigueIBD$MET)
tableMin<-filter(TableFatigueIBD,MET<=2644.682)
tableMax<-filter(TableFatigueIBD,MET>2644.682)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit  + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit  + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#Vita#
tableMin<-filter(TableFatigueIBD,vitamin==0)
tableMax<-filter(TableFatigueIBD,vitamin==1)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication  + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication  + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#Hypertension#
tableMin<-filter(TableFatigueIBD,hypertension==0)
tableMax<-filter(TableFatigueIBD,hypertension==1)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#diabetes#
tableMin<-filter(TableFatigueIBD,diabetes==0)
tableMax<-filter(TableFatigueIBD,diabetes==1)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#depression#
tableMin<-filter(TableFatigueIBD,depression==0)
tableMax<-filter(TableFatigueIBD,depression==1)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#anxiety#
tableMin<-filter(TableFatigueIBD,anxiety==0)
tableMax<-filter(TableFatigueIBD,anxiety==1)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                  depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)


#medication#
tableMin<-filter(TableFatigueIBD,medication==0)
tableMax<-filter(TableFatigueIBD,medication==1)
coxFactor<-coxph(Surv(time,status)~ fatigue ++sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable  + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,status)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable  + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMax)
summary(coxFactor)

#UC####
#年龄#
tableYoung<-filter(TableUC,age<=60)
tableOld<-filter(TableUC,age>60)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol  +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableYoung)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol  +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableOld)
summary(coxFactor)
#性别#
tableFemale<-filter(TableUC,sex==0)
tableMale<-filter(TableUC,sex==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableFemale)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMale)
summary(coxFactor)

#TD
mean(TableUC$TDindex)
tableMin<-filter(TableUC,`TDindex`<=-1.3268)
tableMax<-filter(TableUC,`TDindex`>-1.3268)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#BMI
tableMin<-filter(TableUC, BMI <= 30)
tableMax<-filter(TableUC, BMI > 30)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                 vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#alcohol
tableMin<-filter(TableUC, alcohol == 1 | alcohol == 2)
tableMax<-filter(TableUC, alcohol == 3 | alcohol == 4| alcohol == 5 |alcohol == 6)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking  + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#smoking
tableMin<-filter(TableUC, smoking == 0 )
tableMedian<-filter(TableUC, smoking == 1)
tableMax<-filter(TableUC, smoking == 2)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMedian)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity  + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)


#MET
mean(TableUC$MET)
tableMin<-filter(TableUC,MET<=2645.111)
tableMax<-filter(TableUC,MET>2645.111)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit  + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#Vita#
tableMin<-filter(TableUC,vitamin==0)
tableMax<-filter(TableUC,vitamin==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication  + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication  + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#Hypertension#
tableMin<-filter(TableUC,hypertension==0)
tableMax<-filter(TableUC,hypertension==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#diabetes#
tableMin<-filter(TableUC,diabetes==0)
tableMax<-filter(TableUC,diabetes==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#depression#
tableMin<-filter(TableUC,depression==0)
tableMax<-filter(TableUC,depression==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#anxiety#
tableMin<-filter(TableUC,anxiety==0)
tableMax<-filter(TableUC,anxiety==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue+sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                 depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)



#medication#
tableMin<-filter(TableUC,medication==0)
tableMax<-filter(TableUC,medication==1)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMax)
summary(coxFactor)

#CD####
#年龄#
tableYoung<-filter(TableCD,age<=60)
tableOld<-filter(TableCD,age>60)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol +
                   BMI + vegetable+medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableYoung)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableOld)
summary(coxFactor)
#性别#
tableFemale<-filter(TableCD,sex==0)
tableMale<-filter(TableCD,sex==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue+ Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableFemale)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMale)
summary(coxFactor)

#TD
mean(TableCD$TDindex)
tableMin<-filter(TableCD,`TDindex`<=-1.327296)
tableMax<-filter(TableCD,`TDindex`>-1.327296)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#BMI
tableMin<-filter(TableCD, BMI <= 30)
tableMax<-filter(TableCD, BMI > 30)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#alcohol
tableMin<-filter(TableCD, alcohol == 1 | alcohol == 2)
tableMax<-filter(TableCD, alcohol == 3 | alcohol == 4| alcohol == 5 |alcohol == 6)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking  + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking  + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#smoking
tableMin<-filter(TableCD, smoking == 0 )
tableMedian<-filter(TableCD, smoking == 1)
tableMax<-filter(TableCD, smoking == 2)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMedian)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)


#MET
mean(TableCD$MET)
tableMin<-filter(TableCD,MET<=2644.475)
tableMax<-filter(TableCD,MET>2644.475)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit  + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

#Hypertension#
tableMin<-filter(TableCD,hypertension==0)
tableMax<-filter(TableCD,hypertension==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#diabetes#
tableMin<-filter(TableCD,diabetes==0)
tableMax<-filter(TableCD,diabetes==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#depression#
tableMin<-filter(TableCD,depression==0)
tableMax<-filter(TableCD,depression==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue+sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#anxiety#
tableMin<-filter(TableCD,anxiety==0)
tableMax<-filter(TableCD,anxiety==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable +medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)


#medication#
tableMin<-filter(TableCD,medication==0)
tableMax<-filter(TableCD,medication==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMin)
summary(coxFactor)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMax)
summary(coxFactor)

#四、交互####
#CD####

#年龄#
tableYoung<-filter(TableCD,age<=60)
tableYoung$A<-as.factor(0:0)
tableOld<-filter(TableCD,age>60)
tableOld$A<-as.factor(1:1)

TableAge<-full_join(tableOld,tableYoung)
TableAge$A<-factor(TableAge$A,levels = c("0","1"))
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*A+A+fatigue +sex + Ethnicity +smoking + alcohol +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableAge)
summary(coxFactor)

#性别#
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*sex +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)

#TD
mean(TableCD$TDindex)
tableMin<-filter(TableCD,`TDindex`<=-1.327296)
tableMin$A<-as.factor(0:0)
tableMax<-filter(TableCD,`TDindex`>-1.327296)
tableMax$A<-as.factor(1:1)
tableTD<-full_join(tableMax,tableMin)
tableTD$A<-factor(tableTD$A,levels = c("0","1"))
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET  +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableTD)

summary(coxFactor)

#BMI
tableMin<-filter(TableCD, BMI <= 30)
tableMin$A<-as.factor(0:0)
tableMax<-filter(TableCD, BMI > 30)
tableMax$A<-as.factor(1:1)
tableBMI<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                 medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableBMI)
summary(coxFactor)


#alcohol
tableMin<-filter(TableCD, alcohol == 1 | alcohol == 2)
tableMax<-filter(TableCD, alcohol == 3 | alcohol == 4| alcohol == 5 |alcohol == 6)
tableMin$A<-as.factor(0:0)
tableMax$A<-as.factor(1:1)
tableAL<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*A +sex + Ethnicity +smoking + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableAL)
summary(coxFactor)


#smoking
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*smoking +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)


#MET
mean(TableCD$MET)
tableMin<-filter(TableCD,MET<=2644.475)
tableMax<-filter(TableCD,MET>2644.475)
tableMin$A<-as.factor(0:0)
tableMax$A<-as.factor(1:1)
tableMET<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=tableMET)
summary(coxFactor)
#Vita#

coxFactor<-coxph(Surv(time,statusCD)~ fatigue*vitamin +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)


#Hypertension#

coxFactor<-coxph(Surv(time,statusCD)~ fatigue*hypertension +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)

#diabetes#
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*diabetes +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)


#depression#
tableMin<-filter(TableCD,depression==0)
tableMax<-filter(TableCD,depression==1)
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*depression +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)

#anxiety#
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*anxiety +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)




#medication#
coxFactor<-coxph(Surv(time,statusCD)~ fatigue*medication +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
                 data=TableCD)
summary(coxFactor)


#UC####

#年龄#
tableYoung<-filter(TableUC,age<=60)
tableYoung$A<-as.factor(0:0)
tableOld<-filter(TableUC,age>60)
tableOld$A<-as.factor(1:1)

TableAge<-full_join(tableOld,tableYoung)
TableAge$A<-factor(TableAge$A,levels = c("0","1"))
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*A+A+fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableAge)
summary(coxFactor)

#性别#
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*sex +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)

#TD
mean(TableUC$TDindex)
tableMin<-filter(TableUC,`TDindex`<=-1.3268)
tableMin$A<-as.factor(0:0)
tableMax<-filter(TableUC,`TDindex`>-1.3268)
tableMax$A<-as.factor(1:1)
tableTD<-full_join(tableMax,tableMin)
tableTD$A<-factor(tableTD$A,levels = c("0","1"))
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableTD)

summary(coxFactor)

#BMI
tableMin<-filter(TableUC, BMI <= 30)
tableMin$A<-as.factor(0:0)
tableMax<-filter(TableUC, BMI > 30)
tableMax$A<-as.factor(1:1)
tableBMI<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableBMI)
summary(coxFactor)


#alcohol
tableMin<-filter(TableUC, alcohol == 1 | alcohol == 2)
tableMax<-filter(TableUC, alcohol == 3 | alcohol == 4| alcohol == 5 |alcohol == 6)
tableMin$A<-as.factor(0:0)
tableMax$A<-as.factor(1:1)
tableAL<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableAL)
summary(coxFactor)


#smoking
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*smoking +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)


#MET
mean(TableUC$MET)
tableMin<-filter(TableUC,MET<=2645.111)
tableMax<-filter(TableUC,MET>2645.111)
tableMin$A<-as.factor(0:0)
tableMax$A<-as.factor(1:1)
tableMET<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=tableMET)
summary(coxFactor)
#Vita#

coxFactor<-coxph(Surv(time,statusUC)~ fatigue*vitamin+sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)


#Hypertension#

coxFactor<-coxph(Surv(time,statusUC)~ fatigue*hypertension +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)

#diabetes#
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*diabetes +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)


#depression#

coxFactor<-coxph(Surv(time,statusUC)~ fatigue*depression +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)

#anxiety#
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*anxiety +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)




#medication#
coxFactor<-coxph(Surv(time,statusUC)~ fatigue*medication +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
                 data=TableUC)
summary(coxFactor)


#IBD####

#年龄#
tableYoung<-filter(TableFatigueIBD,age<=60)
tableYoung$A<-as.factor(0:0)
tableOld<-filter(TableFatigueIBD,age>60)
tableOld$A<-as.factor(1:1)

TableAge<-full_join(tableOld,tableYoung)
TableAge$A<-factor(TableAge$A,levels = c("0","1"))
coxFactor<-coxph(Surv(time,status)~ fatigue*A+A+fatigue +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableAge)
summary(coxFactor)

#性别#
coxFactor<-coxph(Surv(time,status)~ fatigue*sex +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)

#TD
mean(TableFatigueIBD$TDindex)
tableMin<-filter(TableFatigueIBD,`TDindex`<= -1.325645)
tableMin$A<-as.factor(0:0)
tableMax<-filter(TableFatigueIBD,`TDindex`> -1.325645)
tableMax$A<-as.factor(1:1)
tableTD<-full_join(tableMax,tableMin)
tableTD$A<-factor(tableTD$A,levels = c("0","1"))
coxFactor<-coxph(Surv(time,status)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableTD)

summary(coxFactor)

#BMI
tableMin<-filter(TableFatigueIBD, BMI <= 30)
tableMin$A<-as.factor(0:0)
tableMax<-filter(TableFatigueIBD, BMI > 30)
tableMax$A<-as.factor(1:1)
tableBMI<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,status)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableBMI)
summary(coxFactor)


#alcohol
tableMin<-filter(TableFatigueIBD, alcohol == 1 | alcohol == 2)
tableMax<-filter(TableFatigueIBD, alcohol == 3 | alcohol == 4| alcohol == 5 |alcohol == 6)
tableMin$A<-as.factor(0:0)
tableMax$A<-as.factor(1:1)
tableAL<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,status)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableAL)
summary(coxFactor)


#smoking
coxFactor<-coxph(Surv(time,status)~ fatigue*smoking +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)


#MET
mean(TableFatigueIBD$MET)
tableMin<-filter(TableFatigueIBD,MET<=2644.682)
tableMax<-filter(TableFatigueIBD,MET>2644.682)
tableMin$A<-as.factor(0:0)
tableMax$A<-as.factor(1:1)
tableMET<-full_join(tableMax,tableMin)
coxFactor<-coxph(Surv(time,status)~ fatigue*A +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=tableMET)
summary(coxFactor)
#Vita#

coxFactor<-coxph(Surv(time,status)~ fatigue*vitamin +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)


#Hypertension#

coxFactor<-coxph(Surv(time,status)~ fatigue*hypertension +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)

#diabetes#
coxFactor<-coxph(Surv(time,status)~ fatigue*diabetes +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)


#depression#

coxFactor<-coxph(Surv(time,status)~ fatigue*depression +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)

#anxiety#
coxFactor<-coxph(Surv(time,status)~ fatigue*anxiety +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)

#medication#
coxFactor<-coxph(Surv(time,status)~ fatigue*medication +sex + Ethnicity +smoking + alcohol + age +
                   BMI +medication +vegetable+ vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
                 data=TableFatigueIBD)
summary(coxFactor)


#PRS作图####
library(magrittr)
phenotype <- fread("Pheno.IBD")
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
covariate <- fread("~/GWAS/covaGrip.csv")
pheno <- merge(phenotype, covariate)
null.r2 <- summary(lm(IBD~., data=pheno[,-c("FID", "IID")]))$r.squared
prs.result <- NULL
setwd("data/")
for(i in p.threshold){
  pheno.prs <- paste0("ibd_prs.", i, ".profile") %>%
    fread(.) %>%
    .[,c("FID", "IID", "SCORE")] %>%
    merge(., pheno, by=c("FID", "IID"))
  model <- lm(IBD~., data=pheno.prs[,-c("FID","IID")]) %>%
    summary
  model.r2 <- model$r.squared
  prs.r2 <- model.r2-null.r2
  prs.coef <- model$coeff["SCORE",]
  prs.result %<>% rbind(.,
                        data.frame(Threshold=i, R2=prs.r2, 
                                   P=as.numeric(prs.coef[4]), 
                                   BETA=as.numeric(prs.coef[1]),
                                   SE=as.numeric(prs.coef[2])))
}
print(prs.result[which.max(prs.result$R2),])
library(ggplot2)

prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 4,
    parse = T
  )  +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)
  ) +
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size =
                                  18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust =
                                 1)
  )

#新冠后分析####
IBDafter<-filter(TableFatigueIBD,status== 1 & endTime>="2020-01-01")
IBDcontrol<-filter(TableFatigueIBD,status==0 &endTime>="2020-01-01")
afterIBD<-full_join(IBDafter,IBDcontrol)
afterIBD$newBegin<-as.Date("2020-01-01")
afterIBD[c("timeN")]<- as.numeric(difftime(afterIBD$endTime,afterIBD$newBegin,
                                                 units = c("days"))/365.25)
afterIBD[varsToFactor]<-lapply(afterIBD[varsToFactor],factor)
afterIBD$fatigue<-factor(afterIBD$fatigue,levels = c("0","1"))

coxFactor<-coxph(Surv(timeN,status)~ fatigue+ sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup, 
                 data=afterIBD)
summary(coxFactor)

afterIBD$tiredness<-factor(afterIBD$tiredness,levels = c("1","2","3","4"))

coxNum<-coxph(Surv(timeN,status)~ tiredness+ sex + Ethnicity +smoking + alcohol + age +
                BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroup,
              data=afterIBD)
summary(coxNum)
testN<-cox.zph(coxNum)
testN

#UC
UCafter<-filter(TableUC,statusUC== 1 & endTime>="2020-01-01")
UCcontrol<-filter(TableUC,statusUC==0 &endTime>="2020-01-01")
afterUC<-full_join(UCafter,UCcontrol)
afterUC$newBegin<-as.Date("2020-01-01")
afterUC[c("timeN")]<- as.numeric(difftime(afterUC$endTime,afterUC$newBegin,
                                           units = c("days"))/365.25)
afterUC[varsToFactor]<-lapply(afterUC[varsToFactor],factor)
afterUC$fatigue<-factor(afterUC$fatigue,levels = c("0","1"))

coxFactor<-coxph(Surv(timeN,statusUC)~ fatigue+ sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC, 
                 data=afterUC)
summary(coxFactor)

afterUC$tiredness<-factor(afterUC$tiredness,levels = c("0","1","2","3"))

coxNum<-coxph(Surv(timeN,statusUC)~ tiredness+ sex + Ethnicity +smoking + alcohol + age +
                BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupUC,
              data=afterUC)
summary(coxNum)
testN<-cox.zph(coxNum)
testN

#CD
CDafter<-filter(TableCD,statusCD== 1 & endTime>="2020-01-01")
CDcontrol<-filter(TableCD,statusCD==0 &endTime>="2020-01-01")
afterCD<-full_join(CDafter,CDcontrol)
afterCD$newBegin<-as.Date("2020-01-01")
afterCD[c("timeN")]<- as.numeric(difftime(afterCD$endTime,afterCD$newBegin,
                                          units = c("days"))/365.25)
afterCD[varsToFactor]<-lapply(afterCD[varsToFactor],factor)
afterCD$fatigue<-factor(afterCD$fatigue,levels = c("0","1"))

coxFactor<-coxph(Surv(timeN,statusCD)~ fatigue+ sex + Ethnicity +smoking + alcohol + age +
                   BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                   anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD, 
                 data=afterCD)
summary(coxFactor)

afterCD$tiredness<-factor(afterCD$tiredness,levels = c("0","1","2","3"))

coxNum<-coxph(Surv(timeN,statusCD)~ tiredness+ sex + Ethnicity +smoking + alcohol + age +
                BMI + vegetable + medication + vitamin + house + fruit + MET + TDindex +sleepness+
                anxiety+depression+hypertension+heartFailure+renalFailure+asthma+copd+diabetes+PRSGroupCD,
              data=afterCD)
summary(coxNum)
testN<-cox.zph(coxNum)
testN
