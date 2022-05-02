## COMPLETE R CODE ##

setwd("C:/Users/rona/Dropbox/Rona/dissertation")

    ### LOAD PACKAGES ###
library(climwin)
library(lubridate)
library(readxl)
library(ncdf4) 
library(raster)
library(ggplot2)
library(writexl)
library(lme4)
library(lmerTest)
library(cowplot)
library(dplyr)

    ### IMPORT ALL DATA ###
#population size 
pop<-read_excel("data/population_1a.xlsx")

#variance in egg laying (2018-2021)
phenology<-read_excel("data/phenology.xlsx")

#date of first egg in each year
fed<-read_excel("data/egg_phenology.xlsx")

#fed in climwin format
egg<-read_excel("data/egg.xlsx")

#mean annual productivity (chicks fledged/breeding pair)
productivity<-read_excel("data/productivity.xlsx")

#mean annual chick diet composition
diet<-read_excel("data/diet.xlsx")

#sandeel fisheries data
ottertrawl<-read.csv("data/ottertrawl.csv")



    ### IDENTIFYING CLIMATE WINDOW ###
## DATA PREPARATION ##
#sst data
sstdata<-nc_open("data/sst2.nc")

#read data from opened file, store as sst
sst<-ncvar_get(sstdata, "analysed_sst")

#get mean value for each timepoint from 4 lat-long squares and convert to celsius
abszeroK<-273.15 ##set absolute zero in kelvin
mean_sst<-apply(sst,3,mean,na.rm=T)-abszeroK
sstdata$mean_sst<-mean_sst #store in sstdata

#create dataframe and export to excel
df <- data.frame(mean_sst = mean_sst)
write_xlsx(df,"data/sst2.xlsx")

## SLIDING WINDOW APPROACH ##

#make sliding window model
eggwin <- slidingwin(xvar = list(Temp = sst$mean_sst),
                     cdate = sst$date,
                     bdate = egg$egg_date,
                     baseline = lm(day_of_year ~ 1, data = egg),
                     cinterval = "week",
                     range = c(20, 0),
                     type = "absolute", refday = c(01,05),
                     stat = "mean",
                     func = "lin")
#find best window
head(eggwin[[1]]$Dataset) 


#visualise windows
summarytab<-eggwin[[1]]$Dataset
plot(NULL,NULL,xlim=c(0,20),ylim=c(-0.5,2),xlab="Weeks before reference date",
     ylab=expression(paste(Delta,"AICc")),las=1)
for(x in 1:dim(summarytab)[1]){
  points(c(summarytab[x,"WindowOpen"],summarytab[x,"WindowClose"]),
         rep(summarytab[x,"deltaAICc"],2),type="l")}
points(c(summarytab[1,"WindowOpen"],summarytab[1,"WindowClose"]),
       rep(summarytab[1,"deltaAICc"],2),type="l",col="red",lwd=8)



#relationship between temp & date in strongest climate window
eggwin[[1]]$BestModel
head(eggwin[[1]]$BestModelData)

#account fo overfitting
#use randwin to re-run 5000 climwin analyses with true signal removed
eggrand <- randwin(repeats = 5000, 
                   xvar = list(Temp = sst$mean_sst),
                   cdate = sst$date,
                   bdate = egg$egg_date,
                   baseline = lm(day_of_year ~ 1, data = egg),
                   cinterval = "week",
                   range = c(20, 0),
                   type = "absolute", refday = c(01,05),
                   stat = "mean",
                   func = "lin")
eggrand[[1]] #store output

#test null hypothesis
pvalue(dataset = eggwin[[1]]$Dataset, datasetrand =
         eggrand[[1]], metric = "C", sample.size = 26)

#visualise climate window models
eggoutput<-eggwin[[1]]$Dataset 
randoutput<-eggrand[[1]]
plothist(dataset=eggoutput,datasetrand=randoutput) #create histogram
plotdelta(dataset=eggoutput) #plot delta AICc values
plotweights(dataset=eggoutput) #plot model weights
plotbetas(dataset=eggoutput) #plot model coefficients
plotwin(dataset=eggoutput) #plot median window as boxplot
eggsingle <- singlewin(xvar = list(Temp = sst$mean_sst),
                       cdate = sst$date,
                       bdate = egg$egg_date,
                       baseline = lm(day_of_year ~ 1, data = egg),
                       cinterval = "day",
                       range = c(7,0),
                       type = "absolute", refday = c(01,05),
                       stat = "mean",
                       func = "lin") #plot best model over bio data
plotbest(dataset = eggoutput,
         bestmodel = eggsingle$BestModel, 
         bestmodeldata = eggsingle$BestModelData) #use output to plot best data

#plot all figures together
plotall(dataset = eggoutput,
        datasetrand = eggrand[[1]],
        bestmodel = eggsingle$BestModel, 
        bestmodeldata = eggsingle$BestModelData)





    ### EFFECT OF POPULATION SIZE ON FED ###
## CALCULATE VARIANCE IN FED FOR 2018, 2019, 2021 ##
#2018 
phenol2018<-subset(phenology,year==2018) #extract 2018 data
#create for loop to record simulated egg laying distribution
dates2018<-c()f
for(x in 1:dim(phenol2018)[1]){
  dates2018<-c(dates2018, rep(
    phenol2018$day_of_year[x],as.numeric(phenol2018$incubating[x])))}
var2018<-var(dates2018) #store variance in egg date
mean2018<-mean(dates2018) #store mean egg date
sd2018<-sd(dates2018) #store std dev egg date

#repeat for 2019 and 2021:
#2019
phenol2019<-subset(phenology,year==2019)
dates2019<-c()
for(x in 1:dim(phenol2019)[1]){
  dates2019<-c(dates2019, rep(
    phenol2019$day_of_year[x],as.numeric(phenol2019$incubating[x])))}
var2019<-var(dates2019)
mean2019<-mean(dates2019)
sd2019<-sd(dates2019)

#2021
phenol2021<-subset(phenology,year==2021)
dates2021<-c()
for(x in 1:dim(phenol2021)[1]){
  dates2021<-c(dates2021, rep(
    phenol2021$day_of_year[x],as.numeric(phenol2021$incubating[x])))}
var2021<-var(dates2021)
mean2021<-mean(dates2021)
sd(dates2021)
sd2021<-sd(dates2021)

#estimate mean var, mean and sd across the 3 years
total_var<-mean(var2018,var2019,var2021)
total_mean<-mean(mean2018,mean2019,mean2021)
meansd<-mean(sd2018,sd2019,sd2021)


#for loop to simulate FED for each year in study period
#based on variance due to population size across 1000 repetitions
#extract slopes and p-value
store.p<-c()
store.slope<-c()

for(x in 1:1000)
{
  null1996<-min(rnorm(round(as.numeric(pop$mean[pop$year==1996])),0,meansd))
  null1997<-min(rnorm(round(as.numeric(pop$mean[pop$year==1997])),0,meansd))
  null1998<-min(rnorm(round(as.numeric(pop$mean[pop$year==1998])),0,meansd))
  null1999<-min(rnorm(round(as.numeric(pop$mean[pop$year==1999])),0,meansd))
  null2000<-min(rnorm(round(as.numeric(pop$mean[pop$year==2000])),0,meansd))
  null2001<-min(rnorm(round(as.numeric(pop$mean[pop$year==2001])),0,meansd))
  null2002<-min(rnorm(round(as.numeric(pop$mean[pop$year==2002])),0,meansd))
  null2003<-min(rnorm(round(as.numeric(pop$mean[pop$year==2003])),0,meansd))
  null2004<-min(rnorm(round(as.numeric(pop$mean[pop$year==2004])),0,meansd))
  null2005<-min(rnorm(round(as.numeric(pop$mean[pop$year==2005])),0,meansd))
  null2006<-min(rnorm(round(as.numeric(pop$mean[pop$year==2006])),0,meansd))
  null2007<-min(rnorm(round(as.numeric(pop$mean[pop$year==2007])),0,meansd))
  null2008<-min(rnorm(round(as.numeric(pop$mean[pop$year==2008])),0,meansd))
  null2009<-min(rnorm(round(as.numeric(pop$mean[pop$year==2009])),0,meansd))
  null2010<-min(rnorm(round(as.numeric(pop$mean[pop$year==2010])),0,meansd))
  null2011<-min(rnorm(round(as.numeric(pop$mean[pop$year==2011])),0,meansd))
  null2012<-min(rnorm(round(as.numeric(pop$mean[pop$year==2012])),0,meansd))
  null2013<-min(rnorm(round(as.numeric(pop$mean[pop$year==2013])),0,meansd))
  null2014<-min(rnorm(round(as.numeric(pop$mean[pop$year==2014])),0,meansd))
  null2015<-min(rnorm(round(as.numeric(pop$mean[pop$year==2015])),0,meansd))
  null2016<-min(rnorm(round(as.numeric(pop$mean[pop$year==2016])),0,meansd))
  null2017<-min(rnorm(round(as.numeric(pop$mean[pop$year==2017])),0,meansd))
  null2018<-min(rnorm(round(as.numeric(pop$mean[pop$year==2018])),0,meansd))
  null2019<-min(rnorm(round(as.numeric(pop$mean[pop$year==2019])),0,meansd))
  null2021<-min(rnorm(round(as.numeric(pop$mean[pop$year==2021])),0,meansd))
  
  nullFED<-c(null1996, null1997,null1998,null1999,null2000,null2001,null2002,
             null2003,null2004,null2005,null2006,null2007,null2008,null2009,
             null2010,null2011,null2012,null2013,null2014,null2015,null2016,
             null2017,null2018,null2019,null2021)
  
  model<-lm(nullFED~pop$year)
  
  #save slope and pvalue
  store.p[x]<-summary(model)$coefficients[2,4]
  store.slope[x]<-summary(model)$coefficients[2,1]
}

hist(store.p) #plot histogram of simulated p-values

length(which(store.p<=0.05))/1000 #calculate type 1 error rate

hist(store.slope) #plot histogram of simulated slopes 
mean(store.slope) #calculate mean simulated slope

## COMPARE SIMULATION TO OBSERVED FEDs ##

#create model of change in FED over time
obsmod<-lm(day_of_year~year,data=fed)
par(mfrow=c(2,2))
plot(obsmod) #check assumptions
par(mfrow=c(1,1))
summary(obsmod) #view model output

#visualise change in FED
obsplot<-ggplot(aes(x=year,y=day_of_year),data=fed)+geom_point()+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Year",y="Day of year of first egg laid")+theme(text = element_text(size=14))
obsplot

#store observed slope from obsmod
cf<-coef(obsmod)
obsslope<-cf[2]

#compare observed change in FED to expectations based on population simulation
#probability of observed shift occurring by chance:
obsprob<-(2*length(which(store.slope>=obsslope)))/10000




    ### EFFECT OF FED ON PRODUCTIVITY ###
## DATA PREPARATION ##
#create dataset with egg date and productivity
eggprod<-data.frame(fed$year,fed$day_of_year,productivity$productivity,pop$mean)
#give sensible column names
names(eggprod)[names(eggprod) == "fed.year"] <- "year"
names(eggprod)[names(eggprod) == "fed.day_of_year"] <- "day_of_year"
names(eggprod)[names(eggprod) == "productivity.productivity"] <- "prod"
names(eggprod)[names(eggprod) == "pop.mean"] <- "pop"

## DATA ANALYSIS ##
#make simple model
prodmod<-lm(prod~day_of_year,data=eggprod)
par(mfrow=c(2,2))
plot(prodmod) #check assumptions of simple model
par(mfrow=c(1,1))
summary(prodmod) #view outputs

#include population in model - control for density dependence
prodpopmod<-lm(prod~pop+day_of_year,data=eggprod)
par(mfrow=c(2,2))
plot(prodpopmod) #check assumptions of density dependence model
par(mfrow=c(1,1))
summary(prodpopmod) #view output

#visualise including outliers
prodplot<-ggplot(aes(x=day_of_year,y=prod),data=eggprod)+geom_point()+
  geom_smooth(method=lm,level=0.95)+theme_classic()+geom_point()+
  labs(x="Day of year of first egg laid",y="Breeding success")+
  theme(text = element_text(size=14))
prodplot

## REPEAT ANALYSIS EXCLUDING THE LOW PROD YEARS ##
#remove 2004, 2011, 2013 from data
prod2<-productivity[-c(9,16,18),]
egg2<-fed[-c(9,16,18),]
pop2<-pop[-c(9,16,18),]

#create dataset with egg date and productivity
eggprod2<-data.frame(egg2$year,egg2$day_of_year,prod2$productivity,pop2$mean)
#give sensible names
names(eggprod2)[names(eggprod2) == "egg2.year"] <- "year"
names(eggprod2)[names(eggprod2) == "egg2.day_of_year"] <- "day_of_year"
names(eggprod2)[names(eggprod2) == "prod2.productivity"] <- "prod"
names(eggprod2)[names(eggprod2) == "pop2.mean"] <- "pop"

#make simple model
prodmod2<-lm(prod~day_of_year,data=eggprod2)
par(mfrow=c(2,2))
plot(prodmod2)#check assumptions of simple model
par(mfrow=c(1,1))
summary(prodmod2) #view output


#include population in model - control for density dependence
prodpopmod2<-lm(prod~pop+day_of_year,data=eggprod2)
par(mfrow=c(2,2))
plot(prodpopmod2) #check assumptions of density dependence model
par(mfrow=c(1,1))
summary(prodpopmod2) #view output

#arrange so model tests if day_of_year sig when pop controlled
#day of year not significant

#visualise excluding outliers
prodplot2<-ggplot(aes(x=day_of_year,y=prod),data=eggprod2)+geom_point()+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Day of year of first egg laid",y="Breeding success")+
  theme(text = element_text(size=14))
prodplot2


#create dataset with productivity columns with and without outliers 
eggprod3<-eggprod
eggprod3$prod2<-eggprod$prod
is.na(eggprod3$prod2)<-eggprod3$year=="2004"
is.na(eggprod3$prod2)<-eggprod3$year=="2011"
is.na(eggprod3$prod2)<-eggprod3$year=="2013"

#plot lines with and without outliers together
combiplot<-ggplot(aes(x=day_of_year),data=eggprod3)+geom_point(aes(y=prod))+
  geom_smooth(method="lm",aes(y=prod),
              color="blue",fill="blue",alpha=.15)+geom_smooth(
                method="lm",aes(y=prod2),
                na.rm=TRUE,color="red",fill="red",alpha=.15)+
  theme_classic()+labs(
    x="Day of year of first egg laid",
    y="Breeding success")+theme(text = element_text(size=14))
combiplot




    ### CHANGE IN SST OVER TIME ###
## FOR WHOLE YEAR ##
#import sst data
sst<-read_excel("data/sst.xlsx")

sst$year<-as.numeric(substring(sst$date,1,4)) #create year column
sst$day<-yday(sst$date) #create day of year column

#make linear model - with quadratic, cubic, year as random factor
sstmod<-lmer(mean_sst~year+day+I(day^2)+I(day^3)+(1|year),data=sst)

#visualise cubic model
dummyday<-seq(1,365,1) #create dummy day
#create predicted values, add cubic, correct for avg year and intercept
sstpred<-dummyday*-5.798e-02+(dummyday^2)*(6.063e-04)+(dummyday^3)*-1.262e-06 - 8.453e+00
plot(sstpred~dummyday,ylab=expression(paste(
  "Predicted SST (",degree,"C)")),xlab="Day of year") #plot sst by dummyday

summary(sstmod) #view model output

#plot sst over time with year trend
sstplot<-ggplot(aes(x=date,y=mean_sst),data=sst)+geom_line(col="grey")+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Year",
       y=expression(paste(
         "SST (",degree,"C)")))+theme(text = element_text(size=14))
sstplot


#FOR CLIMWIN WINDOW 
#create duplicate sst, sst2
sst2<-sst

sst2$year<-as.numeric(substring(sst2$date,1,4)) #create year column
sst2$day<-yday(sst2$date) #create day of year column

#create subset of critical window in each year
sst2<-subset(sst2,day>=96 & day<=126)

#make linear model
sstmod2<-lmer(mean_sst~year+day+I(day^2)+(1|year),data=sst2)
summary(sstmod2) #view output

#visualise cubic model
dummyday2<-seq(102,132,1) #create dummy day
#create predicted values, add quadratic, correct for intercept
sstpred2<-(dummyday2*-1.602e-01)+(dummyday2^2*8.377e-04)-7.436e+00
#check predictions
plot(sstpred2~dummyday2,ylab=expression(paste(
  "Predicted SST (",degree,"C)")),xlab="Day of year")

summary(sstmod2) #view model output

#plot subsetted sst over time with year trend
sstplot2<-ggplot(aes(x=date,y=mean_sst),data=sst2)+geom_line(col="grey")+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Year",
       y=expression(paste(
         "SST (",degree,"C)")))+theme(text = element_text(size=14))
sstplot2

#plot whole year and climwin period SST on single figure
plot_grid(sstplot,sstplot2,labels=c(
  "A","B"),ncol=1)

    ### EVIDENCE OF TROPHIC MISMATCH ###
## PROPORTION OF SANDEELS IN DIET OVER TIME ##

#create linear model
dietmod<-lm(sandeel~year,data=diet)
par(mfrow=c(2,2))
plot(dietmod) #check assumptions
par(mfrow=c(1,1))
summary(dietmod) #view output

#visualise with ggplot
dietplot<-ggplot(aes(x=year,y=sandeel),data=diet)+geom_point()+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Year",y="Proportion of sandeels 
  in chick diet (%)")+theme(text = element_text(size=14))
dietplot

## BREEDING SUCCESS OVER TIME ##

#create linear model
prodyearmod<-glm(productivity~year,data=productivity,family=gaussian)
par(mfrow=c(2,2))
plot(prodyearmod)#check assumptions
par(mfrow=c(1,1))
#assumptions not met - repeat analysis without outliers

#visualise outputs of original mod
summary(prodyearmod)

#visualise original mod with ggplot
prodyearplot<-ggplot(aes(x=year,y=productivity),data=productivity)+geom_point()+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Year",y="Breeding success")+theme(text = element_text(size=14))
prodyearplot

## REPEAT WITHOUT OUTLIERS ##
#create new dataframe without outliers
prod2<-productivity[-c(9,16,18),] #exclude 2004, 2011, 2013

#repeat analysis
prodyearmod<-lm(productivity~year,data=prod2)
par(mfrow=c(2,2))
plot(prodyearmod) #check assumptions
par(mfrow=c(1,1))
summary(prodyearmod) #look at output

#visualise without outliers
prodyearplot2<-ggplot(aes(x=year,y=productivity),data=prod2)+geom_point()+
  geom_smooth(method=lm,level=0.95)+theme_classic()+
  labs(x="Year",y="Breeding Success")+theme(text = element_text(size=14))
prodyearplot2


## POPULATION OVER TIME ##

#create linear model
popmod<-lm(mean~year,data=pop)
par(mfrow=c(2,2))
plot(popmod) #check assumptions
par(mfrow=c(1,1))
summary(popmod) #view output

#based on distribution of residuals, fit non-linear 
population$year2<-population$year^2
popmodnonlin<-lm(mean~year+year2,data=pop)
par(mfrow=c(2,2))
plot(popmodnonlin) #check assumptions
par(mfrow=c(1,1))
summary(popmodnonlin) #view output

#visualise change in population over time
popplot<-ggplot(aes(x=year,y=mean),data=pop)+geom_point()+
  geom_smooth(method="lm",formula=y~x+I(x^2),level=0.95)+theme_classic()+
  labs(x="Year",y="Mean population 
count of adult guillemots")+theme(text = element_text(size=14))
popplot

## COMBINE ALL TROPHIC MISMATCH FIGURES IN ONE PANEL ##
bottom_row<-plot_grid(dietplot,prodplot,ncol=2,labels="AUTO")
toprow<-plot_grid(popplot,labels=c("C"))
plot_grid(bottom_row,toprow,nrow=2)



    ### FISHERIES DATA ANALYSIS ###

## DATA PREPARATION ##
#calculate median for each fisheries data year
avgfish<-ottertrawl%>% group_by(year) %>%
  summarise_at(vars("densityabundance",
                    "densitybiomass","fishnumber",
                    "fishlength"),median)

#create meanfish variable. exclude years missing from chick diet data with one year lag
avgfish2<-avgfish[-c(1:8),]
avgfish2<-avgfish2[-c(5,7),]

#remove 2021 from diet data (no 2021 fish data)
diet2<-diet[-c(13),]

## DATA ANALYSIS ##
#create model - is % sandeel in chick diet explained by sandeel abundance?
sanddietmod<-lm(diet2$sandeel~meanfish2$densityabundance)
par(mfrow=c(2,2))
plot(sanddietmod) #check assumptions
par(mfrow=c(1,1))
summary(sanddietmod) #view output

#create dataframe with temporally matched meanfish and diet data 
sanddiet<-data.frame(diet2$sandeel,avgfish2$densityabundance)
names(sanddiet)[names(sanddiet) == "diet2.sandeel"] <- "sandeel"
names(sanddiet)[names(sanddiet) == "avgfish2.densityabundance"] <- "densityabundance"

#plot proportion of sandeels in diet against their abundance
sanddietplot<-ggplot(aes(x=densityabundance,y=sandeel),data=sanddiet)+geom_point()+
  geom_smooth(method="lm",level=0.95)+theme_classic()+labs(
    x=expression(Density~of~sandeels~(sandeels~per~km^2)),y=
    "Proportion of sandeels in 
  guillemot chick diet")+theme(text = element_text(size=14))
sanddietplot




    ### POST HOC ANALYSES: POPULATION~DIET ###

#match years with 1-year lag 
pop2<-subset(pop,year>2007)
pop2<-pop2[-c(5,7),]
diet2<-diet[-c(11,13),]

#create linear model
popsandmod<-lm(pop2$mean~diet2$sandeel+diet2$year)
par(mfrow=c(2,2))
plot(popsandmod) #check assumptions
par(mfrow=c(1,1))
summary(popsandmod) #visualise output

#make new df of matched pop and diet data
popsand<-data.frame(pop2$mean,diet2$sandeel,diet2$year)
names(popsand)[names(popsand) == "pop2.mean"] <- "pop"
names(popsand)[names(popsand) == "diet2.sandeel"] <- "sandeel"
names(popsand)[names(popsand) == "diet2.year"] <- "year"

#visualise how diet affects population size
popsandplot<-ggplot(aes(x=sandeel,y=pop),data=popsand)+geom_point()+
  geom_smooth(method="lm",level=0.95)+theme_classic()+geom_point()+
  labs(x="Mean annual proportion of 
sandeels in chick diet (%)",y="Mean annual population size")+theme(text = element_text(size=14))
popsandplot



    ### POST HOC ANALYSES: EFFECT OF SANDEEL ABUNDANCE ON FED IN FOLLOWING YEAR ###

#make both cover FED time period 1999-2019
avgfish3<-avgfish[-c(22),]
fed2<-fed[-c(1:3,25),]

#create dataframe of variables of itnerest
fedfish<-data.frame(avgfish3$densityabundance,fed2$day_of_year,fed2$year)
names(fedfish)[names(fedfish) == "avgfish3.densityabundance"] <- "denabun"
names(fedfish)[names(fedfish) == "fed2.day_of_year"] <- "day"
names(fedfish)[names(fedfish) == "fed2.year"]<-"year"

#create model of effect of avgfish on fed in following year
fedfishmod<-lm(day~denabun,data=fedfish)
par(mfrow=c(2,2))
plot(fedfishmod) #check assumptions
par(mfrow=c(1,1))
summary(fedfishmod) #view output

#plot fed against sandeel abundance
fedfishplot<-ggplot(aes(x=denabun,y=day),data=fedfish)+geom_point()+
  geom_smooth(method="lm",level=0.95)+theme_classic()+labs(
    x=expression(Sandeel~density~per~km^2),
    y="Day of year of first egg date")+theme(text = element_text(size=14))
fedfishplot

    ### END CODE ###