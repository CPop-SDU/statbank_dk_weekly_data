# Forecast death counts by epi-year
#
# Marie-Pier Bergeron Boucher
#
# 2020-09-30

# Libraries
require(compositions)
require(demography)
require(tidyverse)
require(ungroup)

#Source of functions
source("./code/MPs_epi_year_death_count_forecast/Functions_projections.R")

#####################################################################################
## -------------------------------- Data ------------------------------------------##
#####################################################################################

### ----- Mortality ----- ###

#Weekly deaths by 5-year age-groups

load("./out/2020-09-30-dk_weekly_deaths_quarterly_population.RData")
dta_mortality<-dk_weekly_deaths_quarterly_population
head(dta_mortality)

#Deaths as matrix
Dx_EpiYear<- tapply(dta_mortality$deaths, INDEX = list(dta_mortality$start_of_age_group,
                                                           dta_mortality$epi_year,
                                                       dta_mortality$sex),
                        FUN=sum)
Dx_EpiYear_Men<-Dx_EpiYear[,,"Men"]
Dx_EpiYear_Women<-Dx_EpiYear[,,"Women"]

# Ungroup up to 105+
Age = unique(dta_mortality$start_of_age_group)
Age_1y<-0:105
nlast = 6

Dx_EpiYear_Women_ungroup <- matrix(0, nrow=length(Age_1y), ncol=dim(Dx_EpiYear_Women)[2])
for (i in 1:dim(Dx_EpiYear_Women_ungroup )[2]){
  Dx_EpiYear_Women_ungroup [,i] <- pclm(x=Age, y=Dx_EpiYear_Women[,i], nlast=nlast)$fitted
}
colnames(Dx_EpiYear_Women_ungroup)<-colnames(Dx_EpiYear_Women)
rownames(Dx_EpiYear_Women_ungroup)<-Age_1y

Dx_EpiYear_Men_ungroup <- matrix(0, nrow=length(Age_1y), ncol=dim(Dx_EpiYear_Men)[2])
for (i in 1:dim(Dx_EpiYear_Men_ungroup )[2]){
  Dx_EpiYear_Men_ungroup [,i] <- pclm(x=Age, y=Dx_EpiYear_Men[,i], nlast=nlast)$fitted
}
colnames(Dx_EpiYear_Men_ungroup)<-colnames(Dx_EpiYear_Men)
rownames(Dx_EpiYear_Men_ungroup)<-Age_1y

### ----- Population at the first day of the 4th quarter  ----- ###

PopDk<-read.csv("code/MPs_epi_year_death_count_forecast/PopQ3_Denmark.csv", header=T, sep=",")
head(PopDk)

Px_Women<- as.matrix(PopDk[PopDk$Sex=="Women",3:ncol(PopDk)])
rownames(Px_Women)<-Age_1y
colnames(Px_Women)<-2008:2019

Px_Men<- as.matrix(PopDk[PopDk$Sex=="Men",3:ncol(PopDk)])
rownames(Px_Men)<-Age_1y
colnames(Px_Men)<-2008:2019

### ----- Calculate exposure and death rates ----- ###

Ex_EpiYear_Men<- (Px_Men[,1:(ncol(Px_Men)-1)]+Px_Men[,2:ncol(Px_Men)])/2
Ex_EpiYear_Women<- (Px_Women[,1:(ncol(Px_Women)-1)]+Px_Women[,2:ncol(Px_Women)])/2

mx_EpiYear_Men<-Dx_EpiYear_Men_ungroup[, 2:(ncol(Dx_EpiYear_Men_ungroup)-2)]/Ex_EpiYear_Men
mx_EpiYear_Women<-Dx_EpiYear_Women_ungroup[, 2:(ncol(Dx_EpiYear_Women_ungroup)-2)]/Ex_EpiYear_Women

matplot(log(mx_EpiYear_Men), type="l",  col=rainbow(ncol(mx_EpiYear_Men)))
matplot(log(mx_EpiYear_Women), type="l",  col=rainbow(ncol(mx_EpiYear_Men)))

### ----- Fertility ----- ###

### Load data
fx<-read.csv("code/MPs_epi_year_death_count_forecast/ASFR.csv", skip=2, header=T, sep=";")
head(fx)
tail(fx)
dim(fx)

#####################################################################################
## ---------------------------Life table and parameters----------------------------##
#####################################################################################

##Number of year ahead to be forecast (at least 2)
h=2

## Parameters
Y1<-2008
Y2<-2012:2018 #last year used for forecast, year epi-year" started
Yfcst<-Y2+1 #year forecast

## ------- Calculate life table 

#Defined r and nx
r<-1
nx<-1

# Females
LT_f<-LifeT_mx(mx_EpiYear_Women, r, nx)

LT_f$ex[1,]

#Males
LT_m<-LifeT_mx(mx_EpiYear_Men, r, nx)
LT_m$ex[1,]

## -------Fertility matrix

Afx<-15:49

#create matrix 
fx_mat<-as.matrix(fx[2:36,2:ncol(fx)]/1000)
rownames(fx_mat)<-Afx
colnames(fx_mat)<-1973:2019

#replace 0
fx_mat[fx_mat==0]<-min(fx_mat[fx_mat>0])/2

# Assume that fertility rates stay constant to the last calendar year forecast

#####################################################################################
## ---------------------------Projections----------------------------##
#####################################################################################

Pop_forecast_Women<- Pop_forecast_Men<-
  Dx_forecast_Women<-Dx_forecast_Men<- matrix(0, length(Age_1y), length(Y2))

#Rename colunms of life table dx 
colnames(LT_f$dx)<-colnames(LT_m$dx)<-2008:2018

#Loop for different starting year

for(i in 1:length(Y2)){
  
  #forecast females mortality independently 
  
  DKfcst_f<-CoDA.fun(dx=LT_f$dx[,as.character(Y1:Y2[i])], h=h, nsim=100)
  
  
  
  #forecast males mortality using females as reference
  
  DKfcst_m<-CoDAC.fun(dx=LT_m$dx[,as.character(Y1:Y2[i])], dx.ref=LT_f$dx[,as.character(Y1:Y2[i])],
                      h=h, nsim1=100, nsim2=100)
  
  
  
  # Population count on Oct.1 of the last year used in the fitting period
  
  PopF<-Px_Women[,as.character(Y2[i])]
  PopM<-Px_Men[,as.character(Y2[i])]
  
  
  ### Projections
  
  #names LT and fx 
  colnames(DKfcst_f$LT$Lx)<-colnames(DKfcst_m$LT$Lx)<-
    colnames(DKfcst_f$LT$lx)<-colnames(DKfcst_m$LT$lx)<-
    colnames(DKfcst_f$LT$Tx)<-colnames(DKfcst_m$LT$Tx)<-
    colnames(DKfcst_f$LT$mx)<-colnames(DKfcst_m$LT$mx)<-
    
    colnames(DKfcst_f$LT.lw$Lx)<-colnames(DKfcst_m$LT.lw$Lx)<-
    colnames(DKfcst_f$LT.lw$lx)<-colnames(DKfcst_m$LT.lw$lx)<-
    colnames(DKfcst_f$LT.lw$Tx)<-colnames(DKfcst_m$LT.lw$Tx)<-
    colnames(DKfcst_f$LT.lw$mx)<-colnames(DKfcst_m$LT.lw$mx)<-
    
    colnames(DKfcst_f$LT.up$Lx)<-colnames(DKfcst_m$LT.up$Lx)<-
    colnames(DKfcst_f$LT.up$lx)<-colnames(DKfcst_m$LT.up$lx)<-
    colnames(DKfcst_f$LT.up$Tx)<-colnames(DKfcst_m$LT.up$Tx)<-
    colnames(DKfcst_f$LT.up$mx)<-colnames(DKfcst_m$LT.up$mx)<- Y2[i]+1:h
  
  fx.fcst<-cbind(fx_mat[,as.character(Y2[i])], fx_mat[,as.character(Y2[i])])
  colnames(fx.fcst)<-Y2[i]+1:h
  rownames(fx.fcst)<-Afx
  
  
  #CCP main forecast
  pop.fcst<-CC_fun(LxF=DKfcst_f$LT$Lx,
                   lxF=DKfcst_f$LT$lx,
                   TxF=DKfcst_f$LT$Tx,
                   
                   LxM=DKfcst_m$LT$Lx,
                   lxM=DKfcst_m$LT$lx,
                   TxM=DKfcst_m$LT$Tx,
                   
                   Fx=fx.fcst, 
                   
                   popF=PopF,
                   popM=PopM)
  
  
  # Population on Sept. 1
  
  PopF_fcst<-cbind(PopF, pop.fcst$projF)
  PopM_fcst<-cbind(PopM, pop.fcst$projM)
  
  
  colnames(PopF_fcst)<-colnames(PopM_fcst)<- Y2[i]+0:h
  
  
  
  # Exposure
  
  ExF_fcst<-(PopF_fcst[,1:(ncol(PopF_fcst)-1)]+PopF_fcst[,2:(ncol(PopF_fcst))])/2
  ExM_fcst<-(PopM_fcst[,1:(ncol(PopM_fcst)-1)]+PopM_fcst[,2:(ncol(PopM_fcst))])/2
  
  
  # Death counts
  
  dxF_fcst<-ExF_fcst*DKfcst_f$LT$mx
  dxM_fcst<-ExM_fcst*DKfcst_m$LT$mx
  
  
  #Save results
  
  Pop_forecast_Women[,i]<- PopF_fcst[,2]
  Pop_forecast_Men[,i]<-PopM_fcst[,2]
  
  Dx_forecast_Women[,i]<-dxF_fcst[,1]
  Dx_forecast_Men[,i]<-dxM_fcst[,1]
  
  
} #end of loop

#Rename columns with proper epi-year
epiyear_fcst<-Y2
for(i in 1:length(Y2)){
  epiyear_fcst[i]<-c(paste0(Y2[i]+1,"/",Y2[i]+2))
}

colnames(Pop_forecast_Women)<-colnames(Pop_forecast_Men)<-Y2+1
colnames(Dx_forecast_Women)<-colnames(Dx_forecast_Men)<-epiyear_fcst

#####################################################################################
## --------------------------- Data set -------------------------------------------##
#####################################################################################

dta<-expand.grid(Age=0:105,
                 EpiYear=epiyear_fcst)

dta$Dx_forecast_Women<-c(Dx_forecast_Women)
dta$Dx_forecast_Men<-c(Dx_forecast_Men)

age.cat<-c(seq(0,100,5),Inf)
dta$Age_5y<-cut(dta$Age, age.cat, right = F)

dta5<- dta[,-1] %>%
  group_by(EpiYear, Age_5y) %>%
  summarise_all(sum)
dta5$Dx_Observed_Women<-c(Dx_EpiYear_Women[,7:13])
dta5$Dx_Observed_Men<-c(Dx_EpiYear_Men[,7:13])

dta5$Error_Women<-dta5$Dx_forecast_Women-dta5$Dx_Observed_Women
dta5$Error_Men<-dta5$Dx_forecast_Men-dta5$Dx_Observed_Men

# harmonize with xSTMF data series
dk_mp_epi_year_death_forecast <-
  dta %>%
  mutate(
    age_group = cut(Age, breaks = c(0, 15, 65, 75, 85, Inf), right = FALSE)
  ) %>%
  select(epi_year = EpiYear, age_group, Dx_forecast_Men, Dx_forecast_Women) %>%
  pivot_longer(
    cols = c(Dx_forecast_Men, Dx_forecast_Women),
    names_to = 'sex',
    values_to = 'mp_epi_year_death_frcst'
  ) %>%
  mutate(
    sex = ifelse(grepl('_Men$', sex), 'Male', 'Female') %>%
      factor(levels = c('Female', 'Male'))
  ) %>%
  group_by(epi_year, sex, age_group) %>%
  summarise(
    deaths_per_epi_year_fcast_mp = sum(mp_epi_year_death_frcst)
  ) %>%
  ungroup() %>%
  mutate(
    country_code = 'DNK',
    epi_year = as.character(epi_year)
  )

save(
  dk_mp_epi_year_death_forecast,
  file = paste0('out/', Sys.Date(), '-dk_mp_epi_year_death_forecast.RData')
)

ggplot(dta5[dta5$EpiYear!="2019/2020",])+
  geom_point(aes(x=Age_5y, y=Dx_Observed_Women))+
  geom_point(aes(x=Age_5y, y=Dx_forecast_Women), col=2, pch=1)+
  facet_wrap(~EpiYear, ncol=2)+
  ylab("Deaths")+
  xlab("Age-group")+
  ggtitle("Death counts observed (black) and forecast (red), Epi-Year, Females")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))

ggplot(dta5[dta5$EpiYear!="2019/2020",])+
  geom_point(aes(x=Age_5y, y=Error_Women))+
  geom_hline(yintercept = 0)+
  ylab("Errors")+
  xlab("Age-group")+
  ggtitle("Errors on death counts observed, Epi-Year, Females")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))

ggplot(dta5[dta5$EpiYear!="2019/2020",])+
  geom_point(aes(x=Age_5y, y=Dx_Observed_Men))+
  geom_point(aes(x=Age_5y, y=Dx_forecast_Men), col=2, pch=1)+
  facet_wrap(~EpiYear, ncol=2)+
  ylab("Deaths")+
  xlab("Age-group")+
  ggtitle("Death counts observed (black) and forecast (red), Epi-Year, Males")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))

ggplot(dta5[dta5$EpiYear!="2019/2020",])+
  geom_point(aes(x=Age_5y, y=Error_Men))+
  geom_hline(yintercept = 0)+
  ylab("Errors")+
  xlab("Age-group")+
  ggtitle("Errors on death counts observed, Epi-Year, Males")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))
