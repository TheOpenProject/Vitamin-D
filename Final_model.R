############################################################################################################
# Programme:  Vitamin D PBPK Final model
# Purpose:    The purpose is to reproduce the simulation results in the following paper (Figure S6 and 
#             Figure S7) to verify the results.
#               
#               Huang, Zhonghui & You, Tao* (2021) Personalise Vitamin D3 Using PBPK Modelling. CPT:PSP.
#               * Correspondence: tao.you@letsgobeyond.co.uk
#               Beyond Consulting Ltd. 14 Tytherington Park Road, Macclesfield, Cheshire, UK  SK10 2EL
# 
# R packages: deSolve, FME, ggplot2
# Input Data: Training set - "25(OH)Dpc10_100.csv"
#             Figure S6. Visual predictive check of 25(OH)D after repeat dose (training set).
#               Panel 1-26: 400IU (10μg) QD
#               Panel 27-43: 4000IU (100μg) QD
#               Number of each panel is listed as “Train No.” in Table S3.
#             Test set - "25(OH)Dtestset.csv"
#             Figure S7. Visual predictive check of 25(OH)D after repeat dose (test set).
#               All doses except 10μg and 100μg)
#               Number of each panel is listed as “Test No.” in Table S3.
# Output:     a. VPC plot: training set (daily repeat dose)
#             b. VPC plot: test set (daily repeat dose)
############################################################################################################
# 0) Load library and data
# 1) Vitamin D PBPK function
# 1.A) Specify differential equations
# 1.B) Specify initial conditions
# 
# 2) Simulate VD PBPK model and output VPC plot
# 2.A) Simulate VD PBPK modela and output VPC plot: training set (daily repeat dose - Figure S6)
# 2.B) Simulate VD PBPK modela and output VPC plot: test set (daily repeat dose - Figure S7)

############################################################################################################
# 0) Load library and data
rm(list=ls())

library(deSolve)  # tested in v1.28
library(ggplot2)  # tested in v3.3.3
library(FME)      # tested in v1.3.6.1
# setwd("../project/VDmodel")  # You need to set the work directory to the path where csv files are located
getwd()

#############################################################################################################
# 1) Vitamin D PBPK function

VitaminD_PBPK<- function (pars,Duration,DOSE) {
  derivs <- function(time, y, pars) {
    with (as.list(c(pars, y)), {
      #                     Species               Compartment
      A_dep   = y[1]   # Vitamin D3         Depot compartment
      A_ven   = y[2]   # Vitamin D3        Venous compartment
      A_l     = y[3]   # Vitamin D3         Liver compartment
      A_rb    = y[4]   # Vitamin D3              Rest of body
      A_art   = y[5]   # Vitamin D3       Aterial compartment
      A_col   = y[6]   # Vitamin D3                Collection (check mass balance)
      A25_ven = y[7]   #   25(OH)D3        Venous compartment
      A25_l   = y[8]   #   25(OH)D3         Liver compartment
      A25_rb  = y[9]   #   25(OH)D3             Rest of body
      A25_art = y[10]  #   25(OH)D3       Aterial compartment
      A25_col = y[11]  #   25(OH)D3                Collection (check mass balance)
      
      # 1.A) Specify differential equations
      # 1 Depot compartment (GI): dA_GI
      dA_dep = -ka*A_dep
      
      # 2 Venous compartment: dA_ven
      dA_ven = Ql*A_l/Vl/Kpl+Qrb*A_rb/Vrb/Kprb-A_ven/Vven*Qco
      
      # 3 liver compartment: dA_l
      dA_l = A_art/Vart*Ql-Ql*A_l/Vl/Kpl-MCLH*A_l/Vl/Kpl+(CL_0+CL_EMAX*D25BASE**CL_GAMMA/(D25BASE**CL_GAMMA+CL_C50**CL_GAMMA))*D25BASE*3+ka*A_dep
      
      # 4 other compartment: dA_rb
      dA_rb = A_art/Vart*Qrb-Qrb*A_rb/Vrb/Kprb
      
      # 5 Aterial compartment: dA_art
      dA_art = A_ven/Vven*Qco-A_art/Vart*Ql-A_art/Vart*Qrb
      
      # 6 VD Collection (for mass balance): dA_COL
      dA_col = MCLH*A_l/Vl/Kpl
      
      # 7 25(OH)D Venous compartment: dA25_ven
      dA25_ven = Ql*A25_l/Vl/Kp25l+Qrb*A25_rb/Vrb/Kp25rb-A25_ven/Vven*Qco
      
      # 8 25(OH)D liver compartment: dA25_l
      dA25_l = A25_art/Vart*Ql-Ql*A25_l/Vl/Kp25l+MCLH*A_l/Vl/Kpl/3-(CL_0+CL_EMAX*(A25_l/Vl/Kp25l)**CL_GAMMA/((A25_l/Vl/Kp25l)**CL_GAMMA +CL_C50**CL_GAMMA))*(A25_l/Vl/Kp25l)
      
      # 9 25(OH)D other compartment: dA25_rb
      dA25_rb = A25_art/Vart*Qrb-Qrb*A25_rb/Vrb/Kp25rb
      
      # 10 25(OH)D Aterial compartment: dA25_art
      dA25_art = A25_ven/Vven*Qco-A25_art/Vart*Ql-A25_art/Vart*Qrb
      
      # 11 25(OH)D Collection (for mass balance): dA25_COL
      dA25_col = (CL_0+CL_EMAX*(A25_l/Vl/Kp25l)**CL_GAMMA/((A25_l/Vl/Kp25l)**CL_GAMMA+CL_C50**CL_GAMMA))*(A25_l/Vl/Kp25l)
      
      return(list(c(dA_dep,dA_ven,dA_l,dA_rb,dA_art,dA_col,
                    dA25_ven,dA25_l,dA25_rb,dA25_art,dA25_col), # ODEs
                  CP = unname(A25_ven/Vven), # 25(OH)D concentration in venous blood
                  V25Dall = unname(A25_ven+A25_l+A25_rb+A25_art+A25_col-A_col/3) # Mass balance
      ))	                  
    })
  }
  
  # 1.B) Specify initial conditions
  # Convert dose: μg => nmol/L
  DOSE<-with(as.list(pars), DOSE/384.64*1000)
  # Set initial conditions according to mass balance
  A_dep_0<-0
  A_ven_0<-with(as.list(pars), (CL_0+CL_EMAX*D25BASE^CL_GAMMA/(D25BASE^CL_GAMMA+CL_C50^CL_GAMMA))*3*D25BASE/MCLH*Vven)
  A_l_0<-with(as.list(pars),(CL_0+CL_EMAX*D25BASE^CL_GAMMA/(D25BASE^CL_GAMMA+CL_C50^CL_GAMMA))*3*D25BASE/MCLH*Vl*Kpl)
  A_rb_0<-with(as.list(pars),(CL_0+CL_EMAX*D25BASE^CL_GAMMA/(D25BASE^CL_GAMMA+CL_C50^CL_GAMMA))*3*D25BASE/MCLH*Vrb*Kprb)
  A_art_0<-with(as.list(pars), (CL_0+CL_EMAX*D25BASE^CL_GAMMA/(D25BASE^CL_GAMMA+CL_C50^CL_GAMMA))*3*D25BASE/MCLH*Vart)
  A_col_0<-with(as.list(pars), 0)
  A25_ven_0<-with(as.list(pars), D25BASE*Vven)
  A25_l_0<-with(as.list(pars), D25BASE*Vl*Kp25l)
  A25_rb_0<-with(as.list(pars),D25BASE*Vrb*Kp25rb)
  A25_art_0<-with(as.list(pars), D25BASE*Vart)
  A25_col_0<-with(as.list(pars), 0)
  # Specify initial conditions
  y <- c(A_dep=A_dep_0,
         A_ven=A_ven_0,
         A_l=A_l_0,
         A_rb=A_rb_0,
         A_art=A_art_0,
         A_col=A_col_0,
         A25_ven=A25_ven_0,
         A25_l=A25_l_0,
         A25_rb=A25_rb_0,
         A25_art=A25_art_0,
         A25_col=A25_col_0) 
  
  # Specify time points to simulate
  times <- c(seq(0,Duration,24))
  times<-times[order(times)]
  
  # ndoses: number of doses
  TIMElast <- max(times)
  ndoses <- TIMElast/24+1
  # Construct a data frame to store dosing information
  DOSEdata <- data.frame(var    = rep(1, times = ndoses),
                         time   = seq(0,TIMElast,24),
                         value  = rep(DOSE, times = ndoses),
                         method = rep("add", times = ndoses))
  
  # Execute numerical simulation
  out <- data.frame(ode(y = y,
                        times = times,
                        parms = pars,
                        func = derivs,
                        events = list(data = DOSEdata),
                        method=c("daspk"),
                        hmax=1,
                        atol = 1e-4,
                        rtol = 1e-4,
                        nout=1,
                        outnames="CP"))
  return(out)
}

#############################################################################################################
# 2) Simulate VD PBPK model and output VPC plot
# 2.A) Simulate VD PBPK model and output VPC plot: training set (daily repeat dose)

MCMC=readRDS(file="MCMC-Run009b.RDS")

VDdataset<-read.csv(file="25OHDpc10_100.csv",stringsAsFactors = FALSE)
rownames(VDdataset)<-seq(1, nrow(VDdataset),1)
VDdataset$PointID = seq(1:dim(VDdataset)[1])
Obsdata_uni<-subset(VDdataset,Time==0)

################################################################################
# “Train No.” in Table S3
i=43
################################################################################
Flag=Obsdata_uni[i,]$Flag
Duration=Obsdata_uni[i,]$Duration*24
m<-Obsdata_uni[i,]$Conc
DOSE<-Obsdata_uni[i,]$Dose

# Specify parameters
param <-c(Qco=312,Ql=70.824,Qrb=241.176,
          Vven=4.2,Vl=1.80,Vrb=62.60,Vart=1.40,
          Kpl=1,Kprb=0.1,MCLH=0.2,
          Kp25l=1,
          Kp25rb=exp(summary(MCMC)["Kp25rb"]["mean",]),
          CL_C50=exp(summary(MCMC)["CL_C50"]["mean",]),
          CL_GAMMA=exp(summary(MCMC)["CL_GAMMA"]["mean",]),
          CL_0=0,
          CL_EMAX=exp(summary(MCMC)["CL_EMAX"]["mean",]),
          ka=0.2,
          D25BASE=m)

# Sample the model with MCMC results
sR1 <- sensRange(func = VitaminD_PBPK,
                 parms = param,
                 parInput = exp(MCMC$pars),
                 Duration=Duration,
                 DOSE=DOSE,
                 sensvar = c("CP"))

out=summary(sR1)
out<-data.frame( out$x/24,out$Mean,out$q05,out$q95)
out$ID=Obsdata_uni[i,]$ID
colnames(out)<-c("Time","Mean","q05","q95","ID")
result<-out
################################################################################
# Annotate simulated "result" data frame with data
result$Conc<-NA
for (i in 1:nrow(VDdataset)) {
  Conc<-VDdataset[i,]$Conc
  Time<-VDdataset[i,]$Time
  ID<-VDdataset[i,]$ID
  for (i in 1:nrow(result)){
    if (result[i,]$ID==ID & result[i,]$Time==Time ){
      result[i,]$Conc=Conc
    }
  }
}

################################################################################
# Plot simulation with data
ggplot(result)+
  geom_point(aes(x = Time, y = Conc), col="darkblue",size = 0.5)+
  geom_line(aes(x = Time, y = Mean), col="darkblue",size = 0.5)+
  geom_ribbon(aes(x = Time, ymin = q05, ymax = q95),fill="darkblue",alpha=0.3)+
  scale_y_continuous(" Plasma 25(OH)D (nmol/L)",limits = c(1,1000),trans = "log10")+
  scale_x_continuous("Days")+
  guides(fill=guide_legend(title="ID"),col=guide_legend(title="ID"))+
  theme_bw()+
  theme(axis.title.x=element_text(vjust=1,  
                                  size=10,
                                  color = "black"),  # X axis title
        axis.title.y=element_text(size=10,
                                  color = "black"),  # Y axis title
        axis.text.x=element_text(size=8, 
                                 angle = 0,
                                 color = "black",
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=8,
                                 angle = 0,
                                 color = "black")) +
  theme(legend.title = element_text(size=10, color = "black", face="bold"),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())+
  
  theme(strip.background =element_rect(fill="#CCCCCC"))+
  theme(strip.text= element_text(size = 8, colour = "black",face="bold"))


################################################################################
# 2.B) Simulate VD PBPK model and output VPC plot: test set (daily repeat dose)
VDdataset<-read.csv(file="25(OH)Dtestset.csv",stringsAsFactors = FALSE)
colnames(VDdataset)<-c("Flag","ID","Dose","Duration","Time","Conc","sd")
rownames(VDdataset)<-seq(1, nrow(VDdataset),1)
VDdataset$PointID = seq(1:dim(VDdataset)[1])
Obsdata_uni<-subset(VDdataset,Time==0)
################################################################################
# “Test No.” in Table S3.
i=83
################################################################################

Flag=Obsdata_uni[i,]$Flag
Duration=Obsdata_uni[i,]$Duration*24
m<-Obsdata_uni[i,]$Conc
DOSE<-Obsdata_uni[i,]$Dose

# Specify parameters
param <-c(Qco=312,Ql=70.824,Qrb=241.176,
          Vven=4.2,Vl=1.80,Vrb=62.60,Vart=1.40,
          Kpl=1,Kprb=0.1,MCLH=0.2,
          Kp25l=1,
          Kp25rb=exp(summary(MCMC)["Kp25rb"]["mean",]),
          CL_C50=exp(summary(MCMC)["CL_C50"]["mean",]),
          CL_GAMMA=exp(summary(MCMC)["CL_GAMMA"]["mean",]),
          CL_0=0,
          CL_EMAX=exp(summary(MCMC)["CL_EMAX"]["mean",]),
          ka=0.2,
          D25BASE=m)

#Simulated data for individual
sR1 <- sensRange(func = VitaminD_PBPK,
                 parms = param,
                 parInput = exp(MCMC$pars),
                 Duration=Duration,
                 DOSE=DOSE,
                 sensvar = c("CP"))

out=summary(sR1)
out<-data.frame( out$x/24,out$Mean,out$q05,out$q95)
out$ID=Obsdata_uni[i,]$ID
colnames(out)<-c("Time","Mean","q05","q95","ID")
result<-out

################################################################################
# Annotate simulated "result" data frame with data
result$Conc<-NA
for (i in 1:nrow(VDdataset)) {
  Conc<-VDdataset[i,]$Conc
  Time<-VDdataset[i,]$Time
  ID<-VDdataset[i,]$ID
  for (i in 1:nrow(result)){
    if (result[i,]$ID==ID & result[i,]$Time==Time ){
      result[i,]$Conc=Conc
    }
  }
}
# Plot simulation with data
ggplot(result)+
  geom_point(aes(x = Time, y = Conc), col="darkblue",size = 0.5)+
  geom_line(aes(x = Time, y = Mean), col="darkblue",size = 0.5)+
  geom_ribbon(aes(x = Time, ymin = q05, ymax = q95),fill="darkblue",alpha=0.3)+
  scale_y_continuous(" Plasma 25(OH)D (nmol/L)",limits = c(1,1200),trans = "log10")+
  scale_x_continuous("Days")+
  guides(fill=guide_legend(title="ID"),col=guide_legend(title="ID"))+
  theme_bw()+
  theme(axis.title.x=element_text(vjust=1,  
                                  size=10,
                                  color = "black"),  # X axis title
        axis.title.y=element_text(size=10,
                                  color = "black"),  # Y axis title
        axis.text.x=element_text(size=8, 
                                 angle = 0,
                                 color = "black",
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=8,
                                 angle = 0,
                                 color = "black")) +
  theme(legend.title = element_text(size=10, color = "black", face="bold"),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())+
  
  theme(strip.background =element_rect(fill="#CCCCCC"))+
  theme(strip.text= element_text(size = 8, colour = "black",face="bold"))
