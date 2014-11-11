#####################################################################################################
#
#   Ameland / Imputation as preparation / Jan-Philipp Kolb  / 09.12.09
#       - Version 1 :   - some variables in the EU-SILC data set do have missing values,
#                         these are imputed with that file
#
#####################################################################################################

library(mice)

# setwd("/eingabe/Ameland/universe/Phase1/2005imputed")

   # Current education activity 
load("PE010.cs.RData")

   # Country
load("PB020.cs.RData")

   # PB050: Personal base weight
load("PB040.cs.RData")

   # PB140: Year of birth
load("PB140.cs.RData")

   # PB150: Sex
load("PB150.cs.RData")

   # PB190: Marital status
load("PB190.cs.RData")

   # Highest ISCED level attained - PE040
load("PE040.cs.RData")

   # PH010: General health
load("PH010.cs.RData")
load("PH030.cs.RData")
load("PH040.cs.RData")
load("PH050.cs.RData")
load("PH060.cs.RData")
load("PH070.cs.RData")

load("PL015.cs.RData")
load("PL015_F.cs.RData")

load("PL020.cs.RData")

load("PL025.cs.RData")

load("PL030.cs.RData")

load("PL035.cs.RData")

load("PL040.cs.RData")

   # income variables
load("PY010G.cs.RData")
load("PY010N.cs.RData")
load("PY010N_F.cs.RData")


   # Do you have a car
load("HS110.cs.RData")

load("HB020.cs.RData")
load("HH010.cs.RData")
load("HH020.cs.RData")
load("HH030.cs.RData")
load("HH040.cs.RData")
load("HH050.cs.RData")
load("HH060.cs.RData")
load("HH080.cs.RData")
load("HH090.cs.RData")

load("HX060.cs.RData")
load("HX070.cs.RData")
load("HX080.cs.RData")


names.hh <- dir()
Hs.vars <- grep("HS",names.hh)
Hs.vars <- names.hh[Hs.vars]
Hs.vars.tmp <- grep(".cs",Hs.vars)
Hs.vars <- Hs.vars[Hs.vars.tmp]  
Hs.vars.tmp <- grep("_F",Hs.vars)
Hs.vars <- Hs.vars[-(Hs.vars.tmp)]
Hs.vars.i <- gsub(".cs.RData","",Hs.vars)


eval(parse(text=paste("load('",Hs.vars,"')",sep="",col=";")))




   ###########################
   # getting data sets sorted
   ###########################


pr <- match(PB030,RB030)

   ###########################
   # functions
   ###########################


xvar <- c("PE010")
vars <- c("PB020","PB140")

imput.EUS <- function(xvar, vars){
   eval(parse(text=paste("dat <- data.frame(",paste(xvar,sep="",collapse=","),",",paste(vars,sep="",collapse=","),")")))
   dat.a <- mice(dat)
   dat.b <- complete(dat.a)

      # to save the imputed variables
   eval(parse(text=paste(xvar,"<-  dat.b$",xvar,sep="")))
   eval(parse(text=paste("save(",xvar,",file='/eingabe/Ameland/universe/Phase1/2005imputed/",xvar,".cs.RData')",sep="",collapse=";")))
}

   ###########################
   # application of imputed income variables - HY010
   ###########################

INC <- HY010
INC[is.na(HY010)] <- HY010_I[is.na(HY010)]
   # no advancement


   ###########################
   # Easy Imputation with sex - PB150
   ###########################

PB150[is.na(PB150)] <- sample(1:2,sum(is.na(PB150)),replace=T)

save(PB150, file="/eingabe/Ameland/universe/Phase1/2005imputed/PB150.cs.RData")

   ###########################
   # Easy Imputation with year of birth - PB140
   ###########################

PB140[is.na(PB140)] <- sample(1964:1954,sum(is.na(PB140)),replace=T)

save(PB140, file="/eingabe/Ameland/universe/Phase1/2005imputed/PB140.cs.RData")

   ###########################
   # Person has ever worked - PL015
   ###########################

   # PL030 - Self-defined current economic status

table(is.na(PL015))

table(is.na(PL015[!(PL030==1 & PL030==2)]))

PL015[PL015_F==-3]<-3
PL015[PL015_F==-2]<-4

imput.EUS(xvar="PL015",vars=c("PL030","PB140","PE010","PB150"))

imput.EUS(xvar=c("PL020","PL025","PL030","PL035","PL040"),vars=c("PL015","PB140","PE010","PB150"))


   ###########################
   # Highest ISCED level attained - PE040
   ###########################

imput.EUS(xvar="PE040",vars=c("PB020","PB140","PE010"))

imput.EUS(xvar="PE020",vars=c("PB020","PB140","PE040"))

imput.EUS(xvar="PE030",vars=c("PB020","PB140","PE040","PE020"))

   ###########################
   # Health
   ###########################

imput.EUS(xvar=c("PH010","PH020"),vars=c("PB020","PB140","PE020"))

imput.EUS(xvar=c("PH030","PH040"),vars=c("PH010","PH020","PB140","PE020"))

   ###########################
   # household configuration
   ###########################

imput.EUS(xvar=c(paste("'",Hs.vars.i,"'",sep="",collapse=",")),vars=c("HB020","HH020","HX060","HX080"))

imput.EUS(xvar=c("HS010", "HS020", "HS030"),vars=c("HB020","HH020","HX060","HX080"))

imput.EUS(xvar=c("HS010", "HS020", "HS030"),vars=c("HB020","HH020","HX060","HX080"))

imput.EUS(xvar=c("HS040","HS050","HS060","HS070","HS080"),vars=c("HS010", "HS020", "HS030"))

imput.EUS(xvar=c("HS090","HS100","HS110","HS120"),vars=c("HS040","HS050","HS060","HS070","HS080"))

imput.EUS(xvar=c("HS130", "HS140", "HS150", "HS160") ,vars=c("HS040","HS050","HS060","HS070","HS080","HS010", "HS020", "HS030"))

imput.EUS(xvar=c("HS170","HS180","HS190"),vars=c("HS070","HS080","HS010", "HS020", "HS030"))

   # HH variables

imput.EUS(xvar=c("HH010","HH020","HH030"),vars=c("HS070","HS080","HS010", "HS020", "HS030"))

imput.EUS(xvar=c("HH040","HH050","HH060"),vars=c("HS070","HS080","HS010", "HS020", "HS030"))

imput.EUS(xvar=c("HH040","HH050","HH060"),vars=c("HH010","HH020","HH030", "HS020", "HS030"))

imput.EUS(xvar=c("HH080","HH090"),vars=c("HH010","HH020","HH030", "HH030", "HH040"))

#####################################################################################################
# Alternative method is to sample from the complete cases

gew=PB040/sum(PB040)
gew <- gew[!is.na(PH010)]
weight="gew"

   # region.p has the length of the person dataset and not the length of the register
region.p <- region[pr]

region.s <- region.p[is.na(PH010)]
sovar <- "region.s"
region.o <- region.p[!is.na(PH010)]
sxvar <- "region.o"

PH010.o <- PH010[!is.na(PH010)]
yvar <- "PH010.o"

ACL <- cut(PB140,breaks=seq(min(PB140,na.rm=T)-1,max(PB140,na.rm=T),age.dist),labels=F)
ACL.s <- ACL[!is.na(PH010)]


ACL.o <- PB140[is.na(PH010)]
PB150.o <- PB150[is.na(PH010)]
PE010.o <- PE010[is.na(PH010)]
PB190.o <- PB190[is.na(PH010)]

ovars <- c("ACL.o","PB150.o","PE010.o","PB190.o")

age.dist <- 3
PB150.s <- PB150[!is.na(PH010)]
PE010.s <- PE010[!is.na(PH010)]
PB190.s <- PB190[!is.na(PH010)]

xvars <- c("ACL.s","PB150.s","PE010.s","PB190.s")

PH010 <- var.mod(sovar, sxvar,yvar,ovars,xvars,weight)

##################################################################
# Imputation mit Haushaltseinkommen          9.03.2010

# es werden für diese Imputation nur Variablen verwendet, die auch schon im Datensatz sind

setwd("/eingabe/Ameland/universe/Phase1/2005")

load("HY010.cs.RData")
load("HY020.cs.RData")

load("HB020.cs.RData")
load("/eingabe/Ameland/universe/Phase1/2005imputed/HH020.cs.RData")
load("HH070.cs.RData")
load("HX080.cs.RData")

load("HY030G.cs.Rdata")
load("HY040G.cs.RData")
load("HY050G.cs.RData")
load("HY080G.cs.RData")
load("HY100G.cs.RData")
load("HY110G.cs.RData")
load("HY120G.cs.RData")
load("HY130G.cs.RData")
load("HY140G.cs.RData")
load("HY145N.cs.RData")



   # alle Variablen heraus finden, die etwas mit dem Einkommen zu tun haben

Namen <- dir()

ind <- grep("HY",Namen)

Namen[ind]

   # Variablen mit einem F hinten dran sind Flags, daran sieht man, ob die Variable imputiert wurde
load("HY010_F.cs.RData")
load("HY010_I.cs.RData")

setwd("/eingabe/Ameland/universe/Phase1/2005imputed")

table(is.na(HY010_I))
table(is.na(HY010))

table(is.na(HY010[!is.na(HY010_I)]))

   # das bringt leider nichts, es sind genauso viele missings da, wie vorher
HY010[is.na(HY010)] <- HY010_I[is.na(HY010)]

na.test <- function(a,b){
   table(is.na(a[!is.na(b)]))
}

na.test(HY010,HY020)
table(is.na(HY010))



 dat.a <- data.frame(HY010,HY020,HH070,HY030G,HY040G,HY050G,HY080G,HY100G,HY110G,HY120G,HY130G,HY140G,HY145N)

   # in folgendem Datensatz ist das Land mit drin, das ist aber ein Faktor
   # HH020 ist der Tenure status, also eine kategoriale Variable

# dat.a <- data.frame(HY010,HY020,HB020,HH020,HH070,HX080,HY030G,HY040G)

# eigentlich sollten nur normalverteilte Daten vorhanden sein

str(dat.a)

   dat <- mice(dat.a,imputationMethod="norm")
   dat.b <- complete(dat)

Namen <- colnames(dat.b)

for (i in Namen){
   eval(parse(text=paste(i,"<- dat.b$",i,sep=""))) 
   cat(i, "\n")
   eval(parse(text=paste("save(",i,",file='/eingabe/Ameland/universe/Phase1/2005imputed/",i,".cs.RData')",sep="")))
}

   # viel Mist abgespeichert, muss wieder gelöscht werden

Namen <- dir()

ind <- grep("cs.RData.cs.RData",Namen)

del1 <- Namen[ind]

for (i in del){
   eval(parse(text=paste("file.remove('",i,"')",sep="")))
}


##################################################################
# Imputation der synthetischen Variablen

imp.vars <- c("PE010","PE040","RB200","RB210","ACL")

eval(parse(text=paste("load('/eingabe/Ameland/universe/results_v1/AML.",imp.vars,".RData')",sep="",collapse=";")))


   # die P Variablen dürfen in manchen Fällen keine Einträge haben

load('/eingabe/Ameland/universe/results_v1/AML.FILE.RData')

eval(parse(text=paste("dat <- data.frame(",paste(imp.vars,collapse=','),")",sep="")))


##################################################################
# Imputation des synthetischen Haushaltseinkommens

load('/eingabe/Ameland/universe/results_v1/AML.AGE.RData')

   # Medianalter

AGE2HID <- tapply(AGE,HID,median,na.rm=T)


   # Anzahl d erwerbstätigen im HH

   # Variable at work

AW <- rep(0,length(RB210))
AW[RB210=="1"] <- 1

AW2HID <- tapply(AW,HID,sum)


   # Retiremnent
RET <- rep(0,length(RB210))
RET[RB210=="3"] <- 1
RET2HID <- tapply(RET,HID,sum)


   # Haushaltsgröße

dat.a <- data.frame(HHG.hh,INC.HY010.d,AW2HID,RET2HID,AGE2HID)

dat <- mice(dat.a,imputationMethod="norm")
dat.b <- complete(dat)

INC.HY010.i <- dat.b$INC.HY010.d

save(INC.HY010.i,file="/eingabe/Ameland/DSGpresentation/DSG.INC.HY010.i.RData")



