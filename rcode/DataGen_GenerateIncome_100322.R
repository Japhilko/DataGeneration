#####################################################################################################
#
# Einkommen generieren - Vorbereitung zur DSG Tagung
#     Projekt: AMELI
#     Autor: Jan-Philipp Kolb
#       - Version 3:   - Einkommen soll generiert werden
#                      - oben nicht zwei verschiedene Sorten von Quantilen
#                      - Erweiterung um Einkommenskomponenten
#
#####################################################################################################

Year <- 2005

   # wenn das Jahr des Surveys rauß gefunden werden soll:

load("HB010.cs.RData")

############################################################
# Einkommen ziehen


   # hier liegt das Einkommen auf dem lokalen Rechner
# load('/home/kolb/Desktop/Ameli/AmeLand/DSGVortrag/data/HY010.cs.RData' )

generate.inc <- function(inc){
 
      # negative Einkommen
   neg.inc <- inc[inc<0]
   quant.neg <- sort(-quantile(-neg.inc, probs = seq(0, 1, 0.25),na.rm=T))
   quant <- quantile(inc, probs = seq(0, 1, 0.005),na.rm=T)

   mpos.inc <- inc[inc>0]
   quant.mpos <- quantile(mpos.inc, probs = seq(0, 1, 0.06),na.rm=T)

   bord <- c(quant.neg,quant.mpos)
   

   INC <- cut(inc,breaks=bord,include.lowest=T,labels=F)
   INC[inc==0]<-0

   min.inc <- tapply(inc,INC,min)
   max.inc <- tapply(inc,INC,max)

   inc.cl <- sort(unique(INC))
   attr(INC,"min.inc")<- min.inc
   attr(INC,"max.inc")<- max.inc
   attr(INC,"inc.cl")<- inc.cl
   attr(INC,"date")<- date()
   attr(INC,"function")<- "generate.inc"
   return(INC)
}


get.dis.Inc3 <- function(inc,inc.basic,dis.basic){

   ff<-function(x){
    lol<-hist(x,nclass=100,plot=F)
    lol<-lm(lol$counts ~ lol$breaks[-1])
    return(function(y){coef(lol)[1] + y*coef(lol)[2]})
   }

   AC.RE1 <- function() {
       ACC <- F
       while(!ACC) { 
           a <- runif(1,min.tinc,max.tinc)
           b<-runif(1,0,max.hist)        
           if (b <= f1(a)) { ACC=T }
       }
       return(a)
   }


   EK <- rep(NA,length(inc))
   EK[inc==0] <- 0

   inc.cl <- attr(inc.basic,"inc.cl")
   max.inc <- attr(inc.basic,"max.inc")
   min.inc <- attr(inc.basic,"min.inc")

   T.INC <- tapply(dis.basic,inc.basic,function(x)x)

   EK[inc==inc.cl[1] & !is.na(inc)] <- -rexp(table(inc==inc.cl[1])[2],rate = 1/6000)

      # für die negativen Einkommen
 
   for (i in 2:(length(inc.cl)-1)){
      n <- length(T.INC[[i]])

      min.tinc <- min.inc[i]
      max.tinc <- max.inc[i]

      max.hist <- max(hist(T.INC[[i]],breaks=100,plot=F)$counts)

      x1 <- rep(NA,n)
      f1 <- ff(T.INC[[i]])

      for (j in 1:n) {
       x1[j] <- AC.RE1()
      }


      EK[inc==inc.cl[i] & !is.na(inc) ] <- x1
   }

      # die Besetzung in der obersten Klasse
   Bes <- table(inc==inc.cl[i+1])
   if (length(Bes)>1){
      EK[inc==inc.cl[i+1]& !is.na(inc)] <- rexp(Bes[2],rate = 1/10000)+ min.inc[i]
   }
   return(EK)

}



   # new household ID

get.ind.reg <- function(x,y){
   x.n <- paste(y,x,sep="")
   lfn.x <- 1:length(x)
   lfn.x <- order(lfn.x)
   x.n <- rep(1:length(table(x.n)),table(x.n))
   x.n <- x.n[lfn.x]
   return(x.n)
}


#####################################################################################################
# Daten einladen
#####################################################################################################

load("/eingabe/Ameland/DSGpresentation/DSG.INC.HY010.i.RData")



INC <- generate.inc(HY010)   

INC.d <- get.dis.Inc3(inc=INC,inc.basic=INC,dis.basic=HY010)



setwd("/eingabe/Ameland/universe/Phase1/")

   # Diese beiden Einkommenskomponenten PY010G und PY050G müssen addiert werden. 
load("PY010G.cs.RData")
load("PY050G.cs.RData")

load("PB020.cs.RData")
load("PE010.cs.RData")
load("PE040.cs.RData")

load("PX030.cs.RData")

load("RX030.cs.RData")
load("RB020.cs.RData")


load("PB030.cs.RData")

load("DB020.cs.RData")
load("DB030.cs.RData")

   # Haushaltseinkommen
load("HY010.cs.RData")
load("HB020.cs.RData")
load("HB030.cs.RData")

   # Zuordnung des Haushaltseinkommens zu einzelnen Personen

RX030 <- get.ind.reg(RX030,RB020)
HB030 <- get.ind.reg(HB030,HB020)
PX030 <- get.ind.reg(PX030,PB020)
HB030 <- get.ind.reg(HB030,HB020)
DB030 <- get.ind.reg(DB030,DB020)

   # Zuordnung des Haushaltseinkommens zu den Personen
rh <- match(RX030,HB030)
ph <- match(PX030,HB030)
hd <- match(HB030,DB030)



HY010.r <- HY010[rh]

   # Das muss mit dem p Datensatz gemacht werden, da alle Einkommenskomponenten für p abgelegt sind.
HY010.p <- HY010[ph]


mod1 <- lm(HY010.p ~ PY050G + PY010G)

summary(mod1)

   # welchen Anteil hat die Einkommenskomponente Einkommen aus Erwerbsarbeit auf das Haushaltseinkommen

   # PY050G schwankt zwischen -290223.3 und 1e+06

PY050G2PX030 <- tapply(PY050G,PX030,function(x)x)

# PY050G2PX030

   # meistens hat nur eine Person im Haushalt "Cash benefits or losses from self-employment"

PY010G2PX030 <- tapply(PY010G,PX030,function(x)x)

   # Gibt es eine Person die Einkünfte aus 

   # gross employee cash or near cash income - 

PY010G2PX030sum <- tapply(PY010G,PX030,function(x)sum(x>0,na.rm=T))

PY010G2PX030sum2 <- tapply(PY010G,PX030,function(x)sum(x,na.rm=T))



mod2 <- lm(HY010 ~ PY010G2PX030sum2)

summary(mod2)

info.INC1 <- PY010G2PX030sum2[PY010G2PX030sum2>0]/HY010[PY010G2PX030sum2>0]

summary(info.INC1)

   # Anzahl der abhängigen Erwerbstätigen ist größer 0.

mod3 <- lm(HY010[PY010G2PX030sum>0] ~ PY010G2PX030sum2[PY010G2PX030sum>0])
summary(mod3)

info.INC2 <- PY010G2PX030sum2[PY010G2PX030sum>0]/HY010[PY010G2PX030sum>0]

summary(info.INC2)

   # Einkommen aus abhängiger Erwerbstätigkeit kann nicht kleiner als null sein.  
summary(PY010G)

mod4 <- lm(HY010[PY010G2PX030sum>0] ~ PY010G2PX030sum2[PY010G2PX030sum>0]*PY010G2PX030sum[PY010G2PX030sum>0])

summary(mod4)

   # mehrere Haushaltstypen werden generiert um das Haushaltseinkommen auf die einzelnen Personen und die Einkommenskomponenten zu verteilen.


   # Haushaltstyp 1: Einpersonenhaushalte mit einem Erwerbstätigen in abhängiger Stellung und nicht in Ausbildung.

HHT <- rep(0,length(HB030))
HHG <- tapply(PX030,PX030,length)
EDU2HHG <- tapply(PE010==1,PX030,sum)


HHT[HHG==1 & PY010G2PX030sum==1 & EDU2HHG==0] <- 1

   # Haushaltstyp 2: Einpersonenhaushalte mit Person in Ausbildung.

HHT[HHG==1 &  EDU2HHG==1] <- 2

   # Haushaltstyp 3: Zweipersonenhaushalte mit zwei Personen in abhängiger Stellung.

HHT[HHG==2 & PY010G2PX030sum==2 & EDU2HHG==0] <- 3





   # Befindet sich die Person momentan in der Ausbildung

sum(PY010G>0 & PE010==1,na.rm=T)

   # Es gibt 10768 Personen die sich in der Ausbildung  und gleichzeitig in abhängiger Erwerbstätigkeit befinden.
   # Allerdings ist der Verdienst dieser Personen minimal.


mean(PY010G>0 & PE010==1,na.rm=T)
max(PY010G>0 & PE010==1,na.rm=T)


mean(PY010G,na.rm=T)

   # PY030G ist erst ab 2007 im Datensatz enthalten

ind <- c("020","050","070","090","100","110","120","130","140")

eval(parse(text=paste("load('PY",ind,"G.cs.RData')",sep="",collapse=";")))

ind <- c("030","040","050","060","070","080","090","110","100")

eval(parse(text=paste("load('HY",ind,"G.cs.RData')",sep="",collapse=";")))


PY010G2HHsum <- tapply(PY010G,PX030,sum,na.rm=T)
PY020G2HHsum <- tapply(PY020G,PX030,sum,na.rm=T)
PY050G2HHsum <- tapply(PY050G,PX030,sum,na.rm=T)
PY070G2HHsum <- tapply(PY070G,PX030,sum,na.rm=T)
PY090G2HHsum <- tapply(PY090G,PX030,sum,na.rm=T)
PY100G2HHsum <- tapply(PY100G,PX030,sum,na.rm=T)
PY110G2HHsum <- tapply(PY110G,PX030,sum,na.rm=T)
PY120G2HHsum <- tapply(PY120G,PX030,sum,na.rm=T)
PY130G2HHsum <- tapply(PY130G,PX030,sum,na.rm=T)
PY140G2HHsum <- tapply(PY140G,PX030,sum,na.rm=T)


P.comp <- sum(PY010G2HHsum,PY020G2HHsum,PY050G2HHsum,PY070G2HHsum,PY090G2HHsum,PY100G2HHsum,PY110G2HHsum,PY120G2HHsum,PY130G2HHsum,PY140G2HHsum,na.rm=T)

# P.comp + eval(parse(text=paste("HY",ind[-length(ind)],"G",sep="",collapse="+")))- HY100G


sum(P.comp,HY030G,HY040G,HY050G,HY060G,HY070G,HY080G,HY090G,HY110G,na.rm=T)- HY100G


   # Modell für HHT==1

mod4 <- lm(PY[HHT==1]~HY010[HHT==1] + DB020[HHT==1])

   # Wenn die Effekte der Länder integriert werden resultiert ein besserer Fit
mod4.1 <- lm(PY[HHT==1]~HY010[HHT==1] )

PE040.HH <- tapply(PE040,PX030,function(x)x[1])


   # das persönliche Einkommen muss auf jeden haushalt summiert werden
PY010G2HHsum <- tapply(PY010G,PX030,sum,na.rm=T)

   # names(PY010G2HHsum)

ind <- match(HB030,names(PY010G2HHsum))

PY <- PY010G2HHsum[ind]


mod5 <- lm(PY[HHT==1]~HY010[HHT==1] + DB020[HHT==1] + PE040.HH[HHT==1])


############################################################
# Einkommenskomponenten auf dem lokalen Rechner visualisieren

setwd("/home/kolb/Desktop/Ameli/AmeLand/data/EU_SILC_2005")

load("DB100.cs.RData")


RX030 <- get.ind.reg(RX030,RB020)

ind <- c("020","050","070","090","100","110","120","130","140")

eval(parse(text=paste("load('PY",ind,"G.cs.RData')",sep="",collapse=";")))

ind <- c("030","040","050","060","070","080","090","110","100")

eval(parse(text=paste("load('HY",ind,"G.cs.RData')",sep="",collapse=";")))



# ord <- order(HY010)
# HY <- HY010[ord]
# PY <- PY[ord]
# plot(PY,ylim=c(0,200000),col="red")
# points(HY,type="l",ylim=c(0,200000))

plot(info.o[,1],xlim=c(0,2500),ylim=c(0,250000), main="HY010 und PY010 für Österreich 2005",col="grey")
DB100.h <- DB100[hd]

info <- cbind(PY[DB100.h==1],HY010[DB100.h==1])
info <- cbind(PY[HB020=="AT"],HY010[HB020=="AT"])

colnames(info) <- c("P.INC","H.INC")

plot(info,xlim=c(0,250000),ylim=c(0,250000), main="HY010 und PY010 für Österreich 2005",col="grey")

info.o <- info[order(info[,1]),]
plot(info.o[,1],xlim=c(0,6000),ylim=c(0,250000), main="HY010 und PY010 für Österreich 2005",col="blue",type="l")
points(info.o[,2],col="red")

   # Haushaltseinkommen in Klassen einteilen

H.INC.cl <- cut(HY010,breaks=quantile(HY010,probs=seq(0,1,0.05),na.rm=T),labels=F)

mean.HY <- tapply(HY010,H.INC.cl,mean)
median.HY <- tapply(HY010,H.INC.cl,median)


boxplot(PY~H.INC.cl,ylim=c(0,250000))
points(median.HY,col="blue")

   # Degree of urbanisation spielt kaum eine Rolle
boxplot(PY~DB100.h,ylim=c(0,250000))

   # Das addierte Einkommen aus abhängiger Erwerbstätigkeit macht ungefähr die Hälfte des Haushaltseinkommens aus. 
   # Dieser Einfluss unterscheidet sich nur geringfügig zwischen den unterschiedlichen Urbanisierungsgraden.
boxplot(PY/HY010~DB100.h,ylim=c(0,2))


   # Dafür ist der Unterschied zwischen den Ländern im EU-SILC Datensatz relevant
jpeg("/home/kolb/Desktop/Doktorarbeit/Text/graphics/boxplot.PY2HB020.jpeg")
   boxplot(PY/HY010~HB020,ylim=c(-1,4),xlab="Land",ylab="PY010.add / HY010",col="blue")
dev.off()





############################################################
# SimPopulation Paket


   # Diese Berechnungen müssen lokal durchgeführt werden, da das Paket SimPopulation nicht auf dem Server ist

   # Erzeugen des eusilcS Datensatzes

setwd("/home/kolb/Desktop/Ameli/AmeLand/data/EU_SILC_2005")

load("RB080.cs.RData")

AGE <- Year-RB080 
ACL <- cut(AGE,breaks=c(-1,16,30,40,65,96),labels=F)

eusilcS <- data.frame(ACL)

load("DB100.cs.RData")
load("RB030.cs.RData")
load("PB030.cs.RData")
load("RX030.cs.RData")
load("DB030.cs.RData")
load("HB030.cs.RData")

rp <- match(RB030,PB030)
rd <- match(RX030,DB030)
rh <- match(RX030,HB030)


DOU <-DB100[rd]

eusilcS$DOU <- DOU

   # Variable FST

  # family status - marital status
load("PB190.cs.RData")
FST <- PB190[rp]
FST[is.na(FST)]<-1
eusilcS$FST <- FST


  # Household size 


get.ind.reg <- function(x,y){
   x.n <- paste(y,x,sep="")
   lfn.x <- 1:length(x)
   lfn.x <- order(lfn.x)
   x.n <- rep(1:length(table(x.n)),table(x.n))
   x.n <- x.n[lfn.x]
   return(x.n)
}


load("RX030.cs.RData")
load("RB020.cs.RData")

RX030 <- get.ind.reg(RX030,RB020)
lfn <- 1:length(RX030)
lfn <- lfn[order(RX030)]
HHG <- tapply(RX030,RX030,length)
HHG <- rep(HHG,HHG)
HHG <- HHG[order(lfn)]

eusilcS$HHG <- HHG

   # VARIABLE PE010

load("PE010.cs.RData")
PE010 <- PE010[rp]

eusilcS$PE010  <- PE010 

   # VARIABLE PE040


# load("PE040.cs.RData")
# PE040 <- PE040[rp]

# eusilcS$PE040  <- PE040 

   # person with highest income erzeugen
load("EUS05.PWHI.RData")
eusilcS$PWHI  <- PWHI

load("RB200.cs.RData")
eusilcS$RB200  <- RB200

load("RB210.cs.RData")
eusilcS$RB210  <- RB210

   # Region erstellen
middle <- c("DE","FR","AT","LU","BE","CH")
nord <- c("DK","NL","IE","FI","SE","UK","IS","NO")
sud <- c("GR","ES","IT","PT","CY","MT")
east <- c("BG","CZ","EE","HU","LV","LT","PL","RO","SK","SI")
REG <- rep(NA,length(RB020))
reg.na <- match(RB020,middle)
REG[!is.na(reg.na)] <- 1
reg.na <- match(RB020,nord)
REG[!is.na(reg.na)] <- 2
reg.na <- match(RB020,sud)
REG[!is.na(reg.na)] <- 3
reg.na <- match(RB020,east)
REG[!is.na(reg.na)] <- 4
eusilcS$REG  <- REG

load("RB090.cs.RData")
SEX <- RB090

eusilcS$SEX  <- SEX

xvars <- colnames(eusilcS)

for (i in 1:length(xvars)){
   eval(parse(text=paste("load('/media/disk/Wien Daten/SynthDat/AML.",xvars[i],".RData')",sep="")))
}   

   # HHG stimmt noch nicht

lfn <- 1:length(HID)
lfn <- lfn[order(HID)]
HHG <- rep(table(HID),table(HID))
HHG <- HHG[order(lfn)]

   # die Variable PE040 ist zu kurz


# load('/media/disk/Wien Daten/SynthDat/AML.pe040.RData' )



eval(parse(text=paste("eusilcP <- data.frame(",xvars[1],")",sep=""))) 


for (i in 2:length(xvars)){
   eval(parse(text=paste("eusilcP$",xvars[i]," <-",xvars[i],sep="")))
}   


load("PY010G.cs.RData")

PY010G <- PY010G[rp]

eusilcS$PY010G <- PY010G

   # die relevante Fkt müsste simContinuous sein


   # eusilcS ist die Sample Distribution, in meinem Fall also EUSILC05

    # multinomial model with random draws
basic1 = colnames(eusilcP1)

eusilcM <- simContinuous(eusilcS1, eusilcP1,strata=REG,basic = basic1, upper = 200000)


hist(HY010,freq=F,ylim=c(0,2e-05),xlim=c(-100000,2e+05),col="green",xlab="Haushaltseinkommen",ylab="Dichte",breaks=1000,main="Histogramm für das Haushaltseinkommen")

###################################################################################

?simContinuous

   # das BSP im Paket:
     set.seed(1234)  # for reproducibility
     data(eusilcS)   # load sample data
     eusilcP <- simStructure(eusilcS)
     eusilcP <- simCategorical(eusilcS, eusilcP)
     basic <- c("age", "rb090", "hsize", "pl030", "pb220a")

     # multinomial model with random draws
     eusilcM <- simContinuous(eusilcS, eusilcP, 
         basic = basic, upper = 200000)

set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
eusilcP <- simStructure(eusilcS)
eusilcP <- simCategorical(eusilcS, eusilcP)
basic <- c("age", "rb090", "hsize", "pl030", "pb220a")

###################################################################################


# Einkommen generieren mit dem Paket der Wiener

   # anderer Ansatz

   # zunächst mal Logit Modell schätzen um zu sagen, ob Einkommen vorhanden oder nicht


par(mfrow=c(1,2))
hist(HY010,freq=F,ylim=c(0,2e-05),xlim=c(-100000,2e+05),col="green",xlab="Haushaltseinkommen",ylab="Dichte",breaks=1000,main="Histogramm für das Haushaltseinkommen")
hist(INC.d,freq=F,ylim=c(0,2e-05),xlim=c(-100000,2e+05),col="green",xlab="Haushaltseinkommen",ylab="Dichte",breaks=1000,main="Histogramm für das Haushaltseinkommen")

plot(density(HY010,na.rm=T),ylim=c(0,5e-05),xlim=c(-100000,3e+05),col="green")
abline(v=attr(INC,"max.inc"),col="purple")


plot(quantile(INC.d,probs=seq(0,1,by=0.001),na.rm=T),col="red",type="l")
points(quantile(HY010[-which(HY010>200000)],probs=seq(0,1,by=0.001),na.rm=T),col="green",type="l")

###############################

plot(density(inc[INC==5]))

dens <- density(inc[INC==5])

fu <- approxfun(x=dens$x,y=dens$y)
x <- seq(min(inc[INC==5]),max(inc[INC==5]),1)

points(x,fu(x),type="l")

