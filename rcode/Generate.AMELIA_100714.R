############################################
#
#  AMELIA
#     Projekt: AMELI
#     Author: Jan-Philipp Kolb
#       - Version 4:   - Function to categorize income 
#                      - HHG is generated
#
# - Update 11.11.2014
############################################

#----------------------------------------------------------#
# Libraries
#----------------------------------------------------------#

library(poLCA)

#----------------------------------------------------------#
# General Settings
#----------------------------------------------------------#

version <- "v1"

Year <- 2005

EUS.path <- paste("/eingabe/Ameland/universe/Phase1/",Year,sep="")

AML.path <- paste("/eingabe/Ameland/universe/results_",version,"/",sep="")

new.path <- "/eingabe/Ameland/universe/results_v2"

#--------------------------------------------------------------------------#
# Load functions
#--------------------------------------------------------------------------#


generate.inc <- function(x,ind,classes){
   quant.mpos <- quantile(x[ind], probs = seq(0, 1, 1/classes),na.rm=T,include.highest=T)
   quant.mpos[1] <- 0
   quant.mpos[length(quant.mpos)] <- max(x,na.rm=T)+1
   INC <- rep(0,length(x))


   INC <- cut(x[ind],breaks=quant.mpos,labels=F)
   # inc.cl <- sort(unique(INC))
   # min.inc <- tapply(x[ind],INC,min)
   # max.inc <- tapply(x[ind],INC,max)

   # attr(INC,"min.inc")<- min.inc
   # attr(INC,"max.inc")<- max.inc
   # attr(INC,"inc.cl")<- inc.cl
   return(INC)
}




get.dis.Inc3 <- function(inc,inc.basic,dis.basic){


   ff<-function(x){
      dens <- density(x)
      return(approxfun(x=dens$x,y=dens$y))
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

      # für die negativen Einkommen
 
#   for (i in 1:(length(inc.cl)-1)){
   for (i in 1:(length(inc.cl))){

      n <- sum(inc==inc.cl[i],na.rm=T)

      min.tinc <- min.inc[i]
      max.tinc <- max.inc[i]

      # max.hist <- max(density(T.INC[[i]],na.rm=T)$y)
      max.hist <- max(density(T.INC[[i]])$y)

      x1 <- rep(NA,n)
      f1 <- ff(T.INC[[i]])

      for (j in 1:n) {
       x1[j] <- AC.RE1()
         #cat(j,"\n")
      }


      EK[inc==inc.cl[i]] <- x1
   }

      # die Besetzung in der obersten Klasse
   # Bes <- table(inc==inc.cl[i+1])
   # if (length(Bes)>1){
   #   EK[inc==inc.cl[i+1]& !is.na(inc)] <- rexp(Bes[2],rate = 1/10000)+ min.inc[i]
   # }
   return(EK)

}


get.ind.reg <- function(x,y){
   x.n <- paste(y,x,sep="")
   lfn.x <- 1:length(x)
   lfn.x <- order(lfn.x)
   x.n <- rep(1:length(table(x.n)),table(x.n))
   x.n <- x.n[lfn.x]
   return(x.n)
}


region.alloc <- function(regions,o.vec){
   eval(parse(text=paste("region <- rep(NA,length(",paste(o.vec),"))",sep="")))

   for (i in 1:length(regions)){   
      eval(parse(text=paste("reg.na <- match(",paste(o.vec),",",regions[i],")",sep="")))
      region[!is.na(reg.na)] <- i
   }
   return(region)
}

   # Eigentliche Erzeugung

Synth.prep <- function(yvarF,varsF,o.indF,s.indF,region,gew)
## yvarF = yvariable
## varsF = X variablen namen
## o.indF = vector MZ Missings
## s.indF = vector Synth Missings
## region=character Kreis
{

   # Hier müssen mindestens zwei Hilfsvariablen eingefügt werden, sonst funktioniert der Algorithmus nicht

      # Mikrozensusvariablen werden in ein Data Frame zusammen gefügt.
      # Achtung bei der EU-SILC Erzeugung wird das anders gemacht
   # eval(parse(text=paste("MZ<-data.frame(",paste(varsF,"=O.dat$",varsF,sep="",collapse=","),")",sep="")))

   if(!missing(s.indF)){
      Synth<-SynthDat[s.indF,]
   }else{
      Synth<-SynthDat
   }


   if (!missing("o.indF")){
      # o.indF = o.indF & as.character(KRS)==region
      # das muss im Vergleich zu der Generierung im Zensus verändert werden
      MZ<-MZ[o.indF,]
      gew <- gew[o.indF]
      yvarF<-yvarF[o.indF]
      #}else{
      #o.indF = as.character(KRS)==region
      #MZ<-MZ[o.indF,]
      #gew <- gew[o.indF]
      #yvarF<-yvarF[o.indF]

   }

   tab3<-data.frame(MZ[,eval(parse(text=paste("c(",paste("'",varsF,"'",collapse=",",sep=""),")",sep="")))],weights=gew,y=yvarF)
   rbi3<-as.data.frame(prop.table(xtabs(paste("formula(weights ~ y + ",paste(names(tab3)[-c(length(names(tab3))-c(0,1))],collapse="+",sep=""),")",sep=" "),data=tab3)))

   #Vorsicht wenn Variable mehr als 9 ausprägungen hat..
   eval(parse(text=paste("res.vars3 <- paste(",paste(paste("rbi3$",varsF,sep="",collapse=","),collapse=","),",sep='')",sep="")))
      # we do have problems with NA's in the samp.res function

     # hier können Stichprobennulen oder strukturelle Nullen auftauchen, Felde$
      # gekennzeichnet werden

      # diese Zeile ist aber notwendig, da sonst nicht gezogen werden kann
   rbi3$Freq[is.na(rbi3$Freq)|rbi3$Freq==0]<-0.00000001
   rbi3$Freq<-   rbi3$Freq/sum(   rbi3$Freq)
   eval(parse(text=paste("res.syn3 <-paste(", paste("Synth$",varsF,sep="",collapse=","),",sep='')",sep="")))

   tmp<-   res.vars3%in%res.syn3
   tmp3 <- rbi3[tmp,]
   tmp3$Freq <- tmp3$Freq/sum(tmp3$Freq)
   tab<-table(res.syn3)
   o<-sapply(as.list(names(tab)),function(x){
      which(x==res.vars3)
   },simplify = F)
   o1<-sapply(o,length)==0
   tab2<-tab[!o1]
   if(length(tab2)==0)return(rep(NA,length(Synth[,1])))
   o2<-o[!o1]

   lol<-mapply(function(a,b){
            sample(rbi3[,"y"][a],b,T,rbi3[,"Freq"][a])
   },o2,as.list(tab2),SIMPLIFY=F)

   y.hat<-as.character(rep(NA,length(Synth[,1])))
   ind<-match(res.syn3,names(tab2))
   ind.na<-!is.na(ind)

   #   for(i in unique(ind[ind.na])){
   #         y.hat[ind%in%i]<-as.vector(lol[[i]])
   #   }
   ####
   res.syn4 <- res.syn3[ind.na]
   lfn <- 1:length(res.syn4)
   lfn <- lfn[order(res.syn4)]
   ord <- order(names(tab2))
   erg <- unlist(lol[ord])
   erg <- erg[order(lfn)]
   y.hat[ind.na==T]<-as.character(erg)
   y.hat <- as.factor(y.hat)
   return(y.hat)
}

   # Das Einkommen wird wieder verstetigt. 

dis.inc.agg <- function(x,x.sbasic,x.ind,x.inc,n.ex){
      # x ist das Einkommen aus dem Originaldatensatz
      # x.sbasic ist das Einkommen aus dem synthetischen Datensatz in kategorialer Form      
      # x.ind sind Regeln für den Originaldatensatz, welche Werte dem synthetischen stetigen Einkommen zu Grunde liegen sollen
      # x.inc ist das kategoriale Einkommen des Originaldatensatzes
      # n.ex ist die Zahl der Extremwerte, die entfernt werden soll

   if(is.logical(x.ind)==F)warning("x.ind should be logical")

   x.sbasic[is.na(x.sbasic)] <- 0

      # Extremwerte werden heraus genommen
   ind <- order(x,decreasing=T)[1:n.ex]

   x.ind[ind] <- F  
      # Nullwerte werden heraus genommen
   x.ind[x==0] <- F

   INC <-x.inc[x.ind]

      # folgende Form wurde so beibehalten, deshalb sind diese Operationen noch etwas komplex
   inc.cl <- sort(unique(INC))
   min.inc <- tapply(x[x.ind],INC,min,na.rm=T)
   max.inc <- tapply(x[x.ind],INC,max,na.rm=T)
   attr(INC,"min.inc")<- min.inc
   attr(INC,"max.inc")<- max.inc
   attr(INC,"inc.cl")<- inc.cl

   xv <- rep(NA,length(x.sbasic))
   xv[x.sbasic!=0] <- get.dis.Inc3(inc=x.sbasic[x.sbasic!=0],inc.basic=INC,dis.basic=x[x.ind])
   xv[x.sbasic==0] <- 0
   return(xv)
}





quantile.korr <- function(x,acuteness){
   xv <- gsub("G","D",x)   
   eval(parse(text=paste("max.xs <- max(",xv,",na.rm=T)")))
   eval(parse(text=paste("xo <- ",x,"[",x,">0 & ",x,"<max.xs & !is.na(",x,")]",sep="")))
   vec1<-cut(xo,breaks=quantile(xo,probs=seq(0,1,1/acuteness)),labels=F)
   eval(parse(text=paste("xs <- ",xv,"[",xv,">0 & ",xv,"<max.xs & !is.na(",xv,")]",sep="")))
   vec2<-cut(xs,breaks=quantile(xs,probs=seq(0,1,1/acuteness)),labels=F,include.lowest=T)
   xom <- tapply(xo,vec1,mean)
   xsm <- tapply(xs,vec2,mean)

      # hier ändert sich die Zusammensetzung je nachdem welcher Wert größer ist
   quot <- (xsm-xom)/xsm

   ind <- match(vec2,names(quot))
   xs3 <- xs - quot[ind]*xs
   eval(parse(text=paste(xv,"[",xv,">0 & ",xv,"<max.xs & !is.na(",xv,")] <- xs3",sep="")))
   eval(parse(text=paste("return(",xv,")")))
}



quantile.korr2 <- function(x){
   xv <- gsub("G","D",x)   
   #eval(parse(text=paste("max.xs <- max(",xv,",na.rm=T)")))
   eval(parse(text=paste("xo <- ",x,"[ !is.na(",x,")]",sep="")))
   vec1<-cut(xo,breaks=quantile(xo,probs=seq(0,1,1/5)),labels=F)
   eval(parse(text=paste("xs <- ",xv,"[!is.na(",xv,")]",sep="")))
   vec2<-cut(xs,breaks=quantile(xs,probs=seq(0,1,1/5)),labels=F,include.lowest=T)
   xom <- tapply(xo,vec1,mean)
   xsm <- tapply(xs,vec2,mean)

      # hier ändert sich die Zusammensetzung je nachdem welcher Wert größer ist
   quot <- (xsm-xom)/xsm

   ind <- match(vec2,names(quot))
   xs3 <- xs - quot[ind]*xs
   eval(parse(text=paste(xv,"[ !is.na(",xv,")] <- xs3",sep="")))
   eval(parse(text=paste("return(",xv,")")))
}


   # diese Funktion dauert sehr lange, und muss noch geprüft werden
gen.vars <- function(x,varsF,MZ){

   eval(parse(text=paste("load('",x,"G.cs.RData')",sep="")))
   
   eval(parse(text=paste("xx <-", x,"G",sep=""))) 
   Synth<- rep(NA,nrow(SynthDat))

   HY.ind <- rep(F,length(xx))
   HY.ind[xx>0 & !is.na(xx)] <- T

   o.indF <-!is.na(xx)

   s.indF <- rep(T,nrow(SynthDat))
   HY.inc<-rep(0,length(HY.ind))
   HY.inc[HY.ind] <- generate.inc(x=xx,ind=HY.ind,classes=15)
  
   gew <- DB090/sum(DB090)

   HY.ic <- HY.inc + 1
   info <- Synth.prep(yvarF=HY.ic,varsF,o.indF,s.indF,region,gew)
   Synth[s.indF] <- as.numeric(as.character(info)) - 1

   conNA <- all(is.na(Synth))==F

   Anz <- length(varsF)

   for(i in Anz:3){

      s.indF <- is.na(Synth)
      varsF <- varsF[-i]
      info <- Synth.prep(yvarF=HY.ic,varsF,o.indF,s.indF,region,gew)
      Synth[s.indF] <- as.numeric(as.character(info)) - 1
   
      conNA <- all(is.na(Synth)==F)
      if(conNA==T)break
   }
   return(Synth)
}

   # Das Haushaltsäquivalenzeinkommen nach Stefan Zins
calc.evs <- function(x){
   y <- c(1,ifelse(x[-1]>14,.5,.3))
   return(y)
}


   # und so wirds wirklich gemacht:
calc.evs2 <- function(x){
   y <- sum(1,ifelse(x[-1]>14,.5,.3))
   return(y)
}



#--------------------------------------------------------------------------#
# Load data
#--------------------------------------------------------------------------#



   # Die schon erzeugten synthetischen Daten einlesen

setwd(AML.path)

   # Übersicht über die Variablen die schon abgespeichert wurden

Namen <- dir()
ind <- grep("AML.",Namen)

load("AML.AGE.RData")
load("AML.ACL.RData")
load("AML.SEX.RData")
load("AML.FST.RData")

load("AML.REG.RData")
# load("AML.CIT.RData")
load("AML.DOU.RData")
load("AML.HHG.RData")
load("AML.PE010.RData")
load("AML.PE040.RData")
# load("AML.PWHI.RData")
load("AML.RB200.RData")
load("AML.RB210.RData")
# load("AML.UEP.RData")
load("AML.FILE.RData")
load("AML.HHG.RData")
load("AML.HID.RData")


   # folgende Variablen wurden erst im Laufe dieses Jobs erstellt

load("AML.PB210.RData")
load("AML.PY010G.RData")
load("AML.PY090G.RData")
load("AML.PY020G.RData")

   # RB200: Residential status


   # RB210: Basic activity status
load("AML.RB210.RData")

SynthDat <- as.data.frame(SEX)
SynthDat$RB210 <- RB210
SynthDat$FST <- FST
SynthDat$ACL <- ACL
SynthDat$PE040 <- PE040
SynthDat$PE010 <- PE010

SynthDat$DOU <- DOU
SynthDat$REG <- REG
SynthDat$HHG <- HHG

SynthDat$AGE <- AGE

rm(HHG,DOU,ACL,FST,SEX)

   # folgende Variablen wurden erst im Laufe dieses Jobs erstellt

SynthDat$PY010G <- PY010G
SynthDat$PY090G <- PY090G
SynthDat$PY020G <- PY020G

   # Die Konstante soll es ermöglichen ein Generierungsmodell auch mit nur einer Hilfsvariable durchzuführen.
SynthDat$const <- rep(1,nrow(SynthDat))

   # Die Daten aus dem EU-SILC einladen

setwd(EUS.path)

load("HY030G.cs.RData")


load("RB050.cs.RData")
   # Geschlecht

load("RB210.cs.RData")


load("DB020.cs.RData")
load("DB030.cs.RData")

   # Das Haushaltsgewicht
load("DB090.cs.RData")

   # degree of urbanisation
load("DB100.cs.RData")

load("HB020.cs.RData")
load("HB030.cs.RData")


load("HY010.cs.RData")
load("HY040G.cs.RData")
load("HY050G.cs.RData")
load("HY070G.cs.RData")
load("HY080G.cs.RData")
load("HY090G.cs.RData") # Brutto Zinsen, Brutto-Dividenden, Brutto-Gewinne aus Kapitalanlagen in Unternehmen ohne eigene Rechtspersönlichkeit
load("HY100G.cs.RData")
load("HY110G.cs.RData") # Von Personen unter 16 Jahren erhaltenes Brutto-Einkommen
load("HY120G.cs.RData")
load("HY130G.cs.RData")


load("HX040.cs.RData")



load("PB020.cs.RData")
load("PB030.cs.RData")

   # Die Variable PB050 existiert nicht.
load("PB040.cs.RData")
load("PB140.cs.RData")
   # Geschlecht
load("PB150.cs.RData")
load("PB190.cs.RData")
   # Nationalität
load("PB210.cs.RData")

load("PE010.cs.RData")
load("PE040.cs.RData")


   # Imputierte Daten werden hier nicht verwendet, damit nicht zu viel von der Struktur verloren geht

load("PH010.cs.RData")
load("PH020.cs.RData")
load("PH030.cs.RData")
load("PH040.cs.RData")
load("PH050.cs.RData")


load("PL080.cs.RData")

load("PX030.cs.RData")

   # Einkommensdaten

load("PY010G.cs.RData")
load("PY020G.cs.RData")
# load("PY030G.cs.RData")

load("PY050G.cs.RData")
load("PY070G.cs.RData")

   # Transfers
load("PY090G.cs.RData")
load("PY100G.cs.RData")
load("PY110G.cs.RData")

load("PY120G.cs.RData")
load("PY130G.cs.RData")

load("PY140G.cs.RData")


load("RB020.cs.RData")
load("RB030.cs.RData")

   # Geburtsjahr
load("RB080.cs.RData")
load("RB090.cs.RData")
   # Alter 
load("RX010.cs.RData")
   # HaushaltsID
load("RX030.cs.RData")

#--------------------------------------------------------------------------#
# Match the data files
#--------------------------------------------------------------------------#

   # new household ID
   # das muss sowohl für die HaushhaltsID als auch für die Personen ID gemacht werden.

DB030 <- get.ind.reg(DB030,DB020)
HB030 <- get.ind.reg(HB030,HB020)
RX030 <- get.ind.reg(RX030,RB020)
PX030 <- get.ind.reg(PX030,PB020)
PB030 <- get.ind.reg(PB030,PB020)
RB030 <- get.ind.reg(RB030,RB020)

   # an erster Stelle stand vorher RB030 aber das dürfte die Personen ID sein
rp <- match(RB030,PB030)
rd <- match(RX030,DB030)
rh <- match(RX030,HB030)

   # an erster Stelle stand vorher PB030 aber das dürfte die Personen ID sein
pd <- match(PX030,DB030)
ph <- match(PX030,HB030)


hd <- match(HB030,DB030)

   # RB030: Personal ID 
   # folgender Befehl ist bewusst ausgeklammert, da das Matching so rum nicht funktioniert

# pr <- match(PB030,RB030)

p <- paste(PB020,PB030,sep="")
r <- paste(RB020,RB030,sep="")

pr <- match(p,r)

ind.ph.info <- tapply(PX030,PX030,function(x)x[1])

ind.ph <- match(HB030,ind.ph.info)

#--------------------------------------------------------------------------#
# Kovariaten vorbereiten
#--------------------------------------------------------------------------#



   # EU_SILC Variablen umbenennen

SEX <- RB090


   # degree of urbanisation
DOU <-DB100[rd]

FST <- PB190[rp]

DOU.p <-DB100[pd]

AGE.p <- Year - PB140

#   save(AGE.p,file="EUS.AGE.p.RData")
 
AGE.r <- Year - RB080

ACL <- cut(AGE.p,breaks=c(-1,16,30,40,65,96),labels=F)

save(ACL,file="EUS.ACL.RData")

HHG.p <- HX040[ph] 

RB210.p <- RB210[pr]

RX010.p <- RX010[pr]

#--------------------------------------------------------------------------#
# Regionsidentifikatoren vorbereiten
#--------------------------------------------------------------------------#


middle <- c("DE","FR","AT","LU","BE","CH")
nord <- c("DK","NL","IE","FI","SE","UK","IS","NO")
sud <- c("GR","ES","IT","PT","CY","MT")
east <- c("BG","CZ","EE","HU","LV","LT","PL","RO","SK","SI")


regions <- c("middle","nord","sud","east")

region <- region.alloc(regions,o.vec="RB020")
region.r <- region
region.p <- region.alloc(regions,o.vec="PB020")
region.h <- region.alloc(regions,o.vec="HB020")

#--------------------------------------------------------------------------#
# Berechnung des Haushaltsäquivalenzeinkommens vorbereiten
#--------------------------------------------------------------------------#



   # das Haushaltseinkommen auf Personenlänge ausdehnen

HID2HID <- tapply(HID,HID,function(x)x[1])

# HID2HID2 <- tapply(HID,HID,function(x)x)

ind.HID <- match(HID,HID2HID)

   # erste Person hat volles Einkommen

Gew2 <- tapply(AGE,HID,calc.evs2)
Gew2 <-unlist(Gew2)


#--------------------------------------------------------------------------#
# Haushaltstypen vorbereiten
#--------------------------------------------------------------------------#

#################
# HHT

# Die Bildung von Haushaltstypen soll dazu dienen eine gewisse Homogenität innerhalb der Haushalte zu erzeugen. Dies basiert auf der Annahme, dass sich die Personen in einem Haushalt gegenseitig beeinflussen. 

   # Die Haushalte werden nach Einkommensumfang eingeteilt, hierfür wird die Variable HY010G verwendet.

HID.u <- tapply(1:length(HID),HID,function(x)x)

   # der Haushaltstyp muss nur für eine Person im Haushalt erstellt werden.

   # an dieser Stelle fehlte das .u, deshalb könnte es hier zu Fehlern gekommen sein.
HID.u2 <- lapply(HID.u,function(x)x[1])

   # Die synthetische Ausgangsmatrix hat für diesen Fall nur den Umfang Anzahl der Haushalte.


ind=unlist(HID.u2)

SynthDat.HH <- as.data.frame(ind)

   # Die Region ist ebenfalls für alle Mitglieder eines Haushalts gleich

SynthDat.HH$REG <- REG[ind]

SynthDat.HH$HHG <- HHG[ind]

SynthDat.HH$DOU <- DOU[ind]


HY010.ind <- rep(F,length(HY010))
HY010.ind[!is.na(HY010)] <- T

HY010.inc <- rep(0,length(HY010))
HY010.inc[HY010.ind] <- generate.inc(x=HY010,ind=HY010.ind,classes=10)

o.indF <- rep(T,length(HY010))
o.indF[is.na(HY010)] <- F

gew=DB090/sum(DB090)

varsF=c("REG","HHG","DOU")
MZ <- data.frame(REG=region.h,DOU=DB100,HHG=HX040)

SynthDat.HH$HHT <- NA

SynthDat <- SynthDat.HH

SynthDat.HH$HHT[s.indF] <- Synth.prep(yvarF=HY010.inc,varsF,o.indF,s.indF,region,gew)


s.indF <- rep(T,nrow(SynthDat.HH))

   # in der Variable HHT 

varsF=c("REG","DOU")
MZ <- data.frame(REG=region.h,DOU=DB100)

s.indF <- is.na(SynthDat.HH$HHT)

SynthDat.HH$HHT[s.indF] <- Synth.prep(yvarF=HY010.inc,varsF,o.indF,s.indF,region,gew)

   # HHT muss jetzt wieder auf den gesamten Datensatz zurück gespielt werden

ind2 <- match(HID,HID.u2)

HHT <- SynthDat.HH$HHT[ind2]


SynthDat$HHT <- HHT

   # HHT muss noch auf Laenge der Personen gebracht werden

HY010.inc.p <- HY010.inc[ph]

#--------------------------------------------------------------------------#
# Weitere Personeninformationen
#--------------------------------------------------------------------------#

#################
# PB210

   # Der Datensatz sollte noch die Nationalität enthalten
   # Die Information über das Geburtsland (PB210) ist leider nur für den P-File vorhanden, für Personen unter 16 Jahren können also keine Angaben gemacht werden.

o.indF <- rep(T,length(PB210))
o.indF[is.na(PB210)] <- F

s.indF <- rep(T,nrow(SynthDat.HH))

gew=PB040/sum(PB040)

varsF=c("ACL","SEX","RB210","PE040","FST","DOU","HHG")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,HHG=HHG.p)

SynthDat$PB210 <- NA

SynthDat$PB210[s.indF] <- Synth.prep(yvarF=PB210,varsF,o.indF,s.indF,region,gew)

   # Nachdem die Erzeugung mit einer Vielzahl von Variablen durchgeführt wurde sind noch NA's enthalten, deshalb wird diese Erzeugung nun noch mal durchgeführt.

s.indF <- FILE=="P" & is.na(SynthDat$PB210)

varsF=c("ACL","SEX","RB210","FST","DOU","HHG")

MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,FST=PB190,DOU=DOU.p,HHG=HHG.p)

SynthDat$PB210[s.indF] <- Synth.prep(yvarF=PB210,varsF,o.indF,s.indF,region,gew)

SynthDat$PB210[SynthDat$PB210==2] <- "EU"
SynthDat$PB210[SynthDat$PB210==3] <- "LOC"
SynthDat$PB210[SynthDat$PB210==4] <- "OTH"
SynthDat$PB210[SynthDat$PB210==1] <- ""

#--------------------------------------------------------------------------#
# Variablen bezüglich Gesundheit
#--------------------------------------------------------------------------#

   # Klassifikationsverfahren durchführen

Hdat <- data.frame(PH010,PH020,PH030,PH040)

ind <- complete.cases(Hdat)

Hdat <- Hdat[ind,]
H.cla=5


f <- cbind(PH010,PH020,PH030,PH040)~1
Hclust <- poLCA(f,Hdat,nclass=H.cla,maxiter=1000)

Hcl <- Hclust$predclass

SynthDat$Hcl <- NA


o.indF <- ind
s.indF <- FILE=="P"
gew=PB040/sum(PB040)
varsF=c("ACL","SEX","RB210","PE040","FST","DOU")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p)
SynthDat$Hcl[s.indF] <- Synth.prep(yvarF=Hcl,varsF,o.indF,s.indF,region,gew)

   # Die Variable PE040 macht oft Probleme, deshalb wird diese hier weg gelassen

s.indF <- FILE=="P" & is.na(SynthDat$Hcl)
varsF=c("ACL","SEX","RB210","FST","DOU")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,FST=PB190,DOU=DOU.p)
SynthDat$Hcl[s.indF] <- Synth.prep(yvarF=Hcl,varsF,o.indF,s.indF,region,gew)


   # Die Variable PE040 macht oft Probleme, deshalb wird diese hier weg gelassen

s.indF <- FILE=="P" & is.na(SynthDat$Hcl)
varsF=c("ACL","SEX","RB210","FST","DOU")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,FST=PB190,DOU=DOU.p)
SynthDat$Hcl[s.indF] <- Synth.prep(yvarF=Hcl,varsF,o.indF,s.indF,region,gew)

#################
# PH010


SynthDat$PH010 <- NA


   # zur Sicherheit wird das Ergebnis der Latent Class Analysis zwischen gespeichert
Hcl1 <- Hcl

Hcl <- rep(NA,length(PH010))
Hcl[ind==T]<- Hcl1


o.indF <- ind
s.indF <- FILE=="P"
gew=PB040/sum(PB040)
varsF=c("ACL","SEX","RB210","FST","DOU","Hcl")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,Hcl)
SynthDat$PH010[s.indF] <- Synth.prep(yvarF=PH010,varsF,o.indF,s.indF,region,gew)


#################
# PH020

SynthDat$PH020 <- NA


o.indF <- ind
s.indF <- FILE=="P"
gew=PB040/sum(PB040)
varsF=c("ACL","SEX","RB210","FST","Hcl","PH010")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,Hcl,PH010)
SynthDat$PH020[s.indF] <- Synth.prep(yvarF=PH020,varsF,o.indF,s.indF,region,gew)

   # Die Generierung der Variablen funktioniert besser, wenn die Zahl der Ausprägungen auf der zu generierenden Variable nicht sehr hoch ist. In diesem Fall können auch viele Hilfsvariablen zur Erzeugung verwendet werden.

#################
# PH030

SynthDat$PH030 <- NA


o.indF <- ind
s.indF <- FILE=="P"
gew=PB040/sum(PB040)
varsF=c("ACL","SEX","RB210","FST","Hcl","PH010","PH020")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,Hcl,PH010,PH020)
SynthDat$PH030[s.indF] <- Synth.prep(yvarF=PH030,varsF,o.indF,s.indF,region,gew)

   # Bei der Erzeugung der Gesundheitsvariablen tauchen relativ wenig fehlende Werte auf.

#################
# PH040


SynthDat$PH040 <- NA


o.indF <- ind
s.indF <- FILE=="P"
gew=PB040/sum(PB040)
varsF=c("ACL","SEX","RB210","FST","Hcl","PH010","PH020","PH030")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,Hcl,PH010,PH020,PH030)
SynthDat$PH040[s.indF] <- Synth.prep(yvarF=PH040,varsF,o.indF,s.indF,region,gew)

#################
# PH050

   # Jan Seger hat folgende Variable auch noch für seine Indikatoren benötigt.


SynthDat$PH050 <- NA


o.indF <- ind
s.indF <- FILE=="P"
gew=PB040/sum(PB040)
varsF=c("ACL","SEX","RB210","FST","PH010","PH020","PH030","PH040")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,PH010,PH020,PH030,PH040)
SynthDat$PH050[s.indF] <- Synth.prep(yvarF=PH040,varsF,o.indF,s.indF,region,gew)



#--------------------------------------------------------------------------#
# Variablen bezüglich Arbeitslosigkeit
#--------------------------------------------------------------------------#

#################
# PL080


SynthDat$PL080 <- NA
ind <- rep(T,length(PL080))

   # für inaktive Personen muss Variable nicht erzeugt werden
ind[RB210%in%c(3,4)] <- F

o.indF <- ind
s.indF <- FILE=="P" & SynthDat$RB210!=3 & SynthDat$RB210!=4
gew=PB040/sum(PB040)


   # Desto mehr Kovariaten zur Erzeugung verwendet werden, desto länger dauert auch die Erzeugung.
   # Eventuell liegt die hohe Zeit, die die Erzeugung dieser Variablen verbraucht auch daran, dass die Zielvariable so viele Ausprägungen hat.

varsF=c("ACL","SEX","RB210","Hcl","REG","PY010G")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,Hcl,REG=region.p,PY010G)
SynthDat$PL080[s.indF] <- Synth.prep(yvarF=PL080,varsF,o.indF,s.indF,region,gew)





#--------------------------------------------------------------------------#
# Die Einkommenskomponenten aus unselbstständiger Erwerbstätigkeit
#--------------------------------------------------------------------------#

#################
# PY010G


SynthDat$PY010G <- NA


   # zuerst mal die Nulleinkommen bestimmen

   # Im synthetischen Datensatz kann eine Person nur ein Einkommen aus Erwerbsarbeit haben, wenn sie auch erwerbstätig ist.

SynthDat$PY010G[SynthDat$RB210!=1] <- 0

SynthDat$PY010G[SynthDat$RB210==1] <- 1


PY010G.ind <- rep(F,length(PY010G))

   # Wichtig 
PY010G.ind[RB210.p==1 & PY010G>1] <- T


PY010G.inc<-rep(0,length(PY010G.ind))
   
   # hier gehen die erzeugten Attribute leider verloren
PY010G.inc[PY010G.ind==1] <- generate.inc(x=PY010G,ind=PY010G.ind,classes=10)




o.indF <- RB210.p==1 & PY010G>1
o.indF[is.na(o.indF)] <- F

s.indF <- SynthDat$PY010G==1


gew=PB040/sum(PB040)

varsF=c("ACL","SEX","RB210","PE040","FST","DOU")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p)

# SynthDat$PY010G <- Synth.prep(yvarF=PY010G.inc,varsF,o.indF,s.indF,region,gew)

SynthDat$PY010G[s.indF] <- Synth.prep(yvarF=PY010G.inc,varsF,o.indF,s.indF,region,gew)

   # Es gibt Kombinationen, die nicht vorhanden sind, die Zahl der fehlenden Werte kann nicht gesenkt werden, wenn 
   # die Anzahl der Einkommensklassen gesenkt wird.

   # Diesen Abschnitt kann man mit dem Prinzip der Mikroaggregation vergleichen

s.indF <- is.na(SynthDat$PY010G) 
# s.indF[SynthDat$PY010G==0] <- F

   # Die NA's können eigentlich nur da auftauchen, wo kein Nulleinkommen vorliegt.

   # Bei der Variable PE040 gibt es Probleme, diese hat für die synthetische Populationen, die bei der originalen 
   # Population nicht vorhanden sind.

varsF=c("ACL","SEX","RB210","FST","DOU")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,FST=PB190,DOU=DOU.p)



SynthDat$PY010G[s.indF] <- Synth.prep(yvarF=PY010G.inc,varsF,o.indF,s.indF,region,gew)

   # Das weglassen von drei Variablen bei der Erzeugung hat nichts verändert an der Zahl der fehlenden Werte.


#################
# PY020G


   # PY020G ist das Non cash income, hier gelten die gleichen Regeln, wie für das cash Income

   # Vektor für die synthetische Population vorbereiten
SynthDat$PY020G <- NA

   # Achtung, wenn das hier noch mal durchgeführt werden soll, dann muss darauf geachtet werden, das die Variable RB210 auch wirklich in der SynthDat Umgebung ist.
SynthDat$PY020G[SynthDat$RB210!=1] <- 0
SynthDat$PY020G[SynthDat$RB210==1] <- 1


PY020G.ind <- rep(F,length(PY020G))
PY020G.ind[RB210.p==1 & PY020G>1] <- T
PY020G.inc<-rep(0,length(PY020G.ind))

PY020G.inc[PY020G.ind==1] <- generate.inc(x=PY020G,ind=PY020G.ind,classes=10)


o.indF <- rep(F,length(PY020G))
o.indF[RB210.p==1] <- T
o.indF[is.na(o.indF)] <- F

s.indF <- SynthDat$PY020G==1
gew=PB040/sum(PB040)



varsF=c("PY010G","const")
MZ <- data.frame(PY010G=PY010G.inc,const=rep(1,length(PB150)))

   # Wenn diese Variable so wie PY010G erzeugt wird resultieren sehr viele Nullzellen. Bei dieser Variable sind die Fallzahlen in den Klassen sehr gering, 
   # werden weniger Klassen verwendet werden dennoch nur fehlende Werte erzeugt. 


   # Was nun folgt ist ein kleiner unschönder Trick.
   # mit der Funktion Synth.prep können keine Ausprägungen von null erzeugt werden
   # bei dieser Variable ist aber die richtige Wiedergabe der Nullkomponenten sehr wichtig

   # Achtung: wenn das zwei mal hintereinander durchgeführt wird muss man auch die 2 wieder abziehen.

PY020G.inc <- PY020G.inc+1

info <- Synth.prep(yvarF=PY020G.inc,varsF,o.indF,s.indF,region,gew)

SynthDat$PY020G[s.indF] <- as.numeric(as.character(info)) -1

   # Pablos Idee war es hier ein multinomiales Modell zu verwenden


#################
# PY030G

#  Employer's social insurance contribution

# SynthDat$PY030G[SynthDat$RB210!=1] <- 0
# SynthDat$PY030G[SynthDat$RB210==1] <- 1


   # Die Höhe dieser Einkommenskomponenten wird später noch bestimmt


   # Bei den Sozialbeiträge der Arbeitgeber (PY030G) wird der Nettobetrag verwendet, da der Bruttobetrag nicht im Datensatz verfügbar ist.
# SynthDat$PY030G <- NA

#--------------------------------------------------------------------------#
# Die Einkommenskomponenten aus Transferleistungen
#--------------------------------------------------------------------------#

#################
# PY090G

   # Negativeinkommen, kommen so selten vor, dass diese im synthetischen Datensatz ausgeschlossen werden.
# sum(PY090G<0,na.rm=T)


   # Arbeitslosenunterstützung
SynthDat$PY090G <- NA

SynthDat$PY090G[SynthDat$RB210==2] <- 1
SynthDat$PY090G[SynthDat$RB210 %in% c(1,3,4)] <- 0


PY090G.ind <- rep(F,length(PY090G))
PY090G.ind[RB210.p==2 & PY090G>1] <- T
PY090G.inc<-rep(0,length(PY090G.ind))

PY090G.inc[PY090G.ind==1] <- generate.inc(x=PY090G,ind=PY090G.ind,classes=10)

o.indF <- PY090G.ind 


s.indF <- SynthDat$PY090G==1
gew=PB040/sum(PB040)

varsF=c("SEX","PB210","PE040","DOU","HHG")
MZ <- data.frame(SEX=PB150,PB210,PE040,DOU=DOU.p,HHG=HHG.p)

SynthDat$PY090G[s.indF] <- Synth.prep(yvarF=PY090G.inc,varsF,o.indF,s.indF,region,gew)

   # Hier resultieren noch viele fehlende Werte, deshalb wird die Zahl der Variablen, die zur Erzeugung genutzt werden reduziert.
   # Zunächst wird die Problemvariable PE040 weggelassen

s.indF <- is.na(SynthDat$PY090G) & SynthDat$RB210==2

varsF=c("SEX","PB210","DOU","HHG")
MZ <- data.frame(SEX=PB150,PB210,DOU=DOU.p,HHG=HHG.p)

SynthDat$PY090G[s.indF] <- Synth.prep(yvarF=PY090G.inc,varsF,o.indF,s.indF,region,gew)

   # Jetzt sind immer noch fehlende Werte enthalten, deshalb wird die Zahl der Erzeugungsvariablen weiter reduziert.

s.indF <- is.na(SynthDat$PY090G) & SynthDat$RB210==2

varsF=c("SEX","PB210","DOU")
MZ <- data.frame(SEX=PB150,PB210,DOU=DOU.p)

SynthDat$PY090G[s.indF] <- Synth.prep(yvarF=PY090G.inc,varsF,o.indF,s.indF,region,gew)

varsF=c("SEX","DOU")
MZ <- data.frame(SEX=PB150,DOU=DOU.p)
SynthDat$PY090G[s.indF] <- Synth.prep(yvarF=PY090G.inc,varsF,o.indF,s.indF,region,gew)

#################
# PY110G

   # Altersleistungen

SynthDat$PY110G <- NA


   # Altersleistungen werden nur bei Personen sind, die auch in Rente sind
SynthDat$PY110G[SynthDat$RB210==3] <- 1
SynthDat$PY110G[SynthDat$RB210 %in% c(1,2,4)] <- 0


varsF=c("ACL","SEX","PE040","FST","PB210")
MZ <- data.frame(ACL,SEX=PB150,PE040,FST=PB190,PB210)

s.indF <- SynthDat$PY110G==1 & SynthDat$RB210==3


PY110G.ind <- rep(F,length(PY110G))
PY110G.ind[RB210.p==3 & PY110G>1] <- T
PY110G.inc<-rep(0,length(PY110G.ind))

PY110G.inc[PY110G.ind==1] <- generate.inc(x=PY110G,ind=PY110G.ind,classes=10)

o.indF <- PY110G.ind 

SynthDat$PY110G[s.indF] <- Synth.prep(yvarF=PY110G.inc,varsF,o.indF,s.indF,region,gew)


s.indF <- is.na(SynthDat$PY110G) & SynthDat$RB210==3

   # Alle zu generierende Werte stehen noch auf fehlend
   # Die hohe Anzahl an NA's hängt hier mit der Variable PB210 zusammen, nachdem diese entfernt wurde resultieren keine NA's mehr
 
varsF=c("ACL","SEX","FST")
MZ <- data.frame(ACL,SEX=PB150,FST=PB190)
SynthDat$PY110G[s.indF] <- Synth.prep(yvarF=PY110G.inc,varsF,o.indF,s.indF,region,gew)


varsF=c("SEX","FST")
MZ <- data.frame(SEX=PB150,FST=PB190)
SynthDat$PY110G[s.indF] <- Synth.prep(yvarF=PY110G.inc,varsF,o.indF,s.indF,region,gew)

#################
# PY120G

   # Krankengeld

   # Das Krankengeld hat vermutlich einen sehr hohen Zusammenhang zu den PH Variablen.

   # 12917 Personen im Datensatz bekommen Krankengeld

SynthDat$PY120G <- NA

ind <- rep(T,length(PY120G))

o.indF <- ind
s.indF <- FILE=="P"

gew=PB040/sum(PB040)

PY120G.ind <- rep(F,length(PY120G))
PY120G.ind[RB210.p==3 & PY120G>1] <- T
PY120G.inc<-rep(0,length(PY120G.ind))
PY120G.inc[PY120G.ind==1] <- generate.inc(x=PY120G,ind=PY120G.ind,classes=10)


   # Es sind hier nicht mehr fehlende Werte erzeugt worden, wenn die Zahl der Klassen höher ist.
varsF=c("ACL","SEX","RB210","FST","Hcl","PH010")
MZ <- data.frame(ACL,SEX=PB150,RB210=RB210.p,PE040,FST=PB190,DOU=DOU.p,Hcl,PH010)
SynthDat$PY120G[s.indF] <- Synth.prep(yvarF=PY120G.inc,varsF,o.indF,s.indF,region,gew)


#################
# PY140G

   # ausbildungsbezogene Leistungen (PY140G)


   # Education allowances refers to grants, scholarships and other education help received by students.

SynthDat$PY140G <- NA

SynthDat$PY140G[SynthDat$PE010==1] <- 1
SynthDat$PY140G[SynthDat$PE010==2] <- 0


varsF=c("ACL","SEX","PB210","PY110G")
MZ <- data.frame(ACL,SEX=PB150,PB210,PY110G)

s.indF <- rep(F,length(SynthDat$PY140G))
s.indF[SynthDat$PY140G==1] <- T
s.indF[SynthDat$PE010==1] <- T


PY140G.ind <- rep(F,length(PY140G))
PY140G.ind[PE010==1 & PY140G>1] <- T

PY140G.inc<-rep(0,length(PY140G.ind))

PY140G.inc[PY140G.ind==T] <- generate.inc(x=PY140G,ind=PY140G.ind,classes=10)

o.indF <- PY140G.ind 

SynthDat$PY140G[s.indF] <- Synth.prep(yvarF=PY140G.inc,varsF,o.indF,s.indF,region,gew)


varsF=c("ACL","SEX","PY110G")
MZ <- data.frame(ACL,SEX=PB150,PY110G)
SynthDat$PY140G[s.indF] <- Synth.prep(yvarF=PY140G.inc,varsF,o.indF,s.indF,region,gew)


varsF=c("ACL","SEX")
MZ <- data.frame(ACL,SEX=PB150)
SynthDat$PY140G[s.indF] <- Synth.prep(yvarF=PY140G.inc,varsF,o.indF,s.indF,region,gew)


#################
# PY130G


   # im EU-SILC Datensatz gibt es die Ausprägung widowed (verwitwet)
   # FST ist hier eine wichtige Variable zu Erzeugung



SynthDat$PY130G <- NA

gew=PB040/sum(PB040)


   # to persons below standard retirement age

PY130G.ind <- rep(F,length(PY130G))
   # dieser Indikator zeigt an, für welche Personen eine solche Einkommenskomponente erzeugt wird

PY130G.ind[RX010.p<65 & PY130G>1 & RB210.p!=1 & RB210.p!=2] <- T
PY130G.inc<-rep(0,length(PY130G.ind))
PY130G.inc[PY130G.ind==1] <- generate.inc(x=PY130G,ind=PY130G.ind,classes=10)


   # Es sind hier nicht mehr fehlende Werte erzeugt worden, wenn die Zahl der Klassen höher ist.


o.indF <- PY130G.ind

   # Die Gruppe der Personen für die die Variable PY130G erzeugt werden muss ist relativ klein. Nur Personen die keine Altersleistungen erhalten, nicht arbeitslos oder erwerbstätig sind kommen für diese Zuwendungen in Frage.
s.indF <- FILE=="P" & SynthDat$RB210!=1 & SynthDat$RB210!=2 & SynthDat$AGE<65

   # Die Region sollte auch mit einbezogen werden
varsF=c("ACL","SEX","PH010","REG")
MZ <- data.frame(ACL,SEX=PB150,PH010,REG=region.p)
SynthDat$PY130G[s.indF] <- Synth.prep(yvarF=PY130G.inc,varsF,o.indF,s.indF,region,gew)

   # In der ersten Version mit der Variable PH010 wurden noch ausschließlich NA's erzeugt


#--------------------------------------------------------------------------#
# Die Einkommenskomponenten aus selbstständiger Erwerbsarbeit
#--------------------------------------------------------------------------#


#################
# PY050G

   # Bei der Einkommenskomponente PY050 sind auch negative Einkommen möglich, deshalb muss die Erzeugung ein wenig umgestellt werden.


   # 1. Schritt: Platz für die Variable in der Matrix SynthDat

SynthDat$PY050G <- NA

   # 2. Schritt: das Gewicht der Ausprägungen in der Ausgangspopulation bestimmen.

gew=PB040/sum(PB040)

   # 3. Schritt: Es wird ein Indikator für den Ausgangsdatensatz erzeugt, der anzeigt welche Personen im Datensatz als Vorlage zur Erzeugung dienen können. Wichtig ist, dass keine fehlenden Werte in den Originaldaten enthalten sind.


   # Hier können außerdem nur Werte berücksichtigt werden, die über null liegen, da sonst das Einteilen nach Quantile nicht mehr richtig klappt.

PY050G.ind <- rep(F,length(PY050G))
PY050G.ind[PY050G>1 & !is.na(PY050G)] <- T

   # Es ist wichtig, dass PY050G und o.indF unterschiedlich sind. Die Nullen dürfen bei der Kategorisierung des Originaleinkommens nicht berücksichtigt werden, bei der Erzeugung des synthetischen Einkommensklassen muss es aber dabei sein.

o.indF <-!is.na(PY050G)

   # 4. Schritt: Für den synthetischen Datensatz wird ebenfalls ein Indikator erstellt, der anzeigt für welche Personen im Datensatz die Variable erzeugt werden soll.

s.indF <- FILE=="P" 

   # 5. Schritt: Sofern es sich um Einkommenskomponenten handelt werden diese diskretisiert.

   # dieses kategorisierte Einkommen darf nicht die Länge von o.indF haben
PY050G.inc<-rep(0,length(PY050G))
PY050G.inc[PY050G.ind] <- generate.inc(x=PY050G,ind=PY050G.ind,classes=15)

   # 6. Schritt: Variablen die als Kovariaten zur Erzeugung dienen sollen werden festgelegt.
   # Dabei ist darauf zu achten, dass die Variablen einen möglichst hohen Zusammenhang zur Zielvariable aufweisen und das Editingregeln berücksichtigt werden.

varsF=c("ACL","SEX","PE040","REG","RB210","HHT")
MZ <- data.frame(ACL,SEX=PB150,PE040,REG=region.p,RB210=RB210.p,HHT=HY010.inc.p)


   # 7. Schritt: Aus den erzeugten Kreuztabellen mit relativen Häufigkeiten wird eine Ausprägung gezogen. Diese Kreuztabellen wurden mit Hilfe der Gewichte erzeugt. 

PY050G.ic <- PY050G.inc+1

   # Hier sollte info auch den Wert 1 enthalten
info <- Synth.prep(yvarF=PY050G.ic,varsF,o.indF,s.indF,region,gew)

SynthDat$PY050G[s.indF] <- as.numeric(as.character(info)) - 1

   # 8. Schritt: Es wird überprüft, wie viele fehlende Werte sich noch in den erzeugten Daten befinden. Die Generierung wird noch einmal mit weniger Kovariaten durchgeführt, falls sich noch fehlende Werte in den Daten befinden.

table(is.na(SynthDat$PY050G[s.indF]))

s.indF <- FILE=="P" & is.na(SynthDat$PY050G)
varsF=c("ACL","SEX","REG","RB210","HHT")
MZ <- data.frame(ACL,SEX=PB150,REG=region.p,RB210=RB210.p,HHT=HY010.inc.p)

info <- Synth.prep(yvarF=PY050G.ic,varsF,o.indF,s.indF,region,gew)

SynthDat$PY050G[s.indF] <- as.numeric(as.character(info)) - 1


#s.indF <- FILE=="P" & is.na(SynthDat$PY050G)
#varsF=c("ACL","SEX")
#MZ <- data.frame(ACL,SEX=PB150)
#SynthDat$PY050G[s.indF] <- Synth.prep(yvarF=PY050G.inc,varsF,o.indF,s.indF,region,gew)

#s.indF <- FILE=="P" & is.na(SynthDat$PY050G)
   
   # const ist eine Möglichkeit nur eine Variable als Hilfsinformation zur Erzeugung zu nutzen
#varsF=c("ACL","const")
#MZ <- data.frame(ACL,const=rep(1,length(PY050G)))
#SynthDat$PY050G[s.indF] <- Synth.prep(yvarF=PY050G.inc,varsF,o.indF,s.indF,region,gew)



#################
# PY070

   # Value of goods produced by own-consumption

   # hier stellt sich auch das Problem, das der Anteil der Nulleinkommen sehr hoch ist.

SynthDat$PY070G <- NA

   # hier können Verluste und Gewinne mit drin sein, deshalb entfällt die Bedingung PY070G>1
PY070G.ind <- rep(F,length(PY070G))
PY070G.ind[PY070G>1 & !is.na(PY070G)] <- T


o.indF <-!is.na(PY070G)

s.indF <- FILE=="P" 
PY070G.inc<-rep(0,length(PY070G.ind))

   # Achtung, hier stand vorher PY070G.inc[PY070G.ind==1] dadurch ist ein Fehler entstanden, weil 
   # PY070G.ind ein logical ist

PY070G.inc[PY070G.ind] <- generate.inc(x=PY070G,ind=PY070G.ind,classes=15)


   # Die produzierten Waren dürften sehr stark von dem Einkommen aus selbstständiger Erwerbsarbeit abhängen
varsF=c("ACL","SEX","PE040","REG","RB210","PY050G")
MZ <- data.frame(ACL,SEX=PB150,PE040,REG=region.p,RB210=RB210.p,PY050G=PY050G.inc)

PY070G.ic <- PY070G.inc + 1
info <- Synth.prep(yvarF=PY070G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$PY070G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- FILE=="P" & is.na(SynthDat$PY070G)
varsF=c("ACL","SEX","REG","RB210","PY050G")
MZ <- data.frame(ACL,SEX=PB150,REG=region.p,RB210=RB210.p,PY050G=PY050G.inc)

info <- Synth.prep(yvarF=PY070G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$PY070G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- FILE=="P" & is.na(SynthDat$PY070G)
varsF=c("SEX","RB210")
MZ <- data.frame(SEX=PB150,RB210=RB210.p)

info <- Synth.prep(yvarF=PY070G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$PY070G[s.indF] <- as.numeric(as.character(info)) - 1



#################
# PY100

   # Old-age benefits

   # hier stellt sich auch das Problem, das der Anteil der Nulleinkommen sehr hoch ist.

   # Die brutto Altersleistungen können nicht negativ sein. Im Originaldatensatz gibt es dennoch 10 Personen, die bei dieser Einkommenskomponente ein negatives Einkommen haben.

SynthDat$PY100G <- NA


PY100G.ind <- rep(F,length(PY100G))
PY100G.ind[PY100G>1 & !is.na(PY100G)] <- T

   # Enthalten sind auch die Gelder aus Vorrauhestandleistungen, im Originaldatensatz sind Personen mit 16 Jahren enthalten, die 

o.indF <-!is.na(PY100G)

s.indF <- FILE=="P" 
PY100G.inc<-rep(0,length(PY100G.ind))
PY100G.inc[PY100G.ind] <- generate.inc(x=PY100G,ind=PY100G.ind,classes=15)


varsF=c("ACL","SEX","PE040","REG","RB210","HHT")
MZ <- data.frame(ACL,SEX=PB150,PE040,REG=region.p,RB210=RB210.p,HHT=HY010.inc.p)

PY100G.ic <- PY100G.inc + 1
info <- Synth.prep(yvarF=PY100G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$PY100G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- FILE=="P" & is.na(SynthDat$PY100G)
varsF <- varsF[-which(varsF=="PE040")]
# MZ <- MZ[,-which(colnames(MZ)=="PE040")]

info <- Synth.prep(yvarF=PY100G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$PY100G[s.indF] <- as.numeric(as.character(info)) - 1

#--------------------------------------------------------------------------#
# Die Einkommenskomponenten verstetigen
#--------------------------------------------------------------------------#

   # INC ist das klassierte Einkommen in der EU_SILC Bevölkerung

INC <- PY010G.inc[PY010G.ind]  

   # folgendes stand sonst immer in der Funktion drin, geht aber verloren, deshalb hier noch mal durchgeführt

inc.cl <- sort(unique(INC))
min.inc <- tapply(PY010G[PY010G.ind],INC,min)
max.inc <- tapply(PY010G[PY010G.ind],INC,max)

attr(INC,"min.inc")<- min.inc
attr(INC,"max.inc")<- max.inc
attr(INC,"inc.cl")<- inc.cl


PY010D <- rep(NA,nrow(SynthDat))

   # hier sind leider noch fehlende Werte enthalten
PY010D[SynthDat$PY010G!=0] <- dis.inc.agg(x=PY010G,x.sbasic=SynthDat$PY010G,x.ind=PY010G.ind,x.inc= PY010G.inc,n.ex=40)




# Erzeugung müsste eigentlich mit dis.inc.agg gemacht werden

info <- quantile.korr(x="PY010G",acuteness=50)





   # PY020G muss Teil der synthetischen Matrix SynthDat sein
   # PY020G aus dem Originaldatensatz muss eingeladen sein
   # das Basiseinkommen muss auch diskretisiert vorliegen
   # felende Werte auf null setzen

   # Es sind jetzt keine NA's mehr im Datensatz. Die Einkommen an den Stellen an denen vorher ein NA stand sind jetzt auf Null gesetzt.



#################
# PY020

   # Bei dieser Einkommenskomponente ist ein Ausreßer drin, der die folgenden Werte total verzerrt

PY020D <- dis.inc.agg(x=PY020G,x.sbasic=SynthDat$PY020G,x.ind=PY020G.ind,x.inc= PY020G.inc,n.ex=2)

   # Quantilskorrektur


   # acuteness gibt die Feinheit der Quantile an, an die angepasst wird. Wenn die Feinheit gering gewählt wird, dann sind gerade in den hohen Einkommen gewisse Stufen zu sehen.
PY020D <- quantile.korr(x="PY020G",acuteness=20)


#################
# Weitere Variablen

SynthDat$PY050G[is.na(SynthDat$PY050G)] <- 0

# dis.inc.agg <- function(x,x.sbasic,x.ind,x.inc,n.ex){
      # x ist das Einkommen aus dem Originaldatensatz
      # x.sbasic ist das Einkommen aus dem synthetischen Datensatz in kategorialer Form      
      # x.ind sind Regeln für den Originaldatensatz, welche Werte dem synthetischen stetigen Einkommen zu Grunde liegen sollen
      # x.inc ist das kategoriale Einkommen des Originaldatensatzes
      # n.ex ist die Zahl der Extremwerte, die entfernt werden soll

PY050D <- dis.inc.agg(x=PY050G,x.sbasic=SynthDat$PY050G,x.ind=PY050G.ind,x.inc= PY050G.inc,n.ex=2)


   # die Verstetigung der Einkommen funktioniert nicht bei PY070D 
   # Checkliste wo Problem liegen könnte:

   # str(PY070G) vergleichen
   # fehlende Werte
   


SynthDat$PY070G[is.na(SynthDat$PY070G)] <- 0
PY070D <- dis.inc.agg(x=PY070G,x.sbasic=SynthDat$PY070G,x.ind=PY070G.ind,x.inc= PY070G.inc,n.ex=2)


SynthDat$PY090G[is.na(SynthDat$PY090G)] <- 0
PY090D <- dis.inc.agg("PY090G")

SynthDat$PY100G[is.na(SynthDat$PY100G)] <- 0
PY100D <- dis.inc.agg(x=PY100G,x.sbasic=SynthDat$PY100G,x.ind=PY100G.ind,x.inc= PY100G.inc,n.ex=2)



SynthDat$PY110G[is.na(SynthDat$PY110G)] <- 0
PY110D <- dis.inc.agg("PY110G")

SynthDat$PY120G[is.na(SynthDat$PY120G)] <- 0
PY120D[is.na(PY120D)] <- 0
PY120D <- dis.inc.agg("PY120G")

SynthDat$PY130G[is.na(SynthDat$PY130G)] <- 0
PY130D <- dis.inc.agg("PY130G")

SynthDat$PY140G[is.na(SynthDat$PY140G)] <- 0
PY140D <- dis.inc.agg(x=PY140G,x.sbasic=SynthDat$PY140G,x.ind=PY140G.ind,x.inc= PY140G.inc,n.ex=2)

PY140D <- quantile.korr(x="PY140G",acuteness=20)


#################################################################################################
# Erzeugten Variablen abspeichern
#################################################################################################

save(HHT, file="AML.HHT.RData")


PB210 <- SynthDat$PB210
save(PB210,file="/eingabe/Ameland/universe/results_v1/AML.PB210.RData")

PL080 <- SynthDat$PL080
save(PL080,file="/eingabe/Ameland/universe/results_v1/AML.PL080.RData")



PY010G <- SynthDat$PY010G
save(PY010G,file="/eingabe/Ameland/universe/results_v1/AML.PY010G.RData")

PY020G <- SynthDat$PY020G
save(PY020G,file="/eingabe/Ameland/universe/results_v1/AML.PY020G.RData")

PY050G <- SynthDat$PY050G
save(PY050G,file="/eingabe/Ameland/universe/results_v1/AML.PY050G.RData")

PY070G <- SynthDat$PY070G
save(PY070G,file="/eingabe/Ameland/universe/results_v1/AML.PY070G.RData")

PY090G <- SynthDat$PY090G
save(PY090G,file="/eingabe/Ameland/universe/results_v1/AML.PY090G.RData")

PY110G <- SynthDat$PY110G
save(PY110G,file="/eingabe/Ameland/universe/results_v1/AML.PY110G.RData")

PY120G <- SynthDat$PY120G
save(PY120G,file="/eingabe/Ameland/universe/results_v1/AML.PY120G.RData")

PY130G <- SynthDat$PY130G
save(PY130G,file="/eingabe/Ameland/universe/results_v1/AML.PY130G.RData")

PY140G <- SynthDat$PY140G
save(PY140G,file="/eingabe/Ameland/universe/results_v1/AML.PY140G.RData")

   # diskretisierte Einkommen abspeichern

save(PY010D,file="/eingabe/Ameland/universe/results_v1/AML.PY010D.RData")
save(PY020D,file="/eingabe/Ameland/universe/results_v1/AML.PY020D.RData")
save(PY050D,file="/eingabe/Ameland/universe/results_v1/AML.PY050D.RData")
save(PY070D,file="/eingabe/Ameland/universe/results_v1/AML.PY070D.RData")
save(PY090D,file="/eingabe/Ameland/universe/results_v1/AML.PY090D.RData")
save(PY100D,file="/eingabe/Ameland/universe/results_v1/AML.PY100D.RData")
save(PY110D,file="/eingabe/Ameland/universe/results_v1/AML.PY110D.RData")
save(PY120D,file="/eingabe/Ameland/universe/results_v1/AML.PY120D.RData")
save(PY130D,file="/eingabe/Ameland/universe/results_v1/AML.PY130D.RData")
save(PY140D,file="/eingabe/Ameland/universe/results_v1/AML.PY140D.RData")

PL080 <- SynthDat$PL080
save(PL080,file="/eingabe/Ameland/universe/results_v1/AML.PL080.RData")


   # Gesundheitsdaten abspeichern

PH010 <- SynthDat$PH010
save(PH010,file="/eingabe/Ameland/universe/results_v1/AML.PH010.RData")

PH020 <- SynthDat$PH020
save(PH020,file="/eingabe/Ameland/universe/results_v1/AML.PH020.RData")

PH030 <- SynthDat$PH030
save(PH030,file="/eingabe/Ameland/universe/results_v1/AML.PH030.RData")


PH040 <- SynthDat$PH040
save(PH040,file="/eingabe/Ameland/universe/results_v1/AML.PH040.RData")



#--------------------------------------------------------------------------#
# Haushaltseinkommen zusammenstellen
#--------------------------------------------------------------------------#

setwd(AML.path)

load("AML.PY010D.RData")
load("AML.PY020D.RData")

   # folgende Einkomemnskomponente fehlt im EU-SILC Datensatz
# load("AML.PY030D.RData")   # fehlt
load("AML.PY050D.RData")
load("AML.PY070D.RData")
load("AML.PY090D.RData")
load("AML.PY100D.RData")   
load("AML.PY110D.RData")
load("AML.PY120D.RData")
load("AML.PY130D.RData")
load("AML.PY140D.RData")


   # Es ist wichtig, das keine Komponente des persönlichen Einkommens fehlende Werte aufweist, da sonst die Summe aller Komponenten für diesen Haushalt auch zu einem fehlenden Wert wird.
   # hier wäre eigentlich auch PY030D enthalten
HY010.p <- PY010D + PY020D + PY050D + PY070D + PY090D + PY100D + PY110D + PY120D + PY130D + PY140D

tisna <- function(x)table(is.na(x))

tisna(PY010D)
tisna(PY020D)
tisna(PY030D)

save(HY010.p,file="/eingabe/Ameland/universe/results_v1/AML.HY010.p.RData")


    # Die Einkompenskomponenten werden auch für den Originaldatensatz aufsummiert

PY010G[is.na(PY010G)] <- 0
PY020G[is.na(PY020G)] <- 0
PY050G[is.na(PY050G)] <- 0
PY070G[is.na(PY070G)] <- 0
PY090G[is.na(PY090G)] <- 0
PY100G[is.na(PY100G)] <- 0
PY110G[is.na(PY110G)] <- 0
PY120G[is.na(PY120G)] <- 0
PY130G[is.na(PY130G)] <- 0
PY140G[is.na(PY140G)] <- 0


   # Die Summe dieser Einkommenskomponenten für den Originaldatensatz


 HY010.pO <- PY010G + PY020G          + PY050G + PY070G + PY090G + PY100G + PY110G + PY120G + PY130G + PY140G
#HY020new <- PY010G + PY020G + PY030G + PY050G + PY070G + PY090G + PY100G + PY110G + PY120G + PY130G + PY140G 

   # HY010 und HY020 sind prinzipiell gleich, Daten zu PY030G fehlen

HY010.phhO <- tapply(HY010.pO,PX030,sum,na.rm=T)
HY010.p.indO <- rep(T,length(HY010.phh))


   # Die Nulleinkommen müssen rauß genommen werden, da sonst keine vernünftigen Quantile gebildet werden können.
HY010.p.indO[is.na(HY010.phh)] <- F
HY010.p.indO[HY010.phh<1] <- F



HY010.p.incO <- rep(0,length(HY010.p.indO))

HY010.p.incO[HY010.p.indO] <- generate.inc(x=HY010.phh,ind=HY010.p.indO,classes=10)

#--------------------------------------------------------------------------#
# Haushaltseinkommenskomponenten generieren
#--------------------------------------------------------------------------#


#################
# 

   # Achtung ind wird weiter oben noch in weiteren Zusammenhängen verwendet.
ind=unlist(HID.u2)

   # Hier sollte das aufsummierte Haushaltseinkommen berücksichtigt werden

SynthDat <- as.data.frame(ind)

   # Die Region ist ebenfalls für alle Mitglieder eines Haushalts gleich

SynthDat$REG <- REG[ind]

SynthDat$HHG <- HHG[ind]

SynthDat$DOU <- DOU[ind]



   # Die Anzahl der Kinder in den Haushalten des Originaldatensatzes muss festgestellt werden.

CPHH <- tapply(RX010<16,RX030,sum)

CPHH[CPHH>4] <- 4

   # Die oberste Klasse wird in eine offene Klasse umgewandelt

CPHHs <- tapply(AGE<16,HID,sum)

CPHHs[CPHHs>4] <- 4

   # hoffentlich ist die Reihenfolge die richtige.
SynthDat.HH$CPHH <- CPHHs

   # Die Anzahl der verheirateten Personen pro Haushalt spielt demgegnüber keine große Rolle 
   # MPHH <- tapply(PB190==1,PX030,sum)

   # Das aufsummierte Einkommen der persönlichen Einkommenskomponenten pro Haushalt wird kategorisiert um dann
   # der Erzeugung weiterer Haushaltseinkommen zu dienen. 

HY010.phh <- tapply(HY010.p,HID,sum,na.rm=T)


HY010.p.ind <- rep(T,length(HY010.phh))

   # ist das aufsummierte Einkommen auch an manchen Stellen null?
   # Auch das aus den persönlichen Nettoeinkommen aufsummierte Haushaltseinkommen ist für manche Haushalte gleich null.


o.indF <- HY010.p.ind
s.indF <- HY010.p.ind
HY010.p.inc <- generate.inc(x=HY010.phh,ind=HY010.p.ind,classes=10)

SynthDat$HY010.p <- HY010.p.inc

   # Die gleichen Einkommensklassen müssen nun noch für den Originaldatrensatz erstellt werden.


HY010.p.incO <- generate.inc(x=HY010.phh,ind=HY010.p.ind,classes=10)

SynthDat <- SynthDat.HH 

#########
# HY030G 

   # unterstellte Miete

SynthDat$HY030G <- NA

HY030G.ind <- rep(F,length(HY030G))
HY030G.ind[HY030G>0 & !is.na(HY030G)] <- T

o.indF <-!is.na(HY030G)

s.indF <- rep(T,nrow(SynthDat))
HY030G.inc<-rep(0,length(HY030G.ind))
HY030G.inc[HY030G.ind] <- generate.inc(x=HY030G,ind=HY030G.ind,classes=15)


varsF=c("REG","HHG","DOU","HY010.p")
MZ <- data.frame(REG=region.h,HHG=HX040,DOU=DB100,HY010.p=HY010.p.incO)

gew <- DB090/sum(DB090)


HY030G.ic <- HY030G.inc + 1
info <- Synth.prep(yvarF=HY030G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY030G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY030G)

varsF=c("REG","HHG","DOU")
info <- Synth.prep(yvarF=HY030G,varsF,o.indF,s.indF,region,gew)
SynthDat$HY030G[s.indF] <- as.numeric(as.character(info)) - 1



#########
# HY040G

   # Einkommen aus Verpachtung und Vermietung (HY030G)


SynthDat$HY040G <- NA

HY040G.ind <- rep(F,length(HY040G))
HY040G.ind[HY040G>0 & !is.na(HY040G)] <- T

o.indF <-!is.na(HY040G)

s.indF <- rep(T,nrow(SynthDat))
HY040G.inc<-rep(0,length(HY040G.ind))
HY040G.inc[HY040G.ind] <- generate.inc(x=HY040G,ind=HY040G.ind,classes=15)


varsF=c("REG","HHG","DOU")
MZ$CPHH <- CPHH 

gew <- DB090/sum(DB090)


HY040G.ic <- HY040G.inc + 1
info <- Synth.prep(yvarF=HY040G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY040G[s.indF] <- as.numeric(as.character(info)) - 1


#########
# HY050G


   # Familienleistungen/ Kindergeld

   # Bei den Familienleistungen spielt die Anzahl der Kinder im Haushalt eine Rolle.


SynthDat$HY050G <- NA

HY050G.ind <- rep(F,length(HY050G))
HY050G.ind[HY050G>0 & !is.na(HY050G)] <- T

o.indF <-!is.na(HY050G)

s.indF <- rep(T,nrow(SynthDat))
HY050G.inc<-rep(0,length(HY050G.ind))
HY050G.inc[HY050G.ind] <- generate.inc(x=HY050G,ind=HY050G.ind,classes=15)


varsF=c("REG","HHG","DOU","CPHH")
MZ <- data.frame(REG=region.h,HHG=HX050,DOU=DB100,CPHH)

gew <- DB090/sum(DB090)


HY050G.ic <- HY050G.inc + 1
info <- Synth.prep(yvarF=HY050G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY050G[s.indF] <- as.numeric(as.character(info)) - 1

   # Wichtig bei dieser Einkommenskomponente ist, das Familienleistungen und Kindergeld zusammengefasst sind. Haushalte können also auch ein Einkommen dieser Form haben, wenn keine Kinder im Haushalt leben. Umgekehrt müssen keine Leistungen gezahlt werden auch wenn sich Kinder im Haushalt befinden.


#########
# HY060G

   # Sonstige Leistungen gegen soziale Ausgrenzung


   # Familienleistungen/ Kindergeld

   # Bei den Familienleistungen spielt die Anzahl der Kinder im Haushalt eine Rolle.


SynthDat$HY060G <- NA

HY060G.ind <- rep(F,length(HY060G))
HY060G.ind[HY060G>0 & !is.na(HY060G)] <- T

o.indF <-!is.na(HY060G)

s.indF <- rep(T,nrow(SynthDat))
HY060G.inc<-rep(0,length(HY060G.ind))
HY060G.inc[HY060G.ind] <- generate.inc(x=HY060G,ind=HY060G.ind,classes=15)


varsF=c("REG","HHG","DOU","CPHH","HY010.p")
MZ <- data.frame(REG=region.h,HHG=HX050,DOU=DB100,CPHH,HY010.p=HY010.p.incO)

gew <- DB090/sum(DB090)


HY060G.ic <- HY060G.inc + 1
info <- Synth.prep(yvarF=HY060G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY060G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY060G[s.indF])

   # Eigentlich müsste es reichen diese Variablenliste zu verkürzen und die MZ unangetastet zu lassen

varsF=c("REG","HHG","DOU")
info <- Synth.prep(yvarF=HY060G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY060G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- is.na(SynthDat$HY060G[s.indF])
varsF=c("HHG","DOU")
info <- Synth.prep(yvarF=HY060G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY060G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY060G[s.indF])
varsF=c("REG","DOU")
info <- Synth.prep(yvarF=HY060G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY060G[s.indF] <- as.numeric(as.character(info)) - 1



# SynthDat$HY060G <- gen.vars(x="HY060",varsF=c("REG","HHG","DOU","HY010.p","CPHH"),MZ)


#########
# HY070G

   # Wohnungsbeihilfen


SynthDat$HY070G <- NA

HY070G.ind <- rep(F,length(HY070G))
HY070G.ind[HY070G>0 & !is.na(HY070G)] <- T

o.indF <-!is.na(HY070G)

s.indF <- rep(T,nrow(SynthDat))
HY070G.inc<-rep(0,length(HY070G.ind))
HY070G.inc[HY070G.ind] <- generate.inc(x=HY070G,ind=HY070G.ind,classes=15)


varsF=c("REG","HHG","DOU","CPHH","HY010.p")
MZ <- data.frame(REG=region.h,HHG=HX050,DOU=DB100,CPHH,HY010.p=HY010.p.incO)

gew <- DB090/sum(DB090)


HY070G.ic <- HY070G.inc + 1
info <- Synth.prep(yvarF=HY070G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY070G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY070G[s.indF])
varsF=c("REG","DOU","HY010.p")
info <- Synth.prep(yvarF=HY070G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY070G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY070G[s.indF])
varsF=c("REG","DOU")
info <- Synth.prep(yvarF=HY070G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY070G[s.indF] <- as.numeric(as.character(info)) - 1



#########
# HY080G

   # Regelmäßig empfangene Geldtransfers zwischen privaten Haushalten



SynthDat$HY080G <- NA

HY080G.ind <- rep(F,length(HY080G))
HY080G.ind[HY080G>0 & !is.na(HY080G)] <- T

o.indF <-!is.na(HY080G)

s.indF <- rep(T,nrow(SynthDat))
HY080G.inc<-rep(0,length(HY080G.ind))
HY080G.inc[HY080G.ind] <- generate.inc(x=HY080G,ind=HY080G.ind,classes=15)


varsF=c("REG","HHG","DOU","HY010.p")


HY080G.ic <- HY080G.inc + 1
info <- Synth.prep(yvarF=HY080G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY080G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- is.na(SynthDat$HY080G[s.indF])
varsF=c("REG","DOU")
info <- Synth.prep(yvarF=HY080G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY080G[s.indF] <- as.numeric(as.character(info)) - 1


#########
# HY090G

   # Brutto Zinsen, Brutto-Dividenden, Brutto-Gewinne aus Kapitalanlagen in Unternehmen ohne eigene Rechtspersönlichkeit



SynthDat$HY090G <- NA

HY090G.ind <- rep(F,length(HY090G))
HY090G.ind[HY090G>0 & !is.na(HY090G)] <- T

o.indF <-!is.na(HY090G)

s.indF <- rep(T,nrow(SynthDat))
HY090G.inc<-rep(0,length(HY090G.ind))
HY090G.inc[HY090G.ind] <- generate.inc(x=HY090G,ind=HY090G.ind,classes=15)

   # Hier wäre eventuell die Information von Interesse, wieviele Personen im Haushalt selbstständig sind
varsF=c("REG","HHG","DOU","HY010.p")


HY090G.ic <- HY090G.inc + 1
info <- Synth.prep(yvarF=HY090G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY090G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- is.na(SynthDat$HY090G[s.indF])
varsF=c("REG","DOU")
info <- Synth.prep(yvarF=HY090G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY090G[s.indF] <- as.numeric(as.character(info)) - 1




#########
# HY110G

   # Von Personen unter 16 Jahren erhaltenes Brutto-Einkommen

Pu16 <- tapply(AGE.r<16,RX030,sum)

   # Diese Variable muss als Bedingung in die Erzeugungmit einfließen.

   # Dann muss aber auch eine Variable mit einer nach oben offenen Klasse eingeführt werden

Pu16clo <- Pu16
Pu16clo[Pu16>5] <- 5


Pu16s <- tapply(AGE<16,HID,sum)

Pu16cl <- Pu16s
Pu16cl[Pu16s>5] <- 5

SynthDat$Pu16cl <- Pu16cl

MZ$Pu16cl <- Pu16clo


   # Eigentliche Erzeugung


SynthDat$HY110G <- NA

HY110G.ind <- rep(F,length(HY110G))
HY110G.ind[HY110G>0 & !is.na(HY110G) & Pu16clo>0] <- T

   # hier dürfen wirklich nur die Werte rauß genommen werden, die fehlend sind
o.indF <-!is.na(HY110G)

s.indF <- rep(T,nrow(SynthDat))
s.indF[Pu16cl>0] <- F

HY110G.inc<-rep(0,length(HY110G.ind))
HY110G.inc[HY110G.ind] <- generate.inc(x=HY110G,ind=HY110G.ind,classes=15)


   # Es gibt nur sehr wenige Haushalte, in denen diese Art von Einkommen vorkommt.

varsF=c("REG","Pu16cl","HHG","DOU","HY010.p")

HY110G.ic <- HY110G.inc + 1
info <- Synth.prep(yvarF=HY110G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY110G[s.indF] <- as.numeric(as.character(info)) - 1

   # Die Zahl der fehlenden Werte hängt hier nicht von der Zahl der Klassen ab, auch nicht von der Zahl der verwendeten Hilfsvariablen. 

s.indF <- s.indF & is.na(s.indF)
varsF=c("REG","Pu16cl","HHG","DOU")
info <- Synth.prep(yvarF=HY110G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY110G[s.indF] <- as.numeric(as.character(info)) - 1



#########
# HY100G

   # Zinsen für Hypothekarkredite


SynthDat$HY100G <- NA

HY100G.ind <- rep(F,length(HY100G))
HY100G.ind[HY100G>0 & !is.na(HY100G)] <- T

o.indF <-!is.na(HY100G)

s.indF <- rep(T,nrow(SynthDat))
HY100G.inc<-rep(0,length(HY100G.ind))
HY100G.inc[HY100G.ind] <- generate.inc(x=HY100G,ind=HY100G.ind,classes=15)

   # Hier wäre eventuell die Information von Interesse, wieviele Personen im Haushalt selbstständig sind
varsF=c("REG","HHG","DOU","HY010.p")


HY100G.ic <- HY100G.inc + 1
info <- Synth.prep(yvarF=HY100G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY100G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- is.na(SynthDat$HY100G[s.indF])
varsF=c("REG","DOU")
info <- Synth.prep(yvarF=HY100G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY100G[s.indF] <- as.numeric(as.character(info)) - 1


MZ$const <- rep(1,nrow(MZ))
SynthDat$const <- rep(1,nrow(SynthDat))

s.indF <- is.na(SynthDat$HY100G[s.indF])
varsF=c("REG","const")
info <- Synth.prep(yvarF=HY100G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY100G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY100G[s.indF])
varsF=c("DOU","const")
info <- Synth.prep(yvarF=HY100G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY100G[s.indF] <- as.numeric(as.character(info)) - 1


#########
# HY120G

   # regelmäßige Vermögenssteuern

SynthDat$HY120G <- NA

HY120G.ind <- rep(F,length(HY120G))
HY120G.ind[HY120G>0 & !is.na(HY120G)] <- T

o.indF <-!is.na(HY120G)

s.indF <- rep(T,nrow(SynthDat))
HY120G.inc<-rep(0,length(HY120G.ind))
HY120G.inc[HY120G.ind] <- generate.inc(x=HY120G,ind=HY120G.ind,classes=15)


varsF=c("REG","HHG","DOU","HY090G")

MZ$HY090G <- HY090G.inc


HY120G.ic <- HY120G.inc + 1
info <- Synth.prep(yvarF=HY120G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY120G[s.indF] <- as.numeric(as.character(info)) - 1

   # Die Variable HY090G ist zwar wichtig in diesem Zusammenhang führt aber auch zu sehr vielen missings

s.indF <- is.na(SynthDat$HY120G[s.indF])
varsF=c("DOU","REG")
info <- Synth.prep(yvarF=HY120G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY120G[s.indF] <- as.numeric(as.character(info)) - 1



#########
# HY130G

   # regelmäßig geleistete Geldtransfers zwischen Haushalten

SynthDat$HY130G <- NA

HY130G.ind <- rep(F,length(HY130G))
HY130G.ind[HY130G>0 & !is.na(HY130G)] <- T

o.indF <-!is.na(HY130G)

s.indF <- rep(T,nrow(SynthDat))
HY130G.inc<-rep(0,length(HY130G.ind))
HY130G.inc[HY130G.ind] <- generate.inc(x=HY130G,ind=HY130G.ind,classes=15)


varsF=c("REG","HHG","DOU","HY010.p")


HY130G.ic <- HY130G.inc + 1
info <- Synth.prep(yvarF=HY130G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY130G[s.indF] <- as.numeric(as.character(info)) - 1

   # Die Variable HY090G ist zwar wichtig in diesem Zusammenhang führt aber auch zu sehr vielen missings

s.indF <- is.na(SynthDat$HY130G[s.indF])
varsF=c("DOU","REG","HHG")
info <- Synth.prep(yvarF=HY130G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY130G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- is.na(SynthDat$HY130G[s.indF])
varsF=c("DOU","REG")
info <- Synth.prep(yvarF=HY130G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY130G[s.indF] <- as.numeric(as.character(info)) - 1


#########
# HY140G


   # Einkommenssteuern und Sozialbeiträge

SynthDat$HY140G <- NA

HY140G.ind <- rep(F,length(HY140G))
HY140G.ind[HY140G>0 & !is.na(HY140G)] <- T

o.indF <-!is.na(HY140G)

s.indF <- rep(T,nrow(SynthDat))
HY140G.inc<-rep(0,length(HY140G.ind))
HY140G.inc[HY140G.ind] <- generate.inc(x=HY140G,ind=HY140G.ind,classes=15)


varsF=c("REG","HHG","DOU","HY010.p")


HY140G.ic <- HY140G.inc + 1
info <- Synth.prep(yvarF=HY140G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY140G[s.indF] <- as.numeric(as.character(info)) - 1

   # Die Variable HY090G ist zwar wichtig in diesem Zusammenhang führt aber auch zu sehr vielen missings

s.indF <- is.na(SynthDat$HY140G[s.indF])
varsF=c("HY010.p","REG")
info <- Synth.prep(yvarF=HY140G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY140G[s.indF] <- as.numeric(as.character(info)) - 1

s.indF <- is.na(SynthDat$HY140G[s.indF])
varsF=c("HY010.p","const")
info <- Synth.prep(yvarF=HY140G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY140G[s.indF] <- as.numeric(as.character(info)) - 1


s.indF <- is.na(SynthDat$HY140G[s.indF])
varsF=c("DOU","REG")
info <- Synth.prep(yvarF=HY140G.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY140G[s.indF] <- as.numeric(as.character(info)) - 1


#########
# HY025

   # Inflationsfaktor für Antwortausfälle im Haushalt

SynthDat$HY025 <- NA

HY025.ind <- rep(F,length(HY025))
HY025.ind[HY025>1 & !is.na(HY025)] <- T

o.indF <-!is.na(HY025)

s.indF <- rep(T,nrow(SynthDat))
HY025.inc<-rep(0,length(HY025.ind))
HY025.inc[HY025.ind] <- generate.inc(x=HY025,ind=HY025.ind,classes=10)


varsF=c("REG","HHG","DOU")


HY025.ic <- HY025.inc + 1
info <- Synth.prep(yvarF=HY025.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY025[s.indF] <- as.numeric(as.character(info))


s.indF <- is.na(SynthDat$HY025[s.indF])
varsF=c("DOU","REG")
info <- Synth.prep(yvarF=HY025.ic,varsF,o.indF,s.indF,region,gew)
SynthDat$HY025[s.indF] <- as.numeric(as.character(info)) 



#--------------------------------------------------------------------------#
# Haushaltseinkommenskomponenten verstetigen
#--------------------------------------------------------------------------#

SynthDat$HY060G[is.na(SynthDat$HY060G)] <- 0


   # bei der Variable HY030G ist etwas schief gelaufen
   # bei HY040G nur Nuller


HY030D <- dis.inc.agg(x=HY030G,x.sbasic=SynthDat$HY030G,x.ind=HY030G.ind,x.inc= HY030G.inc,n.ex=2)

   # hier sind noch fehlende Werte drin, diese werden jetzt auf 0 gesetzt

HY030D[is.na(HY030D)] <- 0


HY030D <- quantile.korr(x="HY030G",acuteness=10)

HY040D <- dis.inc.agg(x=HY040G,x.sbasic=SynthDat$HY040G,x.ind=HY040G.ind,x.inc= HY040G.inc,n.ex=2)
HY050D <- dis.inc.agg(x=HY050G,x.sbasic=SynthDat$HY050G,x.ind=HY050G.ind,x.inc= HY050G.inc,n.ex=2)
HY060D <- dis.inc.agg(x=HY060G,x.sbasic=SynthDat$HY060G,x.ind=HY060G.ind,x.inc= HY060G.inc,n.ex=2)
HY070D <- dis.inc.agg(x=HY070G,x.sbasic=SynthDat$HY070G,x.ind=HY070G.ind,x.inc= HY070G.inc,n.ex=2)

HY080D <- dis.inc.agg(x=HY080G,x.sbasic=SynthDat$HY080G,x.ind=HY080G.ind,x.inc= HY080G.inc,n.ex=2)



HY090D <- dis.inc.agg(x=HY090G,x.sbasic=SynthDat$HY090G,x.ind=HY090G.ind,x.inc= HY090G.inc,n.ex=2)


SynthDat$HY100G[is.na(SynthDat$HY100G)] <- 0

   # hier muss darauf geachtet werden, dass auch die richtige Variable eingeladen ist, bspw gibt es die Variable HY100G zweimal, einmal für den Originaldatensatz und einmal für den synthetischen Datensatz
   # wenn die Variable kategorial ist, kommt eine Fehlermeldung heraus.

HY100D <- dis.inc.agg(x=HY100G,x.sbasic=SynthDat$HY100G,x.ind=HY100G.ind,x.inc= HY100G.inc,n.ex=2)

HY110D <- dis.inc.agg(x=HY110G,x.sbasic=SynthDat$HY110G,x.ind=HY110G.ind,x.inc= HY110G.inc,n.ex=2)

   # Die Zahl der Einkommen die hinten als Ausreißer ausgeschlossen wird könnte erhöht werden.
HY120D <- dis.inc.agg(x=HY120G,x.sbasic=SynthDat$HY120G,x.ind=HY120G.ind,x.inc= HY120G.inc,n.ex=2)
HY130D <- dis.inc.agg(x=HY130G,x.sbasic=SynthDat$HY130G,x.ind=HY130G.ind,x.inc= HY130G.inc,n.ex=2)

   # Der höchste Wert der Variable HY040G liegt bei 2319821.5 Euro, dieser Wert wird bei der Erzeugung ausgespart.
# HY140G[which.max(HY140G)] <- NA
   # auch der zweithöchste Wert sollte eliminiert werden
# HY140G[which.max(HY140G)] <- NA

   # die höchsten 30 Werte werden eliminiert
HY140G[order(HY140G,decreasing=T)][1:30] <- NA


   # Achtung, wenn das so gemacht werden soll, dann muss der Indikator vorne auch noch verändert werden
HY140D <- dis.inc.agg(x=HY140G,x.sbasic=SynthDat$HY140G,x.ind=HY140G.ind,x.inc= HY140G.inc,n.ex=2)

HY025D <- dis.inc.agg(x=HY025,x.sbasic=SynthDat$HY025,x.ind=HY025.ind,x.inc= HY025.inc,n.ex=2)

#--------------------------------------------------------------------------#
# Haushaltseinkommenskomponenten mit Einkommen aus persönlichen Komponenten addieren
#--------------------------------------------------------------------------#


HY010new <- HY010.phh + HY030D + HY040D + HY050D + HY060D + HY070D + HY080D + HY090D + HY110D - HY100D


   # dieses Einkommen unterscheidet sich nur geringfügig in den HH komponenten
HY020new <- HY010.phh + HY030D + HY040D + HY050D + HY060D + HY070D + HY080D + HY090D + HY110D - HY100D - HY120D - HY130D - HY140D


   # das tatsächliche HHeinkommen, wenn man die Komponenten aufsummiert
HY020tat <- HY010.phhO + HY030G + HY040G + HY050G + HY060G + HY070G + HY080G + HY090G + HY110G - HY100G - HY120G - HY130G - HY140G

save(HY010new,file="/eingabe/Ameland/universe/results_v1/AML.HY010new.RData")
save(HY020new,file="/eingabe/Ameland/universe/results_v1/AML.HY020new.RData")


#--------------------------------------------------------------------------#
# Haushaltseinkommenskomponenten abspeichern
#--------------------------------------------------------------------------#

HY030G <- SynthDat$HY030G
save(HY030G,file="/eingabe/Ameland/universe/results_v1/AML.HY030G.RData")


HY040G <- SynthDat$HY040G
save(HY040G,file="/eingabe/Ameland/universe/results_v1/AML.HY040G.RData")

HY060G <- SynthDat$HY060G
save(HY060G,file="/eingabe/Ameland/universe/results_v1/AML.HY060G.RData")

HY070G <- SynthDat$HY070G
save(HY070G,file="/eingabe/Ameland/universe/results_v1/AML.HY070G.RData")

HY080G <- SynthDat$HY080G
save(HY080G,file="/eingabe/Ameland/universe/results_v1/AML.HY080G.RData")

HY090G <- SynthDat$HY090G
save(HY090G,file="/eingabe/Ameland/universe/results_v1/AML.HY090G.RData")

HY110G <- SynthDat$HY110G
save(HY110G,file="/eingabe/Ameland/universe/results_v1/AML.HY110G.RData")


HY100G <- SynthDat$HY100G
save(HY100G,file="/eingabe/Ameland/universe/results_v1/AML.HY100G.RData")

HY120G <- SynthDat$HY120G
save(HY120G,file="/eingabe/Ameland/universe/results_v1/AML.HY120G.RData")

HY130G <- SynthDat$HY130G
save(HY130G,file="/eingabe/Ameland/universe/results_v1/AML.HY130G.RData")

HY140G <- SynthDat$HY140G
save(HY140G,file="/eingabe/Ameland/universe/results_v1/AML.HY140G.RData")

HY025 <- SynthDat$HY025
save(HY025,file="/eingabe/Ameland/universe/results_v1/AML.HY025.RData")

save(HY030D,file="/eingabe/Ameland/universe/results_v1/AML.HY030D.RData")
save(HY040D,file="/eingabe/Ameland/universe/results_v1/AML.HY040D.RData")
save(HY050D,file="/eingabe/Ameland/universe/results_v1/AML.HY050D.RData")
save(HY060D,file="/eingabe/Ameland/universe/results_v1/AML.HY060D.RData")
save(HY070D,file="/eingabe/Ameland/universe/results_v1/AML.HY070D.RData")
save(HY080D,file="/eingabe/Ameland/universe/results_v1/AML.HY080D.RData")
save(HY090D,file="/eingabe/Ameland/universe/results_v1/AML.HY090D.RData")
save(HY100D,file="/eingabe/Ameland/universe/results_v1/AML.HY100D.RData")
save(HY110D,file="/eingabe/Ameland/universe/results_v1/AML.HY110D.RData")
save(HY120D,file="/eingabe/Ameland/universe/results_v1/AML.HY120D.RData")
save(HY130D,file="/eingabe/Ameland/universe/results_v1/AML.HY130D.RData")
save(HY140D,file="/eingabe/Ameland/universe/results_v1/AML.HY140D.RData")


#--------------------------------------------------------------------------#
# SynthDat.HH vervollständigen
#--------------------------------------------------------------------------#

   # Die Region sollte in diesem File enthalten sein

SynthDat$CIT <- CIT[SynthDat$ind]



#--------------------------------------------------------------------------#
# SynthDat abspeichern
#--------------------------------------------------------------------------#

setwd(AML.path)

save(SynthDat, file="/eingabe/Ameland/universe/results_v1/SynthDat.RData")

   # Das gleiche auch für die Haushalts Datei
save(SynthDat, file="/eingabe/Ameland/universe/results_v1/SynthDat.HH.RData")

#--------------------------------------------------------------------------#
# Haushalte korrigieren
#--------------------------------------------------------------------------#

   # Es gibt sehr viele Haushalte, in denen die älteste Person jünger als 16 ist.

min.AGE <- tapply(AGE,HID,max)

ind <- which(min.AGE<18)

ind2 <- names(min.AGE)[ind]

info <- match(HID,ind2)

   # nur die Haushalte die es betrifft
AGE.i <- AGE[!is.na(info)]
HID.i <- HID[!is.na(info)]

   # wie groß sind die Problemhaushalte

   # Das höchste Alter ist 80
   # es werden gleichverteilte Daten gezogen

tab.AGE <- table(AGE[AGE>18])



   # hier soll nur eine Person über die Altersgrenze gehoben werden

indi <- tapply(1:length(HID.i),HID.i,function(x)x[1])

AGE.i[indi] <- as.integer(sample(names(tab.AGE),length(AGE.i[indi]),prob=prop.table(tab.AGE),replace=T))

AGE[!is.na(info)] <- AGE.i

save(AGE,file="/eingabe/Ameland/universe/results_v1/AML.AGE.RData")

#--------------------------------------------------------------------------#
# Einkommen korrigieren
#--------------------------------------------------------------------------#

   # es gibt sehr viele Einkommen, die sehr stark ins negative gehen

ind <- which.min(HY020new)

HY020G <- HY020

HY020D <- HY020new

# HY020D <- quantile.korr(x="HY020G",acuteness=20)

HY020D <- quantile.korr2(x="HY020G")

HY020D[is.na(HY020D)] <- 0
HY020D[HY020D<0] <- .05*HY020D[HY020D<0]

# HY020D[HY020D>60000] <- .08*HY020D[HY020D>60000]

HY020Dn <- HY020D

save(HY020Dn,file="AML.HY020Dn.RData")

   # Das Gesamteinkommen ist nun korrigiert und muss auf die Huashalte und Personen zurück gespielt werden

   # Differenz um die es geht:

   # positive Einkommenskomponenten
# Einkommen <- cbind(HY010.phh,HY030D,HY040D,HY050D,HY060D,HY070D,HY080D,HY090D,HY110D)

# Eink.prop <- apply(Einkommen,1,function(x)x/sum(x))

Diff <- HY020D - HY020Dn

diff_ <- function(x,Differenz,GE){
      # Differenz
      # GE ist das Gesamteinkommen
   x[is.na(x)]<-0
   x <- x - x/GE * Differenz
   return(x)
}

HY110Dn <- diff_(x=HY110D,Differenz=Diff,GE=HY020D)
HY090Dn <- diff_(x=HY090D,Differenz=Diff,GE=HY020D)
HY080Dn <- diff_(x=HY080D,Differenz=Diff,GE=HY020D)
HY070Dn <- diff_(x=HY070D,Differenz=Diff,GE=HY020D)
HY060Dn <- diff_(x=HY060D,Differenz=Diff,GE=HY020D)
HY050Dn <- diff_(x=HY050D,Differenz=Diff,GE=HY020D)
HY040Dn <- diff_(x=HY040D,Differenz=Diff,GE=HY020D)
HY030Dn <- diff_(x=HY030D,Differenz=Diff,GE=HY020D)
HY010.phhn<- diff_(x=HY010.phh,Differenz=Diff,GE=HY020D)

   # jetzt ist die Frage, einen wie großen Anteil jede Person am HHeinkommen hat.

Ant.p <- tapply(HY010.p,HID,function(x)x/sum(x,na.rm=T))

   # welche Person betrifft das jeweils
   # Hier könnte man noch schauen, wie die Einkommen in den einzelnen Haushalten verteilt sind

Ant.p[is.na(Ant.p)] <- 0

HY010.phhn.list <- as.list(HY010.phhn)

for (i in 1:length(Ant.p)){
   Ant.p[[i]] <- Ant.p[[i]]*HY010.phhn.list[[i]]
   # cat(i, "\n")
}

   # Diese Differenz die für jede einzelne Person noch begezogen werden muss, wird jetzt als Einzeldatensatz dargestellt.

ind.p <- tapply(1:length(HID),HID,function(x)x)

   # Die beiden listen werden in vektoren verwandelt
ind.p.u <- unlist(ind.p)

Ant.p.u <- unlist(Ant.p)

   # und miteinander gematcht
ind.HH <- match(1:length(HID),ind.p.u)

Ant.pp <- Ant.p.u[ind.HH]

   # jetzt wird wieder obige Funktion verwendet

Diff.p <- HY010.p - Ant.pp
PY140Dn <- diff_(x=PY140D,Differenz=Diff.p,GE=HY010.p)
PY130Dn <- diff_(x=PY130D,Differenz=Diff.p,GE=HY010.p)
PY120Dn <- diff_(x=PY120D,Differenz=Diff.p,GE=HY010.p)
PY110Dn <- diff_(x=PY110D,Differenz=Diff.p,GE=HY010.p)
PY100Dn <- diff_(x=PY100D,Differenz=Diff.p,GE=HY010.p)
PY090Dn <- diff_(x=PY090D,Differenz=Diff.p,GE=HY010.p)
PY070Dn <- diff_(x=PY070D,Differenz=Diff.p,GE=HY010.p)
PY050Dn <- diff_(x=PY050D,Differenz=Diff.p,GE=HY010.p)
PY020Dn <- diff_(x=PY020D,Differenz=Diff.p,GE=HY010.p)
PY010Dn <- diff_(x=PY010D,Differenz=Diff.p,GE=HY010.p)

   # jetzt wird das neu errechnete Einkommen wieder zusammen gezählt und geschaut ob das gewünschte Ergebnis resultiert.

HY010.pn <- PY010Dn + PY020Dn + PY050Dn + PY070Dn + PY090Dn + PY100Dn + PY110Dn + PY120Dn + PY130Dn + PY140Dn

HY010.phhn <- tapply(HY010.pn,HID,sum,na.rm=T)

   # Der ANteil wurde nicht auf die negativen Komponenten umgerechnet

HY020new2 <- HY010.phhn + HY030Dn + HY040Dn + HY050Dn + HY060Dn + HY070Dn + HY080Dn + HY090Dn + HY110Dn - HY100D - HY120D - HY130D - HY140D

   # Da die Quintile Share Ratio nun immer noch negativ ist, werden nun noch einmal die negativen Komponenten verändert

HY100Dn <- quantile.korr(x="HY100G",acuteness=10)
HY120Dn <- quantile.korr(x="HY120G",acuteness=10)
HY130Dn <- quantile.korr(x="HY130G",acuteness=10)
HY140Dn <- quantile.korr(x="HY140G",acuteness=10)


   # Wenn das abzuziehende Einkommen pauschal mit einem Faktor runter gewichtet ist, dann sinkt der Wert für die Quintile Share ratio noch weiter ab, deshalb müssen die Steuern nur bei den niedrigen Einkommen verringert werden

HY020new2 <- HY010.phhn + HY030Dn + HY040Dn + HY050Dn + HY060Dn + HY070Dn + HY080Dn + HY090Dn + HY110Dn - HY100Dn - HY120Dn - HY130Dn - HY140Dn

HY110Dn[is.na(HY110Dn)] <-  0 
HY090Dn[is.na(HY090Dn)] <-  0 
HY080Dn[is.na(HY080Dn)] <-  0 
HY060Dn[is.na(HY060Dn)] <-  0 
HY040Dn[is.na(HY040Dn)] <-  0 
HY030Dn[is.na(HY030Dn)] <-  0 

HY020new2[is.na(HY020new2)] <- 0
HY100Dnn <- HY100Dn

abs<-quantile(HY020new2,c(0.2,0.8),na.rm=T)
# HY100Dnn[HY020new2<abs[1] & HY100Dn>0] <- HY100Dnn[HY020new2<abs[1] & HY100Dn>0]*.5
HY100Dnn[HY020new2<abs[1]] <- HY100Dnn[HY020new2<abs[1]]*.1

HY020new2 <- HY010.phhn + HY030Dn + HY040Dn + HY050Dn + HY060Dn + HY070Dn + HY080Dn + HY090Dn + HY110Dn - HY100Dnn - HY120Dn - HY130Dn - HY140Dn

   # jetzt werden einzelne negative Ausreißer bearbeitet

ind <- which.min(HY020new2)

HY140Dnn <- HY140Dn
HY140Dnn[order(HY020new2)][1:20] <- 0

   # wer ein negatives Gesamteinkommen hat zahlt keine Einkommenssteuer
HY140Dnn[HY020new2<0] <- 0

   # wer niedriges Einkommen hat zahlt weniger Einkommenssteuer
HY140Dnn[HY020new2<abs[1]] <- 0.1*HY140Dnn[HY020new2<abs[1]]

HY020new2 <- HY010.phhn + HY030Dn + HY040Dn + HY050Dn + HY060Dn + HY070Dn + HY080Dn + HY090Dn + HY110Dn - HY100Dnn - HY120Dn - HY130Dn - HY140Dnn

   # in diesem Fall werden viel zu hohe Einkommenssteuern bezahlt
HY140Dn[ind] <- 0


   # Zur Kontrolle

abs<-quantile(HY020new2,c(0.2,0.8),na.rm=T)
t1s <- mean(HY020new2[HY020new2<=abs[1]],na.rm=T)
t2s <-mean(HY020new2[HY020new2>abs[2]],na.rm=T)
t2s/t1s

   # bei der Variable HY020 sollten keine negativen Einkommen enthalten sein, da sich sonst der Gini nicht vernünftig berechnen lässt

ind.b0 <- which(HY020new2<0)

###############################################
# ab hier Version 2
#
#
#
#

save(AGE,file="/eingabe/Ameland/universe/results_v2/AML.AGE.RData")

HY040D <- HY140Dnn

HY020 <- HY020new2
save(HY020,file="/eingabe/Ameland/universe/results_v2/AML.HY020.RData")

HY010.p <- HY010.pn
save(HY010.p,file="/eingabe/Ameland/universe/results_v2/AML.HY010.p.RData")


   # Die Quintile Share Ratio ist für den Datensatz noch zu hoch
   # diese liegt in Europa zischen 3 und 8, im Datensatz momentan ca bei 9


QSR <- function(INC){
   abs<-quantile(INC,c(0.2,0.8),na.rm=T)
   t1s <- mean(INC[INC<=abs[1]],na.rm=T)
   t2s <-mean(INC[INC>abs[2]],na.rm=T)
   return(t2s/t1s)
}

QSR(AML.HY020)

#--------------------------------------------------------------------------#
# Trunkierung der obersten Einkommen
#--------------------------------------------------------------------------#

Th <- quantile(PY010D,.995)

ind.Th <- as.numeric(PY010D>Th)


ff<-function(x){
   dens <- density(x)
   return(approxfun(x=dens$x,y=dens$y))
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

min.tinc <- quantile(PY010D,.99)
max.tinc <- quantile(PY010D,.995)


T.INC <- PY010D[PY010D>min.tinc & PY010D<max.tinc]
max.hist <- max(density(T.INC)$y)
n <- sum(ind.Th)

x1 <- rep(NA,n)
f1 <- ff(T.INC)
for (j in 1:n) {
   x1[j] <- AC.RE1()
}

PY010D[ind.Th==1] <- x1

   # jetzt sollte noch die Zahl der Nulleinkommen verringert werden


   # Es dürfen hier nur die Personen ausgewählt werden, die auch tatsächlich arbeiten gehen
# NullInc <- (1:length(PY010D))[PY010D==0 & RB210==1]

   # Das Einkommen PY010G ist nur bei den Personen null, die keiner Erwerbstätigkeit nachgehen.


HY010.pn <- PY010D + PY020D + PY050D + PY070D + PY090D + PY100D + PY110D + PY120D + PY130D + PY140D
HY010.phhn <- tapply(HY010.pn,HID,sum,na.rm=T)
HY020a <- HY010.phhn + HY030D + HY040D + HY050D + HY060D + HY070D + HY080D + HY090D + HY110D - HY100D - HY120D - HY130D # - HY140D


#--------------------------------------------------------------------------#
# Berechnung des Haushaltsäquivalenzeinkommens
#--------------------------------------------------------------------------#
region.hh <- SynthDat$REG
   # Die Variable HY140G wird zur Korrektur verwendet

   # Das C steht für correction

HY140C <- rep(NA,length(HY020a)) 
HY140C[HY010.phhn<0] <- 0
HY140C[HY010.phhn>0 & HY010.phhn<quantile(HY010.phhn,.2)] <- 0
HY140C[is.na(HY140C) & HY010.phhn<quantile(HY010.phhn,.8)] <- HY010.phhn[is.na(HY140C) & HY010.phhn<quantile(HY010.phhn,.8)]*.15
HY140C[HY010.phhn>0 & HY010.phhn>quantile(HY010.phhn,.5)] <- HY010.phhn[HY010.phhn>quantile(HY010.phhn,.5)]*.3
HY140C[HY010.phhn>0 & HY010.phhn>quantile(HY010.phhn,.6)] <- HY010.phhn[HY010.phhn>quantile(HY010.phhn,.6)]*.4
HY140C[HY010.phhn>0 & HY010.phhn>quantile(HY010.phhn,.8)& region.hh==1] <- HY010.phhn[HY010.phhn>quantile(HY010.phhn,.8)& region.hh==1]*.4
HY140C[HY010.phhn>0 & HY010.phhn>quantile(HY010.phhn,.8)& region.hh==2] <- HY010.phhn[HY010.phhn>quantile(HY010.phhn,.8)& region.hh==2]*.5
HY140C[HY010.phhn>0 & HY010.phhn>quantile(HY010.phhn,.8)& region.hh==3] <- HY010.phhn[HY010.phhn>quantile(HY010.phhn,.8)& region.hh==3]*.45
HY140C[HY010.phhn>0 & HY010.phhn>quantile(HY010.phhn,.8)& region.hh==4] <- HY010.phhn[HY010.phhn>quantile(HY010.phhn,.8)& region.hh==4]*.6


   # Bei der Variable HY110 darf es nur Einkommen geben, wenn auch Personen unter 16 im Haushalt sind
# pu16 <- tapply(AGE<16,HID,sum)

# HY110C <- HY110D
# HY110C[pu16==0] <- 0


   # HY080 sind erhaltene Einkommenstransfers

HY080C <- HY080D

HY080C[HY010.phhn>quantile(HY010.phhn,.4)] <- 0
# HY080C[HY010.phhn>quantile(HY010.phhn,.2)] <- jitter(rep(5000,length(HY080C[HY010.phhn>quantile(HY010.phhn,.2)])))
HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==1] <- jitter(rep(5000,length(HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==1])))
HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==2] <- jitter(rep(10000,length(HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==2])))
HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==3] <- jitter(rep(8000,length(HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==3])),factor=20)
HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==4] <- jitter(rep(15000,length(HY080C[HY010.phhn>quantile(HY010.phhn,.2) & region.hh==4])),factor=5)


   # In der Variable HY090 sind auch noch extreme Ausreißer enthalten


Th <- quantile(HY090D,.995)
ind.Th <- as.numeric(HY090D>Th)
min.tinc <- quantile(HY090D,.99)
max.tinc <- quantile(HY090D,.995)
T.INC <- HY090D[HY090D>min.tinc & HY090D<max.tinc]
max.hist <- max(density(T.INC)$y)
n <- sum(ind.Th)
x1 <- rep(NA,n)
f1 <- ff(T.INC)
for (j in 1:n) {
   x1[j] <- AC.RE1()
}
HY090C <- HY090D
HY090C[ind.Th==1] <- x1




HY020a <- HY010.phhn + HY030D + HY040D + HY050D + HY060D + HY070D + HY080C + HY090C + HY110C - HY100D - HY120D - HY130D - HY140C

HY120C <- HY120D
HY120C[HY020a<0] <- 0

HY130C <- HY130D
HY130C[HY020a<0] <- 0

HY100C <- HY100D
HY100C[HY020a<0] <- 0

HY140C <- HY140D
HY140C[HY020n<0] <- 0

HY020n <- HY010.phhn + HY030D + HY040D + HY050D + HY060D + HY070D + HY080C + HY090C + HY110C - HY100C - HY120C - HY130C - HY140C

HY140C <- HY140D
HY140C[HY020n<0] <- 0
   # dieses muss eventuell wiederholt werden


HY020n <- HY010.phhn + HY030D + HY040D + HY050D + HY060D + HY070D + HY080C + HY090C + HY110C - HY100C - HY120C - HY130C - HY140C

EDI <- (HY020n/Gew2)[ind.HID]

Lq <- quantile(EDI,.2)
Hq <- quantile(EDI,.8)

sum(EDI[EDI>Hq])/sum(EDI[EDI<Lq& EDI>0])

save(EDI,file="AML.EDI.RData")

   # andere korrigierte Variablen ebenfalls abspeichern

setwd("/eingabe/Ameland/universe/results_v2")

save(HY010.phhn,file="AML.HY010.phhn.RData")
save(HY080C,file="AML.HY080C.RData")
save(HY090C,file="AML.HY090C.RData")
save(HY110C,file="AML.HY110C.RData")
save(HY100C,file="AML.HY100C.RData")
save(HY120C,file="AML.HY120C.RData")
save(HY130C,file="AML.HY130C.RData")
save(HY140C,file="AML.HY140C.RData")

   # die Variable PY010C wird auch abgespeichert, vielleicht kann man sie ja noch mal brauchen
setwd("/eingabe/Ameland/universe/results_v2")

# PY010C <- PY010D
# save(PY010C,file="AML.PY010C.RData")

#--------------------------------------------------------------------------#
# Verbesserung der regionalen Heterogenitäten
#--------------------------------------------------------------------------#


library(maptools)
library(spdep)


   # Der EU-SILC Datensatz weist deutlich mehr Heterogenität als AMELIA auf, dies muss noch verbessert werden



yvar <- "PY010C"

   # Welche Variablen wurden zur Erzeugung verwendet

xvars = c("ACL","SEX","RB210","PE040","FST","DOU")

   # diese Variablen werden eingeladen

eval(parse(text=paste("load('/eingabe/Ameland/universe/results_v2/AML.",yvar,".RData')",sep="",collapse=";")))

eval(parse(text=paste("load('/eingabe/Ameland/universe/results_v1/AML.",xvars,".RData')",sep="",collapse=";")))

   # zusätzlich muss noch die Information über die Stadt eingeladen werden

load("/eingabe/Ameland/universe/results_v1/AML.CIT.RData")


eval(parse(text=paste("IDENTIFIER <- paste(",paste(xvars,collapse=",",sep=""),",sep='')",sep="")))

gr.ident <- tapply(IDENTIFIER,IDENTIFIER,length)

gr.ident.o <- sort(gr.ident,decreasing=T)

   # nun kommt die Karte ins Spiel.

   # Im vorliegenden Fall soll Heterogenität und räumliche Muster gleichzeitig erzeugt werden.

load("/eingabe/Ameland/pluria.RData")

coords <- coordinates(pluria)

IDs <- pluria@data$CIT

DAT.nb <- knn2nb(knearneigh(coords,k=4),row.names=IDs)

DAT.W <- nb2listw(DAT.nb,style="W")

eval(parse(text=paste("tab.y <- tapply(",yvar,",CIT,sum)",sep="")))

ind <- match(IDs, names(tab.y))

moran <- moran.test(tab.y[ind],listw=DAT.W)

eval(parse(text=paste("yvar.n <- ",yvar,sep="")))
eval(parse(text=paste("yvar.a <- ",yvar,sep="")))


moran.v <- vector()
moran.v[1] <- moran.test(tab.y[ind],listw=DAT.W)$statistic


for (i in 1:length(gr.ident.o)){
   yvar.a <- yvar.n
   yvar.n[IDENTIFIER==names(gr.ident.o)[i]] <- sort(yvar.n[IDENTIFIER==names(gr.ident.o)[i]])
   tab.y <- tapply(yvar.n,CIT,sum)
   moran.v[i+1] <- moran.test(tab.y[ind],listw=DAT.W)$statistic
   if (as.numeric(moran.v[i+1])<as.numeric(moran.v[i])){
      yvar.n<-yvar.a
   }
   cat(i, "\n")
}


   # Das S hinter der Variable steht für spatial

PY010S <- yvar.n

save(PY010S,file="/eingabe/Ameland/universe/results_v2/AML.PY010S.RData")

   # jetzt muss die Variable wieder auf das Äquivalenzeinkommen gespielt werden


HY010.ps <- PY010S + PY020D + PY050D + PY070D + PY090D + PY100D + PY110D + PY120D + PY130D + PY140D
HY010.phhs <- tapply(HY010.ps,HID,sum,na.rm=T)
HY020s <- HY010.phhs + HY030D + HY040D + HY050D + HY060D + HY070D + HY080C + HY090C + HY110C - HY100C - HY120C - HY130C - HY140C


EDIS <- (HY020s/Gew2)[ind.HID]

save(EDIS,file="AML.EDIS.RData")

#--------------------------------------------------------------------------#
# Haupteinkommensbezieher
#--------------------------------------------------------------------------#


A <- split(PY010C,HID)
B <- mapply(function(x)which.max(x)==1:length(x),A)
C <- unsplit(B,HID)

HEB <- C

save(HEB,file="AML.HEB.RData")


tapply(1:length(HID),HID,function(PY010C)


#--------------------------------------------------------------------------#
# Steuerklasse
#--------------------------------------------------------------------------#

TCL<-rep(NA,length(PY010C))
TCL[PY010C==0] <- 0

TCL[PY010C!=0] <- cut(PY010C[PY010C!=0],breaks=quantile(PY010C[PY010C!=0],probs=seq(0,1,1/5)),labels=F,include.lowest=T)

save(TCL,file="AML.TCL.RData")


#--------------------------------------------------------------------------#
# Variable INB
#--------------------------------------------------------------------------#


   # Koordinaten in x Richtung heraus finden

ab <- which.min(Ko.x.min)

   # Punkt der am weitesten links liegt
plot(pluria)
plot(pluria[ab,],col="black",add=T)
ab <- which.min(Ko.y.min)
   # Punkt der am weitesten südlich liegt

plot(pluria[pluria@data$DIS==4,],col="yellow")
CIT4<-pluria@data[pluria@data$DIS==4,"CIT"]

ind <- match(CIT4,pluria@data$CIT)

Ko.x <- Ko.x.min[ind]

DIS4 <- pluria[pluria@data$DIS==4,]

plot(DIS4[order(Ko.x) [1:20],],col="black",add=T)

pluria@data$SDIS<-NA
pluria@data$SDIS[pluria@data$DIS%in%c(1:2)]<-"A"
for(i in 3:40) {
   CITi<-pluria@data[pluria@data$DIS==i,"CIT"]
   ind <- match(CITi,pluria@data$CIT)
   Ko.x <- Ko.x.min[ind]
   DISi <- pluria[pluria@data$DIS==i,]
   pluria@data$SDIS[ind[order(Ko.x)][1:round(length(Ko.x)/2)]]<-"A"
   pluria@data$SDIS[ind][is.na(pluria@data$SDIS[ind])] <- "B"
}

DIS4 <- pluria[pluria@data$DIS==11,]
co<-c("red","green")
names(co)<-c("A","B")
plot(DIS4,col=co[DIS4@data$SDIS])

i<-40
DIS11<-pluria[pluria@data$DIS==i,];CIT11<-DIS11@data[,"CIT"];plot(DIS11,col=co[DIS11@data$SDIS]);i<-i+1;

i<-1
cat(i,"\n");plot(DIS11[DIS11@data[,"CIT"]==CIT11[i],],col="black",add=T);i<-i+1;
cat(i,"\n");plot(DIS11[DIS11@data[,"CIT"]==CIT11[i],],col="black",add=T);i<-i-1;
plot(DIS11[DIS11@data[,"CIT"]==CIT11[23],],col="orange",add=T)



pluria@data[pluria@data$DIS==11,"SDIS"][16]<-"B"
pluria@data[pluria@data$DIS==11,"SDIS"][41]<-"B"

pluria@data[pluria@data$DIS==16,"SDIS"][6]<-"B"
pluria@data[pluria@data$DIS==18,"SDIS"][16]<-"B"
pluria@data[pluria@data$DIS==27,"SDIS"][43]<-"A"
pluria@data[pluria@data$DIS==27,"SDIS"][34]<-"A"
pluria@data[pluria@data$DIS==33,"SDIS"][31]<-"B"
pluria@data[pluria@data$DIS==39,"SDIS"][23]<-"B"
pluria@data[pluria@data$DIS==40,"SDIS"][5:6]<-"A"

pluria@data$INB<-paste(pluria@data$DIS,pluria@data$SDIS,sep="")
ab <- as.factor(pluria@data$INB)
levels(ab) <- 1:length(unique(ab))
pluria@data$INB<- ab

load('/home/kolb/Desktop/Ameli/AmeLand/data/SynthDat/AML.CIT.RData' )

ind2 <- match(CIT,pluria@data$CIT)

INB <- pluria@data$INB[ind2]

setwd("/home/kolb/Desktop/Ameli/AmeLand/data/SynthDat")

save(INB,file="AML.INB.RData")

#--------------------------------------------------------------------------#
# Stratas einteilen
#--------------------------------------------------------------------------#

get.reord <- function(x){
	lfn <- 1:length(x)
	lfn <- lfn[order(x)]
	return(lfn)
}


INB <- as.numeric(INB)
INB2HH <- tapply(INB,HID,function(x)x[1])
HH2HH <- tapply(HID,HID,function(x)x[1])

SIND <- rep(NA,length(INB2HH))
SINDI <- rep(NA,length(INB2HH))
Len <- length(table(INB))

for (i in 1:Len){
   inc2 <- HY020[INB2HH==i]
   hid <- HH2HH[INB2HH==i]
   inc <- inc[order(inc2)]
   hid2 <- hid[order(inc2)]   
   sinda <- cut(1:length(inc),quantile(1:length(inc),probs=seq(0,1,1/strat)),labels=F,include.lowest=T)
   ind <- match(hid,hid2)
   SINDI[INB2HH==i] <- sinda[ind]
   cat(i, "\n")
}


   # noch mal neu, weil die Zuordnung oben nicht funktioniert hat

for (i in 1:Len){
   info <- length(SIND[INB2HH==i])
   SIND[INB2HH==i] <-  cut(1:info,quantile(1:info,probs=seq(0,1,1/strat)),labels=F,include.lowest=T)
   cat(i, "\n")
}

SIND.rand <- SIND
save(SIND,file="AML.SIND_rand.RData")


SIND <- SINDI
save(SIND,file="AML.SIND_sort.RData")


save(INB2HH,file="AML.INB2HH.RData")
#####

   # In einer Schicht soll die Varianz besonders groß sein, in mindestens einer anderen Schicht soll die Varianz demgegenüber besonders klein sein. 

   # jedes fünfte Element wird verwendet


SIND <- rep(NA,length(INB2HH))

for (i in 1:Len){
   inc2 <- HY020[INB2HH==i]

      # aufsteigende Ordnung 
   oi <- order(inc2)
      # jedes fünfte Element wird heraus genommen
   u25 <- sum(inc2>200000)

   indfirst <- seq(1,length(inc2),by=5)
   indfirst <- head(indfirst,length(indfirst)-u25)
   indfi     <- unique(c(oi[indfirst],which(inc2>200000)))
   SIND[INB2HH==i][indfi] <- 1
#[which(inc2>250000)]
      # nun kommen zei eher homogene Schichten
   oi <- oi[!(oi%in%indfi)]
   inds <- cut(1:length(oi),quantile(1:length(oi),probs=seq(0,1,1/4)),labels=F,include.lowest=T)
   SIND[INB2HH==i][oi[inds==1]] <- 2
   SIND[INB2HH==i][oi[inds==2]] <- 3
   oi <- oi[-which(inds%in%c(1,2))]

   indt <- sample(oi,round(length(oi)/2),replace=F)
   SIND[INB2HH==i][oi[which(oi%in%indt)]] <- 4
   oi <- oi[-which(oi%in%indt)]
   SIND[INB2HH==i][oi] <- 5
   cat(i, "\n")
}
summary(t(tapply(HY020,list(SIND,INB2HH),sd)))
      # der Rest wird in drei Schichten aufgeteilt

   


#--------------------------------------------------------------------------#
# Imputation
#--------------------------------------------------------------------------#



   # Das Einkommen der obersten 20 % muss verringert werden

ind <- PY010D>quantile(PY010D,.8)

fehl <- rep(NA,sum(ind==T))

PY010I <- c(fehl,PY010G)




   # Variablen die zur Imputation verwendet werden
inc.vars <- c("PY010","PY020","PY050","PY070","PY090","PY100","PY110","PY120","PY130","PY140")

setwd(AML.path)
eval(parse(text=paste("load('AML.",inc.vars,"D.RData')",sep="",collapse=";")))

setwd(EUS.path)
eval(parse(text=paste("load('",inc.vars,"G.cs.RData')",sep="",collapse=";")))


   # die erste Variable wird anders behandelt
inc.vars <- c("PY020","PY050","PY070","PY090","PY100","PY110","PY120","PY130","PY140")

   # zu ersetzende Einheiten und Originalvariablen werden aneinander gehängt
eval(parse(text=paste(inc.vars[-1],"I <- c(",inc.vars[-1],"D[ind],",inc.vars[-1],"G)",sep="")))

eval(parse(text=paste("IDat <- cbind(PY010I,",paste(inc.vars[-1],"I",sep="",collapse=","),")",sep="")))

IDat.m <- mice(IDat)
IDat.c <- complete(IDat.m)

   # Achtung die Imputation dauert sehr lange


save(HY080C,file="AML.HY080C.RData") 
save(HY090C,file="AML.HY090C.RData") 
save(HY100C,file="AML.HY100C.RData")  
save(HY110C,file="AML.HY110C.RData")    
save(HY120C,file="AML.HY120C.RData")  
save(HY130C,file="AML.HY130C.RData")   
save(HY140C,file="AML.HY140C.RData")   
save(PY010C,file="AML.PY010C.RData")  


scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY080C.RData AML.HY080C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY090C.RData AML.HY090C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY100C.RData AML.HY100C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY110C.RData AML.HY110C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY120C.RData AML.HY120C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY130C.RData AML.HY130C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.HY140C.RData AML.HY140C.RData
scp -r kolb@192.168.1.48:/eingabe/Ameland/universe/results_v2/AML.PY010C.RData AML.PY010C.RData

   # Es wird eine neue Variable eingeführt - die Verschuldung

save(HY100C,file="AML.HY100C.RData")

######################################
# es muss eine homogenere Variable als CIT erzeugt werden


gleiche_laenge<- function(x){
   maxanz<-max(nchar(x))
   pst <- function(x){if(nchar(x)<=maxanz){paste(paste(rep(0,maxanz-nchar(x)),collapse=""),x,sep="")}else{x}}
   y<- sapply(x,pst)
   return(y)
}

load("AML.CIT.RData")
load("AML.DIS.RData")
load("AML.HID.RData")
load("AML.DOU.RData")
load("AML.CITh.RData")

CIT2 <- gleiche_laenge(CIT)

CIT2DIS <- tapply(CIT,DOU,table)

CITh <- rep(0,length(CIT))
#CITh[CIT==1213] <- 1
#CITh[CIT==721] <- 2

HID2CIT <- tapply(HID[CIT%in%c(1213,721)],CIT[CIT%in%c(1213,721)],table)

CITh.i2 <- cut(1:length(HID2CIT[[1]]),breaks = quantile(1:length(HID2CIT[[1]]),seq(0,1,length.out=5)), labels=FALSE,include.lowest=TRUE)
CITh.i2[CITh.i2==5] <- 4
ind <- match(HID,names(HID2CIT[[1]]))
CITh[is.na(ind)==F] <- CITh.i2[ind][is.na(ind)==F]
k <- max(CITh.i2)

CITh.i3 <- cut(1:length(HID2CIT[[2]]),breaks = quantile(1:length(HID2CIT[[2]]),seq(0,1,length.out=5)), labels=FALSE,include.lowest=TRUE)
CITh.i3[CITh.i3==5] <- 4
ind <- match(HID,names(HID2CIT[[2]]))
CITh[is.na(ind)==F] <- CITh.i3[ind][is.na(ind)==F]+k
k <- max(CITh)

CIT2DIS <- tapply(CIT,paste(DIS,DOU,sep=""),table)

for (i in 3:length(CIT2DIS)){
   comi <- CIT2DIS[[i]]
   if (length(comi)>4){
      CITh.i <- cut(sort(comi),breaks = quantile(comi,seq(0,1,length.out=round(sum(comi)/15000))), labels=FALSE,include.lowest=TRUE)
      ind <- match(CIT,names(sort(comi)))
      CITh[is.na(ind)==F] <- CITh.i[ind][is.na(ind)==F]+k
   }else{
      ind <- match(CIT,names(comi))
      CITh[is.na(ind)==F] <- (1:length(comi))[ind][is.na(ind)==F]+k
      CITh.i <- length(comi) + k
   }
     k <- k+max(CITh.i)
     cat( i, "\n")
}

lapply(CIT2DIS,length)
lapply(CIT2DIS,sum)
table(CITh==0)

table(CITh[DOU==1&DIS==10])

CITh[CITh==0&DOU==1&DIS==10] <- 3889

table(CITh[DOU==2&DIS==10])

CITh[CITh==0&DOU==2&DIS==10] <- 3890

DIS[CITh==1803]
table(DOU[CITh==1803])

table(CITh[DOU==1&DIS==3])

table(CITh[DOU==3&DIS==10])

CITh[CITh==1803] <- 1804
CITh[CITh==1803] <- 1804
CITh[CITh==1805] <- 1804

CITh[CITh==41&DOU==2] <- 3891


tab.CITh <- table(CITh)

   # 200 Gemeinden

sd(tab.CITh)/mean(tab.CITh)

CIT2HID <- tapply(CIT,HID,table)

DIS[CITh==3875]
DOU[CITh==3875]
table(CITh[DOU==1&DIS==9])
CITh[CITh==3875] <- 3876

DIS[CITh==3850]
DOU[CITh==3850]
table(CITh[DOU==1&DIS==7])
CITh[CITh==3850] <- 3851

DIS[CITh==3761]
DOU[CITh==3761]
table(CITh[DOU==1&DIS==38])
CITh[CITh==3761] <- 3762

DIS[CITh==122]
DOU[CITh==122]
table(CITh[DOU==3&DIS==14])
CITh[CITh==122] <- 123

DIS[CITh==1692]
DOU[CITh==1692]
table(CITh[DOU==1&DIS==22])
CITh[CITh==1692] <- 123


save(CITh,file="AML.CITh.RData")

k <- 1
for (i in 1:66){
   ind <- names(sort(tab.CITh)[k])
   #ind <- names(which.min(tab.CITh))
   DIS.ind <- unique(DIS[CITh==ind])
   DOU.ind <- unique(DOU[CITh==ind])
   tabDOUDIS <- table(CITh[DOU==DOU.ind&DIS==DIS.ind])
   tabDOUDIS <- tabDOUDIS[-which(names(tabDOUDIS)==ind)]
   if(length(tabDOUDIS)>0){
      CITh[CITh==ind] <- names(which.min(tabDOUDIS))
      tab.CITh <- table(CITh)
   }else{k<-k+1}   
   cat("Gemeinden: ",length(tab.CITh),"Wert: ",sd(tab.CITh)/mean(tab.CITh),"\n")
}


proof <- tapply(DOU,CITh,function(x)length(unique(x)))
proof2 <- tapply(DIS,CITh,function(x)length(unique(x)))


save(CITh,file="AML.CITh.RData")

