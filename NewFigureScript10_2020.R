library(reshape2)
library(ggplot2)
theme_set(theme_bw())
theme_eg=theme_update(panel.grid.minor=element_blank(), 
                      panel.grid.major=element_blank(), 
                      strip.background=element_blank())
library(car)
library(vegan)
library(lme4)
library(lmtest)
library(zoo) 
library(emmeans)
library(effsize)
library(multcompView)
library(ggplotgui)
library(dplyr)

#setwd("G:/My Drive/HeatingRings") #or sometimes D: 
raw <- read.csv("HeatingRingsAGBData.csv")

#subset of temps, starting from 6/16 (beginning of expt) to 8/31 (end)
temp_data <- read.csv("HeatingRingsTempData2016_subset.csv")


#for analysis (will need to relabel midwest later)
m.raw <- melt(raw, id=c("pot", "heating.ring", "trt", "ring", "population", "soil", "Inocula.sour"), 
              measure=c("CORLAN", "DALPUR", "ECHPUR", "LESCAP", "RUDHIR", "RATPIN", "BOUCUR", "ANDGER", 
                        "SCHSCO", "KOEMAC"))
m.raw$ring <- as.factor(m.raw$ring)
levels(m.raw$soil) <- c("Autoclaved", "Local", "Midwest", "Oldfield", "Southern")

levels(m.raw$population) <- c("Local seed", "Midwest seed", "Southern seed")

#if you want to include the autoclaved treatment, run this line:
#m.raw=m.raw[m.raw$soil %in% c("autoclaved", "local remnant", "midwest remnant", "southern remnant"),]

#if you want to exclude the autoclaved treatment, run this line: 
m.raw <- m.raw[m.raw$soil %in% c("Local", "Midwest", "Southern"),]
#line no longer works (it's removing all lines not considering assigned levels)
#this works 
m.raw <- m.raw[m.raw$soil %in% c("local remnant","midwest remnant","southern remnant"),]


m.raw.a <- dcast(m.raw, heating.ring+trt+ring+population+soil+variable~., mean) 
names(m.raw.a)[7] <- "biomass"
m.raw.a$logbiomass <- log(m.raw.a$biomass+0.001)
names(m.raw.a)[8] <- "logbiomass"
m.raw.a$species <- m.raw.a$variable
levels(m.raw.a$species) <- c("Coreopsis", "Dalea", "Echinacea", "Lespedeza", "Rudbeckia", 
                          "Ratibida", "Bouteloua", "Andropogon", "Schizachyrium", "Koeleria")






######--------community figures--------#######
#averaging across pseudoreplicates within a block
matrixform.a <- dcast(m.raw, heating.ring+trt+ring+population+soil~variable, fun=mean)
#beta div 
newcol <- interaction(matrixform.a$trt,matrixform.a$population)
matrixform.a$trtXpop <- newcol
#dist matrix 
bray <- vegdist(matrixform.a[6:14], method="bray",bias.adjust=T)
#trt by pop 
dis.bray <- betadisper(bray,matrixform.a$trtXpop)
#get a warning for squared distances being negative and changed to zero 
#ignore for now 
permutest(dis.bray,pairwise=T)
dis.bray
dis.bray$distances
interaction.short <- data.frame(trt=matrixform.a$trt,population=matrixform.a$population,trtXpop=matrixform.a$trtXpop)
interaction.short$median <- dis.bray$distances

#make dataframe
com <- matrixform.a[,c("heating.ring", "trt", "ring", "population", "soil")]
com$simpsons <- diversity(matrixform.a[,c("CORLAN", "DALPUR", "ECHPUR", "LESCAP", "RUDHIR", 
                                       "RATPIN", "BOUCUR", "ANDGER", "SCHSCO")], 
                       index="simpson")
com$prod <- apply(matrixform.a[,c("CORLAN", "DALPUR", "ECHPUR", "LESCAP", "RUDHIR", 
                               "RATPIN", "BOUCUR", "ANDGER", "SCHSCO")], 
               MARGIN=1, FUN <- sum)

com$median <- interaction.short$median
com$trtXpop <- interaction.short$trtXpop


#melt 
m.com <- melt(com, id=c("trt", "ring", "population", "soil","trtXpop"),measure=c("simpsons","median","prod"))
comfig <- dcast(m.com, trt+population+trtXpop+variable~., fun=mean)
names(comfig)[5] <-"average"
comfig.se <- dcast(m.com, trt+population+trtXpop+variable~., fun=sd)
comfig.se$se <- comfig.se$./sqrt(12)
comfig.se$. <- NULL
comfig <- merge(comfig, comfig.se, by=c("trt", "population", "trtXpop","variable"), na.rm=TRUE)

names(m.com)[7] <- "average"


qplot(trt, average, data=comfig, fill=trt, shape=population, xlab="", size=I(4.2)) +
  facet_grid(variable~population, scales="free") +
  geom_errorbar(aes(ymin=average-se, ymax=average+se, width=0.2)) +
  theme(legend.position="none") +
  scale_shape_manual(values=c(21,22,24)) +
  scale_fill_manual(values=c("gray", "black")) +
  geom_jitter(mapping=aes(x=trt, y=average, shape=population),data=m.com,size=I(1),width=0.04,height = 0.04)
ggsave("FigureImages/comfig.png")


##supplemnetary figs on nonsig effects of soil on comm variables 
m.com <- melt(com, id=c("trt", "ring", "population", "soil"), measure=c("prod"))
comfig <- dcast(m.com, soil+variable~., fun=mean)
names(comfig)[3] <- "average"
comfig.se <- dcast(m.com, soil+variable~., fun=sd)
comfig.se$se <- comfig.se$./sqrt(12)
comfig.se$. <- NULL
comfig <- merge(comfig, comfig.se, by=c("soil","variable"))

names(m.com)[6] <- "average"

qplot(soil, average, data=comfig, shape=soil, xlab="", size=I(4.2)) +
  facet_grid(~variable, scales="free") + 
  geom_errorbar(aes(ymin=average-se, ymax=average+se, width=0.2)) + 
  theme(legend.position="none") +
  scale_shape_manual(values=c(21,22,24)) +
  geom_jitter(mapping=aes(x=soil, y=average, shape=soil),data=m.com,size=I(1),width=0.04,height = 0.04)

ggsave("FigureImages/NonSigSoilProd.png")


m.com <- melt(com, id=c("trt", "ring", "population", "soil"), measure=c("simpsons"))
comfig <- dcast(m.com, soil+variable~., fun=mean)
names(comfig)[3] <- "average"
comfig.se <- dcast(m.com, soil+variable~., fun=sd)
comfig.se$se <- comfig.se$./sqrt(12)
comfig.se$.NULL
comfig <- merge(comfig, comfig.se, by=c("soil","variable"))

names(m.com)[6]="average"

qplot(soil, average, data=comfig, shape=soil, xlab="", size=I(4.2)) + 
  facet_grid(~variable, scales="free") +
  geom_errorbar(aes(ymin=average-se, ymax=average+se, width=0.2)) + 
  theme(legend.position="none") + 
  scale_shape_manual(values=c(21,22,24)) +
  geom_jitter(mapping=aes(x=soil, y=average, shape=soil),data=m.com,size=I(1),width=0.04,height = 0.04)


ggsave("FigureImages/NonSigSoilSimpsons.png")



#supplemental soil beta div plot
bray <- vegdist(matrixform.a[6:14], method="bray",bias.adjust=T)
dis.braysoil <- betadisper(bray,matrixform.a$soil, bias.adjust = T)
matrixform.a$median <- dis.braysoil$distances

melt.intx <- melt(matrixform.a,id=c("soil"),measure="median")

cast.intx1 <- dcast(melt.intx,soil~.,fun=mean)
names(cast.intx1) <- c("soil","median")
cast.intx2 <- dcast(melt.intx, soil~.,fun=sd)
names(cast.intx2) <- c("soil","sd.disp")
cast.intx3 <- dcast(melt.intx, soil~., fun=length)
names(cast.intx3) <- c("soil", "n")
cast.intx4 <- merge(cast.intx1, cast.intx2, by=c("soil"))
cast.intx <- merge(cast.intx4, cast.intx3, by=c("soil"))
cast.intx$se <- cast.intx$sd/sqrt(cast.intx$n)

#ylim=c(0.26,0.355)
qplot(soil, median, data=cast.intx, shape=soil, size=I(4.2)) +
  geom_errorbar(aes(ymin=median-se, ymax=median+se, width=0.2)) +
  theme(legend.position="none", panel.background=element_rect()) +
  scale_shape_manual(values=c(21,22,24)) + 
  geom_jitter(mapping=aes(x=soil, y=median, shape=soil),data=matrixform.a,size=I(1),width=0.04,height = 0.04) 

ggsave("FigureImages/NonSigBetaSoil.png")


#####---temp by provenance plots----######
forfig3 <- dcast(m.raw.a[m.raw.a$variable %in% c("CORLAN", "ECHPUR", "RATPIN", "RUDHIR", "DALPUR", "LESCAP","ANDGER","BOUCUR","SCHSCO"),], 
              trt+population+species~., fun=mean, value.var="logbiomass")
names(forfig3)[4] <- "logmean"
forfig3.se <- dcast(m.raw.a[m.raw.a$variable %in% c("CORLAN", "ECHPUR","RATPIN", "RUDHIR", "DALPUR", "LESCAP","ANDGER","BOUCUR","SCHSCO"),], 
                 trt+population+species~., fun=sd, value.var="logbiomass")
forfig3.se$logse <- forfig3.se$./2
forfig3.se$. <- NULL
forfig3 <- merge(forfig3, forfig3.se, by=c("trt", "population", "species"))
forfig3$backmean <- exp(forfig3$logmean)
forfig3$lowerse <- exp(forfig3$logmean-forfig3$logse)
forfig3$upperse <- exp(forfig3$logmean+forfig3$logse)

#subset for raw values of spp to include in each plot 
multspp <- m.raw[m.raw$variable  %in% c("CORLAN", "ECHPUR", "RUDHIR", "RATPIN","DALPUR", "LESCAP","ANDGER","BOUCUR","SCHSCO"),]
names(multspp)[9] <- "biomass"
multspp$species[multspp$variable=="CORLAN"]<-"Coreopsis"
multspp$species[multspp$variable=="ECHPUR"]<-"Echinacea"
multspp$species[multspp$variable=="RATPIN"]<-"Ratibida"
multspp$species[multspp$variable=="RUDHIR"]<-"Rudbeckia"
multspp$species[multspp$variable=="DALPUR"]<-"Dalea"
multspp$species[multspp$variable=="LESCAP"]<-"Lespedeza"
multspp$species[multspp$variable=="ANDGER"]<-"Andropogon"
multspp$species[multspp$variable=="BOUCUR"]<-"Bouteloua"
multspp$species[multspp$variable=="SCHSCO"]<-"Schizachyrium"

multspp$logbiomass <- log(multspp$biomass+0.001)
multspp$backmean <- exp(multspp$logbiomass)

#straight raw values 
# qplot(trt, backmean, data=forfig3, ylab="biomass (g)", xlab="", fill=trt, shape=population, size=I(2.5)) +
#   facet_grid(species~population, scales="free") + 
#   geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.15)) + 
#   theme(strip.text.y=element_text(face="italic"), legend.position="none") + 
#   scale_fill_manual(values=c("gray", "black")) +
#   scale_shape_manual(values=c(21,22,24)) +
#   geom_point(mapping=aes(x=trt, y=backmean, shape=population),data=multspp,size=I(0.4)) 


#a little jitter, going with this
qplot(trt, backmean, data=forfig3, ylab="biomass (g)", xlab="", fill=trt, shape=population, size=I(3)) +
  facet_grid(species~population, scales="free") + 
  geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.25)) + 
  theme(strip.text.y=element_text(face="italic"), legend.position="none") + 
  scale_fill_manual(values=c("gray", "black")) +
  scale_shape_manual(values=c(21,22,24)) +
  geom_jitter(mapping=aes(x=trt, y=backmean, shape=population),data=multspp,size=I(0.55),width=0.04, height = 0.04) 

# #boxplot? 
# qplot(trt, backmean, data=forfig3, ylab="biomass (g)", xlab="", fill=trt, size=I(4.2)) +
#   facet_grid(species~population, scales="free") + 
#   geom_boxplot(color="black") +
#   #geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.2)) + 
#   theme(strip.text.y=element_text(face="italic"), legend.position="none") + 
#   scale_fill_manual(values=c("gray", "black")) +
#   #scale_shape_manual(values=c(21,22,24)) +
#   geom_jitter(mapping=aes(x=trt, y=backmean, shape=population),data=multspp,size=I(0.7),width=0.04,height = 0.04) 
# 
# 
# ggplot(forfig3,aes(x=trt,y=backmean, fill=trt)) +
#   geom_boxplot(color="black") +
#   facet_grid(species~population,scales="free") +
#   scale_fill_manual(values=c("gray", "black")) +
#   geom_jitter(mapping=aes(x=trt, y=backmean, shape=population),data=multspp,size=I(0.7),width=0.04,height = 0.04)
##nope 
ggsave("FigureImages/Fig3.png")

#boxplot
# qplot(trt, backmean, data=forfig3, ylab="biomass (g)", xlab="",fill=trt) +
#   facet_grid(species~population, scales="free") +
#   #geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.15)) +
#   theme(strip.text.y=element_text(face="italic"), legend.position="none") +
#   scale_fill_manual(values=c("gray", "black")) +
#   scale_shape_manual(values=c(21,22,24)) +
#   geom_boxplot(mapping=aes(x=trt, y=backmean, shape=population),data=multspp,size=I(0.4),width=0.2,alpha=0.5,outlier.size = -1) +
#   geom_point(mapping=aes(x=trt, y=backmean, shape=population),data=multspp,size=I(0.4))


# ###fig 4
# forfig4=dcast(m.raw.a[m.raw.a$variable %in% c("BOUCUR", "RATPIN"),], 
#               trt+population+species~., fun=mean, value.var="logbiomass")
# names(forfig4)[4]="logmean"
# forfig4.se=dcast(m.raw.a[m.raw.a$variable %in% c("BOUCUR", "RATPIN"),], 
#                  trt+population+species~., fun=sd, value.var="logbiomass")
# forfig4.se$logse=forfig4.se$./2
# forfig4.se$.=NULL
# forfig4=merge(forfig4, forfig4.se, by=c("trt", "population", "species"))
# forfig4$backmean=exp(forfig4$logmean)
# forfig4$lowerse=exp(forfig4$logmean-forfig4$logse)
# forfig4$upperse=exp(forfig4$logmean+forfig4$logse)
# 
# #subset for raw values of spp to include in each plot 
# RB=m.raw[m.raw$variable  %in% c("RATPIN","BOUCUR"),]
# names(RB)[9]="biomass"
# RB$species<-ifelse(RB$variable=="RATPIN","Ratibida","Bouteloua")
# RB$logbiomass=log(RB$biomass+0.001)
# RB$backmean=exp(RB$logbiomass)
# 
# qplot(trt, backmean, data=forfig4, ylab="biomass (g)", xlab="", fill=trt, shape=population, size=I(4)) +
#   facet_grid(species~population, scales="free") + 
#   geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.2)) + 
#   theme(strip.text.y=element_text(face="italic"), legend.position="none") + 
#   scale_fill_manual(values=c("gray", "black")) +
#   scale_shape_manual(values=c(21,22,24)) +
#   geom_jitter(mapping=aes(x=trt, y=backmean, shape=population),data=RB,size=I(0.7),width=0.04,height = 0.04)
# 
# ggsave("FigureImages/Fig4.png")


#####--------soil legume fig-------#####
forfig4 <- dcast(m.raw.a[m.raw.a$variable %in% c("LESCAP", "DALPUR"),], 
              soil+species~., fun=mean, value.var="logbiomass")
names(forfig4)[3] <- "logmean"
forfig4.se <- dcast(m.raw.a[m.raw.a$variable %in% c("LESCAP", "DALPUR"),], 
                 soil+species~., fun=sd, value.var="logbiomass")
forfig4.se$logse <- forfig4.se$./sqrt(24)
forfig4.se$. <- NULL
forfig4 <- merge(forfig4, forfig4.se, by=c("soil", "species"))
forfig4$backmean <- exp(forfig4$logmean)
forfig4$lowerse <- exp(forfig4$logmean-forfig4$logse)
forfig4$upperse <- exp(forfig4$logmean+forfig4$logse)

#subset for raw values of spp to include in each plot 
LD <- m.raw[m.raw$variable  %in% c("LESCAP","DALPUR"),]
names(LD)[9] <- "biomass"
LD$species<-ifelse(LD$variable=="LESCAP","Lespedeza","Dalea")
LD$logbiomass <- log(LD$biomass+0.001)
LD$backmean <- exp(LD$logbiomass)

qplot(soil, backmean, data=forfig4, ylab="biomass (g)", xlab="Soil inoculum", shape=soil,size=I(4)) + 
  facet_wrap(~species, scales="free") + 
  geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.15)) + 
  scale_shape_manual(values=c(21,22,24)) + 
  theme(strip.text.x=element_text(face="italic"),legend.position="none") +
  geom_jitter(mapping=aes(x=soil, y=backmean, shape=soil),data=LD,size=I(0.5),width=0.04, height = 0.04)
#think about outlier for N dalea in future 

ggsave("FigureImages/Fig4.png") 



#####-------andger soil----#####
forfig6=dcast(m.raw.a[m.raw.a$variable %in% c("ANDGER"),], 
              population+soil+species~., fun=mean, value.var="logbiomass")
names(forfig6)[4]="logmean"
forfig6.se=dcast(m.raw.a[m.raw.a$variable %in% c("ANDGER"),], 
                 population+soil+species~., fun=sd, value.var="logbiomass")
forfig6.se$logse=forfig6.se$./sqrt(8)
forfig6.se$.=NULL
forfig6=merge(forfig6, forfig6.se, by=c("population", "soil", "species"))
forfig6$backmean=exp(forfig6$logmean)
forfig6$lowerse=exp(forfig6$logmean-forfig6$logse)
forfig6$upperse=exp(forfig6$logmean+forfig6$logse)

#subset for raw values of spp to include in each plot 
A=m.raw[m.raw$variable  %in% c("ANDGER"),]
names(A)[9]="biomass"
A$species[A$variable=="ANDGER"]<-"Andropogon"
A$logbiomass=log(A$biomass+0.001)
A$backmean=exp(A$logbiomass)


qplot(soil, backmean, data=forfig6, ylab="biomass (g)", xlab="Soil inoculum", shape=population, size=I(3.5)) + 
  facet_grid(species~population, scales="free") + 
  geom_errorbar(aes(ymin=lowerse, ymax=upperse, width=0.2)) + 
  theme(strip.text.y=element_text(face="italic"), legend.position="none") +
  scale_shape_manual(values=c(21,22,24)) + 
  geom_jitter(mapping=aes(x=soil, y=backmean, shape=population),data=A,size=I(0.7),width=0.04,height = 0.04)

ggsave("FigureImages/fig6.png")


########----- supplemental temp data----#### 
hist(temp_data$h_temp_Avg) #mean elevated temp readings in addition to ambient?
hist(temp_data$u_temp_Avg) #mean ambient avg temps? 



#differenced to get true elevated temp readings 
temp_data$temp_elev <- temp_data$h_temp_Avg - temp_data$u_temp_Avg
hist(temp_data$temp_elev,
     breaks = 40)
qplot(temp_data$temp_elev,
      geom = "histogram",
      binwidth = 0.5,
      col = I("black"),
      fill = I("grey"),
      main = "Hourly Temperature Readings, Summer 2016",
      xlab = "Temperature (°C)",
      ylab = "Count") +
      #xlim = c(-6,6),
      #breaks = seq(-6,6,0.5)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(-6,6, by = 1), limits=c(-6,6)) +
  geom_vline(aes(xintercept = median(temp_data$temp_elev)), col ='red',size=1) +
  geom_vline(aes(xintercept = mean(temp_data$temp_elev)), col ='black',size=1, linetype = "dashed")

#just do ggplot instead 

target <- subset(temp_data, temp_elev > 2.5 & temp_elev < 3.5)
nrow(target)/nrow(temp_data)
# about 50% of the time within 0.5°C within target temp
# at least 1°C above 
target1 <- subset(temp_data, temp_elev > 1)
nrow(target1)/nrow(temp_data)
