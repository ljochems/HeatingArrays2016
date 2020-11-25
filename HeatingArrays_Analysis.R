################################################################
######### libraries for analyses and plots #######
library(reshape2) 
library(ggplot2)
library(car)
library(vegan)
library(lme4) 
library(lmtest)
library(MVN)
library(zoo) 
library(emmeans)
library(effsize)
library(multcompView)
library(MVN)
library(sjPlot)


#setwd() to directory where HeatingArraysABG_Data.csv is stored 
#read in raw data
raw <- read.csv("HeatingArraysAGB_Data.csv")


#set up dataframe for multivariate analyses 
m.raw <- melt(raw, id=c("pot", "heating.ring", "trt", "ring", "population", "soil", "Inocula.sour"), 
              measure=c("CORLAN", "DALPUR", "ECHPUR", "LESCAP", "RUDHIR", "RATPIN", "BOUCUR", "ANDGER", 
                        "SCHSCO", "KOEMAC"))
m.raw$ring <- as.factor(m.raw$ring)
m.raw$soil <- as.factor(m.raw$soil) #soil provenance 
m.raw$population <- as.factor(m.raw$population) #population refers to seed provenance 
levels(m.raw$population) <- c("Local seed", "Midwest seed","Southern seed")


#####-------------MULTIVARIATE Analyses---------######
#averaging across pseudoreplicates within a block
matrixform.a <- dcast(m.raw, heating.ring+trt+ring+population+soil~variable, fun=mean)
as.data.frame(matrixform.a)

#test for multivariate normality using library(MVN)
#subset to spp responses 
spp <- matrixform.a[6:14]
#mardia for skewness and kurtosis, and low sample size 
#indicates whether dataset follows mvn at sig level of 0.05 
result_mard <- mvn(spp, mvnTest = "mardia", multivariatePlot = "qq")
result_hz <- mvn(spp, mvnTest ="hz", multivariatePlot = "qq")
result_roys <- mvn(spp,mvnTest="royston", multivariatePlot = "qq")
result_dh <- mvn(spp,mvnTest = "dh", multivariatePlot = "qq")
result_en <- mvn(spp,mvnTest="energy", multivariatePlot = "qq")

all <- data.frame(result_mard$multivariateNormality,result_hz$multivariateNormality,
                  result_roys$multivariateNormality,result_dh$multivariateNormality,result_en$multivariateNormality)
#mvn assumption rejected by all tests 


#adonis2 (from vegan package) allows marginal significance tests (not sequential) 
#response are the species columns 
adonis2(matrixform.a[6:14]~trt*population*soil, data=matrixform.a, 
        method="bray", by="margin")
#removing nonsignificant interactions:
adonis2(matrixform.a[6:14]~trt + population + soil + trt:population + population:soil + trt:soil, 
        data=matrixform.a, method="bray", by="margin")
adonis2(matrixform.a[6:14]~trt + population + soil, 
                 data=matrixform.a, method="bray", by="margin")

#dbrda (same as capscale except handles negative eigenvalues better and nonparamatric) w/ matrixform.a
#should be same as adonis2 (and it is!)
db <- dbrda(matrixform.a[6:14]~trt*population*soil+Condition(ring),
            data=matrixform.a, dist="bray"); Anova(db, by="margin")
#removing NS interactions
db <- dbrda(matrixform.a[6:14]~trt+population+soil + trt:population + trt:soil + population:soil + Condition(heating.ring), 
            data=matrixform.a,dist="bray");Anova(db, by="margin")
db <- dbrda(matrixform.a[6:14]~trt+population+soil + Condition(ring), 
            data=matrixform.a, dist="bray"); Anova(db, by="margin")

#how to plot constrained ordinations in qplot, by extracting cap1 & 2 coordinates (best shows population and trt clustering)
dbrdaresult <- scores(db)
matrixform.a$dbRDA1 <- dbrdaresult$sites[,1]
matrixform.a$dbRDA2 <- dbrdaresult$sites[,2]

#plot first two axes 
qplot(dbRDA1, dbRDA2, data=matrixform.a, fill=trt,shape=population, size=I(2.1)) + 
  scale_fill_manual(values=c("gray", "black")) +
  scale_shape_manual(values=c(21,22,24)) +
  theme(legend.title=element_blank())
ggsave("emily figures/fig1.pdf", width=4.5, height=3)
##for population legend 
qplot(dbRDA1, dbRDA2, data=matrixform.a, shape=I(21), fill=trt) +
  scale_fill_manual(values=c("gray", "black")) +
  theme(legend.title=element_blank())
ggsave("emily figures/fig1legend.pdf", width=4.5, height=3)

qplot(dbRDA1, dbRDA2, data=matrixform.a, color=trt, size=I(2.1)) + 
  scale_color_manual(values=c("gray", "black")) +
  #scale_shape_manual(values=c(21,22,24)) +
  theme(legend.title=element_blank())

#does warming reduce differences among pops (ie Beta diversity)? 
bray <- vegdist(matrixform.a[6:14], method="bray")
#trt 
dis.bray <- betadisper(bray,matrixform.a$trt, bias.adjust=T)
permutest(dis.bray, pairwise=T) #overall p=0.011
TukeyHSD(dis.bray)
#population
dis.braypop <- betadisper(bray,matrixform.a$population, bias.adjust=T)
permutest(dis.braypop, pairwise=T) #overall p=0.011
TukeyHSD(dis.braypop)
#soil
dis.braysoil <- betadisper(bray,matrixform.a$soil, bias.adjust = T)
permutest(dis.braysoil,pairwise=T)
TukeyHSD(dis.braysoil)


#betadisper for trt by population as group 
#interx of factors (since were both sig individually)
newcol <- interaction(matrixform.a$trt,matrixform.a$population)
matrixform.a$trtXpop <- newcol
#dist matrix 
bray <- vegdist(matrixform.a[7:15], method="bray",bias.adjust=T)
#trt by pop 
dis.braytXpop <- betadisper(bray,matrixform.a$trtXpop)
permutest(dis.braytXpop,pairwise=T)

#beta dispersion on all treatment interactions 
trtXpopXsoil <- interaction(matrixform.a$trt,matrixform.a$population,matrixform.a$soil)
matrixform.a$trtXpopXsoil <- trtXpopXsoil
dis.bray3way <- betadisper(bray,matrixform.a$trtXpopXsoil)
permutest(dis.bray3way,pairwise=T)

#warming by soil region 
trtXsoil <- interaction(matrixform.a$trt,matrixform.a$soil)
matrixform.a$trtXsoil <- trtXsoil
dis.braytXs <- betadisper(bray,matrixform.a$trtXsoil)
permutest(dis.braytXs,pairwise=T)

#pop by soil 
popXsoil <- interaction(matrixform.a$population,matrixform.a$soil)
matrixform.a$popXsoil <- popXsoil
dis.braypXs <- betadisper(bray,matrixform.a$popXsoil)
permutest(dis.braypXs,pairwise=T)


#example plot for reduced dispersion across seed provenance 
dis.bray
dis.bray$distances
interaction.short <- data.frame(trt=matrixform.a$trt,population=matrixform.a$population,trtXpop=matrixform.a$trtXpop)
interaction.short$median=dis.bray$distances
melt.intx <- melt(interaction.short,id=c("trt","population","trtXpop"),measure="median")
cast.intx1 <- dcast(melt.intx,trt+population+trtXpop~.,fun=mean)
names(cast.intx1) <- c("trt","population", "trtXpop", "median.disp")
cast.intx2 <- dcast(melt.intx, trt+population+ trtXpop~.,fun=sd)
names(cast.intx2) <- c("trt","population","trtXpop", "sd.disp")
cast.intx3 <- dcast(melt.intx, trt+population+trtXpop~., fun=length)
names(cast.intx3) <- c("trt", "population", "trtXpop", "n")
cast.intx4 <- merge(cast.intx1, cast.intx2, by=c("trt","population","trtXpop"))
cast.intx <- merge(cast.intx4, cast.intx3, by=c("trt","population","trtXpop"))
cast.intx$se <- cast.intx$sd/sqrt(cast.intx$n)
#plot with seed provenance panels 
qplot(trtXpop, median.disp, data=cast.intx, ylim=c(0, 0.5), fill=trt, shape=population, size=I(2.5)) +
  facet_grid(~population, scales="free") +
  geom_errorbar(aes(ymin=median.disp-se, ymax=median.disp+se, width=0.1)) +
  theme(legend.position="none", panel.background=element_rect()) +
  scale_shape_manual(values=c(21,22,24,21,22,24))+
  scale_fill_manual(values=c("white", "black"))



###########--------COMMUNITY---------############
#new community dataframe for simpsons and total productivity 
com <- matrixform.a[,c("heating.ring", "trt", "ring", "population", "soil")]
#calculate diversity for each averaged mesocosm 
com$simpsons <- diversity(matrixform.a[,c("CORLAN", "DALPUR", "ECHPUR", "LESCAP", "RUDHIR", 
                                          "RATPIN", "BOUCUR", "ANDGER", "SCHSCO")], 
                          index="simpson")
#calculate total producitivty for each averaged mesocosm 
com$prod <- apply(matrixform.a[,c("CORLAN", "DALPUR", "ECHPUR", "LESCAP", "RUDHIR", 
                                  "RATPIN", "BOUCUR", "ANDGER", "SCHSCO")], 
                  MARGIN=1, FUN=sum)

#test for normality for diversity and producivity 
com_sw <- mvn(com[7:8], univariateTest = "SW", desc=TRUE)
com_cvm <- mvn(com[7:8], univariateTest = "CVM", desc=TRUE)
com_lil <- mvn(com[7:8], univariateTest = "Lillie", desc=TRUE)
com_sf <- mvn(com[7:8], univariateTest = "SF", desc=TRUE)
com_ad <- mvn(com[7:8], univariateTest = "AD", desc=TRUE)

com_results <- data.frame(com_sw$univariateNormality, com_cvm$univariateNormality,com_lil$univariateNormality,
                          com_sf$univariateNormality,com_ad$univariateNormality)

qq <- mvn(com[7:8], mvnTest = "mardia", univariatePlot = "qqplot")
hist <- mvn(com[7:8], mvnTest = "mardia", univariatePlot = "histogram")

#Mixed effects models with ring/trt as random intercept (meaning that intercept varies among rings and trt, 
#but within each ring (8 levels total bc 8 rings (4 in each trt)))
#Analysis of Variance (Anova() with type III )
#simpon's diversity 
simp <- lmer(simpsons ~ trt*population*soil + (1|ring/trt), data=com); Anova(simp, type="3", test.statistic="F")
#removing NS interx 
simp <- lmer(simpsons ~ trt+ population+ soil + (1|ring/trt), data=com); Anova(simp, type="3", test.statistic="F")

#total community productivity 
prod <- lmer(prod~trt*population*soil + (1|ring/trt), data=com); Anova(prod, type="3",test.statistic="F")
#removing NS interx
prod <- lmer(prod~trt+population+soil + (1|ring/trt), data=com); Anova(prod, type="3",test.statistic="F")

prod.em <- emmeans(prod,~trt|population)
pwpp(prod.em)
prod.em <- emmeans(prod,~population|trt)
pwpp(prod.em)

#Beta dispersion
com$median <- dis.braysoil$distances
beta <- lmer(median~trt*population*soil + (1|ring/trt), data=com); Anova(beta, type="3",test.statistic="F")

beta <- lmer(median~trt+population+soil + (1|ring/trt), data=com); Anova(beta, type="3",test.statistic="F")



#######--------UNIVARIATE aka indiviudal species responses---------#############
#test for univariate normality in spp dataset 
# since our data does not follow a mvn distribution, we cannot necesarily assume that all spp dist are not normally distributed 
result_mard <- mvn(spp, mvnTest = "mardia", univariatePlot = "qqplot")

#histograms
result_mard <- mvn(spp, mvnTest = "mardia", univariatePlot = "histogram")

#quantitative assessment 
result_sw <- mvn(spp, univariateTest = "SW", desc=TRUE)
result_cvm <- mvn(spp, univariateTest = "CVM", desc=TRUE)
result_lil <- mvn(spp, univariateTest = "Lillie", desc=TRUE)
result_sf <- mvn(spp, univariateTest = "SF", desc=TRUE)
result_ad <- mvn(spp, univariateTest = "AD", desc=TRUE)

uni_results <- data.frame(result_sw$univariateNormality, result_cvm$univariateNormality,result_lil$univariateNormality,
                          result_sf$univariateNormality,result_ad$univariateNormality)


#Set up dataframe for species models 
#averaging across pseudoreplicates within a ring/trt (ignoring site)
m.raw.a <- dcast(m.raw, heating.ring+trt+ring+population+soil+variable~., mean) 
names(m.raw.a)[7] <- "biomass"
m.raw.a$logbiomass <- log(m.raw.a$biomass+0.001)
names(m.raw.a)[8] <- "logbiomass"
m.raw.a$species <- m.raw.a$variable
levels(m.raw.a$species) <- c("Coreopsis", "Dalea", "Echinacea", "Lespedeza", "Rudbeckia", 
                             "Ratibida", "Bouteloua", "Andropogon", "Schizachyrium", "Koeleria")

###Species with main effect of warming and/or warming:seed provenance 
#Bouteloua curtipendula 
boucur <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="BOUCUR",]); Anova(boucur, type="3",test.statistic="F")
boucur <- lmer(logbiomass~trt+population+soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="BOUCUR",]); Anova(boucur, type="3", test.statistic="F")
summary(boucur)
sjPlot::plot_model(boucur)

#post hoc 
emmeans(boucur, pairwise ~ trt|population)
emmeans(boucur, pairwise ~ population|trt)
CLD(boucur.emm) 
pwpp(boucur.emm)
pairs(boucur.emm)

##Ratibida pinnata 
ratpin <- lmer(logbiomass~trt*population*soil + (1|ring/trt),
            data=m.raw.a[m.raw.a$variable=="RATPIN",]); Anova(ratpin, type="3",test.statistic="F") 
ratpin <- lmer(logbiomass~trt+population+soil + trt:population +(1|ring/trt), 
               data=m.raw.a[m.raw.a$variable=="RATPIN",]); Anova(ratpin,type="3",test.statistic="F")
summary(ratpin)
sjPlot::plot_model(ratpin)

#post hoc
cld(emmeans(ratpin, list(pairwise~population*trt), adjust = c("Tukey")))
emmeans(ratpin, pairwise ~ trt|population)
emmeans(ratpin, pairwise ~ population|trt)
CLD(ratpin.emm) 
pwpp(ratpin.emm)


###species with main effect of population
##Coreopsis laceolata
corlan <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="CORLAN",]); Anova(corlan,type="3",test.statistic="F")
corlan <- lmer(logbiomass~trt+population+soil + (1|ring/trt), 	
            data=m.raw.a[m.raw.a$variable=="CORLAN",]); Anova(corlan,type="3",test.statistic="F")
summary(corlan)
sjPlot::plot_model(corlan)

#post hoc
cld(emmeans(corlan, list(pairwise~population), adjust = c("Tukey")))
emmeans(corlan, pairwise ~ trt|population)
emmeans(corlan, pairwise ~ population|trt)


##Echinacea purpurea
echpur <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="ECHPUR",]); Anova(echpur, type="3", test.statistic="F")
echpur <- lmer(logbiomass~trt+population+soil+ (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="ECHPUR",]); Anova(echpur, type="3",test.statistic="F")
summary(echpur)
sjPlot::plot_model(echpur)

#post hoc
emmeans(echpur, pairwise ~ trt|population)
emmeans(echpur, pairwise ~ population|trt)


##Rudbeckia hirta
rudhir <- lmer(logbiomass~trt*population*soil + (1|ring/trt),
            data=m.raw.a[m.raw.a$variable=="RUDHIR",]); Anova(rudhir, type="3", test.statistic="F")
rudhir <- lmer(logbiomass~trt+population+soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="RUDHIR",]); Anova(rudhir, type="3", test.statistic="F")
sjPlot::plot_model(rudhir)

#post hoc
cld(emmeans(rudhir, list(pairwise~population*soil), adjust = c("Tukey")))


##Schizachyrium scoparium
schsco <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 	
            data=m.raw.a[m.raw.a$variable=="SCHSCO",]); Anova(schsco, type="3", test.statistic="F")
schsco <- lmer(logbiomass~trt+population+soil+trt:population+trt:soil+trt:population:soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="SCHSCO",]); Anova(schsco, type="3", test.statistic="F")
schsco <- lmer(logbiomass~trt+population+soil+trt:population+trt:soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="SCHSCO",]); Anova(schsco, type="3", test.statistic="F")
schsco <- lmer(logbiomass~trt+population+soil+(1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="SCHSCO",]); Anova(schsco, type="3", test.statistic="F")
summary(schsco)
sjPlot::plot_model(schsco)

#post hoc
emmeans(schsco, pairwise ~ trt|population)
emmeans(schsco, pairwise ~ population|trt)

schsco.emm <- emmeans(schsco,~trt|population)
CLD(schsco.emm) 
pwpp(schsco.emm)


###Species with soil main effects
##Lespedeza capitata 
lescap <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 	
            data=m.raw.a[m.raw.a$variable=="LESCAP",]); Anova(lescap, type="3",test.statistic="F")
lescap <- lmer(logbiomass~trt+population+soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="LESCAP",]); Anova(lescap, type="3",test.statistic="F")
summary(lescap)
sjPlot::plot_model(lescap)

#post hoc
emmeans(lescap, pairwise ~ trt|population)
emmeans(lescap, pairwise ~ population|trt)

lescap.emm <- emmeans(lescap,~soil)
pwpp(lescap.emm)


##Dalea purpurea
dalpur <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="DALPUR",]); Anova(dalpur,type="3",test.statistic="F")
dalpur <- lmer(logbiomass~trt+population+soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="DALPUR",]);Anova(dalpur, type="3", test.statistic="F")

summary(dalpur)
sjPlot::plot_model(dalpur)

#post hoc
CLD(emmeans(dalpur, list(pairwise~population*trt), adjust = c("Tukey")))
dalpur.emm <- emmeans(dalpur,~population|trt)
CLD(dalpur.emm) 
emmeans(dalpur, pairwise ~ trt|population)
cld(emmeans(dalpur, list(pairwise~soil), adjust = c("Tukey"))) 



###Species with population:soil interaction
##Andropogon gerardii 
andger <- lmer(logbiomass~trt*population*soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="ANDGER",]); Anova(andger,type="3", test.statistic="F")
andger <- lmer(logbiomass~trt+population+soil + population:soil + (1|ring/trt), 
            data=m.raw.a[m.raw.a$variable=="ANDGER",]); Anova(andger, type="3",test.statistic="F")
sjPlot::plot_model(andger)

#post hoc
cld(emmeans(andger, list(pairwise~population*soil), adjust = c("Tukey")))

andger.emm <- emmeans(andger,~population|trt)
pwpp(andger.emm)

andger.emm <- emmeans(andger,~population|soil)
pwpp(andger.emm)

andger.emm <- emmeans(andger,~soil|population)
pwpp(andger.emm)
