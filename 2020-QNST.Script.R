##############################Queen nutrient storage timing experiment 2020#######################################

###Set Working Directory

####all packages used here####
require("ggpubr")
require("MuMIn")
require(multcomp)
require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)
library("survival")
library("survminer")
require(RColorBrewer)
require(car)
require(readr)
library(tidyr)
library(lme4)
library(lmerTest)
library(lmtest)
require(outliers)


#####Color palette for QNST figures#####
BluesPalette <- c("#252525","#a8ddb5","#d8b365","#3182BD","#bcbddc")
#Making a palette with black for CTL, then 4 for NS, 4 for PSL:
BluOranPalette <- c("#252525","#a8ddb5","#d8b365","#3182BD","#08519C","##dd1c77","#a6bddb","#E6550D","#A63603")
#Making a palette with black for CTL, 3  for NS, 4 for PS:
BluOran8Palette <- c("#252525","#a8ddb5","#d8b365","#3182BD","#dd1c77","#a6bddb","#E6550D","#31a354")
#Making a palette with black, 2 (NS3d, NS6d), 4 for PS:
BlOrReducedPal<-c("#252525","#a8ddb5","#d8b365","#dd1c77","#a6bddb","#E6550D","#31a354")
#Making a pallete with black for CTL + 4 for PS treatments
OranPalette<- c("#252525","#dd1c77","#a6bddb","#E6550D","#31a354") 


####Survival data####
#this dataset has 5 timepoint observations for each individual bee (0, 3, 6, 9, 12)
Surv_cons <- read.csv("Surv_consol_plus_QNST.csv")
print(Surv_cons, 1:10)
Surv_cons$SampleID = as.factor(Surv_cons$SampleID)
Surv_cons$Treatment = as.factor(Surv_cons$Treatment)
Surv_cons$Colony = as.factor(Surv_cons$Colony)
Surv_cons$Died = as.logical(Surv_cons$died)
Surv_cons$timepoint = as.numeric(Surv_cons$timepoint)
Surv_cons$WingAvgmm = as.numeric(Surv_cons$WingAvgmm)
Surv_cons$StartWt = as.numeric(Surv_cons$StartWt)
##reordering Colony
Surv_cons <- Surv_cons[with(Surv_cons, order(Colony, Treatment)),]
head(Surv_cons)
summary(Surv_cons$Colony)
#Need to drop colonies with <3 bees represented (which shows up as <15 colony counts)
#make a list of colonies with >10 occurrences in the dataset, then keep as a new dataframe
keep <- names(table(Surv_cons$Colony))[table(Surv_cons$Colony)>10]
SC2<-filter(Surv_cons,Colony %in% keep)
#this drops colonies BS02,BS03,MP3,QP1,WP7 because they only have 1-2 representative bees in the dataset
#colonies kept = ABS05,BS01,DR2,MP10,QP2,QP3,WP8,WP9
summary(SC2$Colony)
summary(SC2$Treatment)
#create the response variable
#Death column - whether or not bees died before end of expt
S <- Surv(time = SC2$timepoint, event=SC2$Died)

#building the model using cox PH regression
#the exp (coef) variable is the Hazard ratio
#HR=1: no effect; HR<1: reduction in hazard; HR>1: increase in hazard
#Using Gemma Baron's Survival Cox Regression code 
surv.cox<-coxph(S~1, data = SC2, robust = TRUE) #added robust=TRUE to get robust standard errors & CIs
surv.cox1<-coxph(S~Treatment, data=SC2, robust = TRUE)
surv.cox2<-coxph(S~Colony, data=SC2, robust = TRUE)
surv.cox3<-coxph(S~StartWt, data = SC2, robust = TRUE)
surv.cox4<-coxph(S~WingAvgmm, data = SC2, robust = TRUE)
surv.cox5<-coxph(S~Treatment+Colony, data=SC2, robust = TRUE)
surv.cox6<-coxph(S~Treatment*Colony, data=SC2, robust = TRUE) 
surv.cox7<-coxph(S~Treatment+StartWt, data=SC2, robust = TRUE)
surv.cox8<-coxph(S~Treatment*StartWt, data=SC2, robust = TRUE)
surv.cox9<-coxph(S~Treatment+WingAvgmm, data=SC2, robust = TRUE)
surv.cox10<-coxph(S~Treatment*WingAvgmm, data=SC2, robust = TRUE)
surv.cox11<-coxph(S~Treatment+Colony+StartWt, data=SC2, robust = TRUE)
surv.cox12<-coxph(S~Treatment*Colony+StartWt, data=SC2, robust = TRUE)
surv.cox13<-coxph(S~Treatment+Colony*StartWt, data=SC2, robust = TRUE)
surv.cox14<-coxph(S~Treatment*StartWt+Colony, data=SC2, robust = TRUE)
surv.cox15<-coxph(S~Treatment+StartWt+WingAvgmm, data=SC2, robust = TRUE)
surv.cox16<-coxph(S~Treatment*StartWt+WingAvgmm, data=SC2, robust = TRUE)
surv.cox17<-coxph(S~Treatment+StartWt*WingAvgmm, data=SC2, robust = TRUE)
surv.cox18<-coxph(S~Treatment*WingAvgmm+StartWt, data=SC2, robust = TRUE)
surv.cox19<-coxph(S~Treatment+WingAvgmm+StartWt+Colony, data=SC2, robust = TRUE) 
surv.cox20<-coxph(S~Treatment*WingAvgmm+StartWt+Colony, data=SC2, robust = TRUE) 
surv.cox21<-coxph(S~Treatment+WingAvgmm*StartWt+Colony, data=SC2, robust = TRUE) 
surv.cox22<-coxph(S~Treatment+WingAvgmm+StartWt*Colony, data=SC2, robust = TRUE) 
surv.cox23<-coxph(S~Treatment+WingAvgmm*StartWt+Colony, data=SC2, robust = TRUE) 
surv.cox24<-coxph(S~Treatment*StartWt+Colony+WingAvgmm, data=SC2, robust = TRUE) 
surv.cox25<-coxph(S~Treatment*Colony+StartWt+WingAvgmm, data=SC2, robust = TRUE) 
surv.cox26<-coxph(S~Treatment*WingAvgmm+Colony*StartWt, data=SC2, robust = TRUE) 
#model selection
model.sel(surv.cox,surv.cox1,surv.cox2,surv.cox4,surv.cox5, surv.cox6,surv.cox7,surv.cox8,surv.cox9,surv.cox10,surv.cox11,surv.cox12,
          surv.cox13, surv.cox14, surv.cox15, surv.cox16, surv.cox17, surv.cox18, surv.cox19, surv.cox20, surv.cox21, surv.cox22, surv.cox23, surv.cox24, surv.cox25, surv.cox26 )
summary(surv.cox17)
test.ph <- cox.zph(surv.cox17) #checking for assumptions of proportional hazards
test.ph
Anova(surv.cox17)


####FIGURE - Survival####
#need data that have proportion of bees in each treatment still alive at each timepoint
SurvSummProp3 <- read_csv("SurvSummProp3.csv")
print(SurvSummProp3)
ggplot(SurvSummProp3, aes(x=Days, y=PropAlive, group=Treatment))+geom_point()+ geom_line()+ xlab("Days") + ylab("Proportion Survived")
ggplot(SurvSummProp3, aes(x=Days,y=PropAlive, group=Treatment, color=Treatment))+ geom_step() + xlab("Days")
SurvSP3<- SurvSummProp3[-c(26:45),] #no PS bees - just by excluding rows.
#reorder Treatment to have NS12d last
SurvSP3$Treatment <- factor(SurvSP3$Treatment,levels = c("CTL","NS3d","NS6d","NS9d","NS12d"))
#figure for survival curve
sur <- ggplot(SurvSP3, aes(x=Days,y=PropAlive, group=Treatment, color=Treatment))+ geom_step(size=1) + xlab("Days") + ylab("Proportion of queens survived") +
  theme_classic()+scale_x_continuous(breaks = c(0,3,6,9,12))+scale_color_manual(values = BluesPalette, breaks = c("CTL","NS3d","NS6d","NS9d","NS12d"), labels=c("Control","NS3d","NS6d","NS9d","NS12d"))+
  theme(text = element_text(family = "Arial", size = 12))
sur

#figure 3A and 3B to export 
fig1 <- ggarrange(
  sur,ncol = 1, nrow =  1,
  widths = c(4, 4)
)
fig1

####Incremental Weight change####
WtIncrem2 <- read.csv("WtIncQNSTnew.csv")
head(WtIncrem2)
WtIncrem2$SampleID = as.factor(WtIncrem2$SampleID)
WtIncrem2$Colony = as.factor(WtIncrem2$Colony)
WtIncrem2$Treatment = as.factor(WtIncrem2$Treatment)
WtIncrem2$timepoint = as.factor(WtIncrem2$timepoint)
WtIncrem2$weight = as.numeric(WtIncrem2$weight)
WtIncrem2$WingAvgmm = as.numeric(WtIncrem2$WingAvgmm)

#check for normality
ggdensity(WtIncrem2$weight)
ggqqplot(WtIncrem2$weight)
shapiro.test(WtIncrem2$weight) #weight data are not normally distributed. p=0.016
#need to figure out transformations, how we handle the data. More possible code below#
#Need to drop colonies with <3 bees represented 
#make a list of colonies with >10 occurrences in the dataset (>2 bees), then keep as a new dataframe
summary(WtIncrem2$Colony)
keep3 <-names(table(WtIncrem2$Colony))[table(WtIncrem2$Colony)>10]
WtIncRedx<-filter(WtIncrem2,Colony %in% keep3)
#same colonies kept/dropped as in survival dataset
#this drops colonies BS01,BS02,BS03,MP3,QP1,WP7 because they only have 1-2 representative bees in the dataset
#colonies kept = BS05,DR2,MP10,QP2,QP3,WP8
summary(WtIncRedx$Colony)

#separate in two dataset one for Pollen and another for Nectar
WtIncRedxN = filter(WtIncRedx, !(Treatment %in% c("PS12d", "PS3d", "PS6d", "PS9d"))) #nectar and control samples
WtIncRedxP = filter(WtIncRedx, !(Treatment %in% c("NS12d", "NS3d", "NS6d", "NS9d"))) #pollen and control samples

#redo analysis with reduced data for Nectar
WtIncRedxN<-WtIncRedxN %>% 
  filter(Treatment != "NS9d") #removed 1 outlier, NS9d has few ind. 
head(WtIncRedxN)
dim(WtIncRedxN) 


lm.nullN <- lmer(weight ~ 1 + (1|SampleID) + (1|Colony), data = WtIncRedxN, REML=FALSE)
lm.nullN0 <- lmer(weight ~ 1 + (1|Colony), data = WtIncRedxN, REML=FALSE)
lm.nullN1 <- lmer(weight ~ 1 + (1|SampleID), data = WtIncRedxN, REML=FALSE)
lm.wtIN1 <- lmer(weight ~ Treatment + (1|SampleID) + (1|Colony), data = WtIncRedxN, REML=FALSE)
lm.wtIN2 <- lmer(weight ~ timepoint + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML=FALSE)
lm.wtIN3 <- lmer(weight ~ WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML=FALSE)
lm.wtIN5 <- lmer(weight ~ Treatment + timepoint + (1|SampleID) + (1|Colony), data = WtIncRedxN, REML = FALSE)
lm.wtIN6 <- lmer(weight ~ Treatment * timepoint + (1|SampleID) + (1|Colony), data = WtIncRedxN, REML=FALSE)
lm.wtIN7 <- lmer(weight ~ Treatment + WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML = FALSE)
lm.wtIN8 <- lmer(weight ~ Treatment * WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML = FALSE)
lm.wtIN9 <- lmer(weight ~ Treatment + timepoint + WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML=FALSE)
lm.wtIN10 <- lmer(weight ~ Treatment * timepoint + WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML=FALSE)
lm.wtIN11 <- lmer(weight ~ Treatment + timepoint * WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxN, REML=FALSE)
#model selection
model.sel(lm.nullN,lm.nullN1, lm.nullN0, lm.wtIN1, lm.wtIN2, lm.wtIN3, lm.wtIN5,lm.wtIN6, lm.wtIN7,lm.wtIN8, lm.wtIN9, lm.wtIN10, lm.wtIN11)
summary(lm.wtIN6) #best fit model
Anova(lm.wtIN6)
lrtest(lm.wtIN6,lm.nullN)
summary(glht(lm.wtIN6, linfct=mcp(Treatment="Tukey")))
summary(glht(lm.wtIN6, linfct=mcp(timepoint="Tukey")))

#redo analysis with reduced data for Pollen
lm.nullP <- lmer(weight ~ 1 + (1|SampleID) + (1|Colony), data = WtIncRedxP, REML=FALSE)
lm.nullP0 <- lmer(weight ~ 1 + (1|Colony), data = WtIncRedxP, REML=FALSE)
lm.nullP1 <- lmer(weight ~ 1 + (1|SampleID), data = WtIncRedxP, REML=FALSE)
lm.wtIP1 <- lmer(weight ~ Treatment + (1|SampleID) + (1|Colony), data = WtIncRedxP, REML=FALSE)
lm.wtIP2 <- lmer(weight ~ timepoint + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML=FALSE)
lm.wtIP3 <- lmer(weight ~ WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML=FALSE)
lm.wtIP5 <- lmer(weight ~ Treatment + timepoint + (1|SampleID) + (1|Colony), data = WtIncRedxP, REML = FALSE)
lm.wtIP6 <- lmer(weight ~ Treatment * timepoint + (1|SampleID) + (1|Colony), data = WtIncRedxP, REML=FALSE)
lm.wtIP7 <- lmer(weight ~ Treatment + WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML = FALSE)
lm.wtIP8 <- lmer(weight ~ Treatment * WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML = FALSE)
lm.wtIP9 <- lmer(weight ~ Treatment + timepoint + WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML=FALSE)
lm.wtIP10 <- lmer(weight ~ Treatment * timepoint + WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML=FALSE)
lm.wtIP11 <- lmer(weight ~ Treatment + timepoint * WingAvgmm + (1|SampleID) + (1|Colony), data=WtIncRedxP, REML=FALSE)
#model selection
model.sel(lm.nullP,lm.nullP1, lm.nullP0, lm.wtIP1, lm.wtIP2, lm.wtIP3, lm.wtIP5,lm.wtIP6, lm.wtIP7,lm.wtIN8, lm.wtIP9, lm.wtIP10, lm.wtIP11)
summary(lm.wtIP6) #best fit model
Anova(lm.wtIP6)
lrtest(lm.wtIP6,lm.nullP)
summary(glht(lm.wtIP6, linfct=mcp(Treatment="Tukey")))
summary(glht(lm.wtIP6, linfct=mcp(timepoint="Tukey")))


####FIGURE - Weight change####

#Weight change
WtTotalChange <- read_csv("WeightChange1.csv")
print(WtTotalChange)
WtTotalChange$SampleID = as.factor(WtTotalChange$SampleID)
WtTotalChange$Colony = as.factor(WtTotalChange$Colony)
WtTotalChange$Treatment = as.factor(WtTotalChange$Treatment)
WtTotalChange$StartWeight = as.numeric(WtTotalChange$`StartWeight(g)`)
WtTotalChange$Weight_d3 = as.numeric(WtTotalChange$Weight_d3)
WtTotalChange$Weight_d6 = as.numeric(WtTotalChange$Weight_d6)
WtTotalChange$Weight_d9 = as.numeric(WtTotalChange$Weight_d9)
WtTotalChange$EndWeight = as.numeric(WtTotalChange$`EndWeight(g)_d12`)
WtTotalChange$WingAvgmm = as.numeric(WtTotalChange$WingAvgmm)
WtTotalChange$WtCh = WtTotalChange$EndWeight-WtTotalChange$StartWeight #total weight change for a bee
#removing empty NA rows
WeightCh<-WtTotalChange[!(rowSums(is.na(WtTotalChange))==NCOL(WtTotalChange)),] 
print(WeightCh)
WeightCh$Colony = as.factor(WeightCh$Colony)
#Need to drop colonies with <3 bees represented 
#make a list of colonies with >2 occurrences in the dataset, then keep as a new dataframe
keep2 <-names(table(WeightCh$Colony))[table(WeightCh$Colony)>2]
TotWtRedx<-filter(WeightCh,Colony %in% keep2)
#same colonies kept/dropped as in survival dataset
#this drops colonies BS02,BS03,MP3,QP1,WP7 because they only have 1-2 representative bees in the dataset
#colonies kept = ABS05,BS01,DR2,MP10,QP2,QP3,WP8,WP9
#check for normality
ggdensity(TotWtRedx$WtCh)
ggqqplot(TotWtRedx$WtCh)
shapiro.test(TotWtRedx$WtCh) #weight data are NOT normally distributed. p=0.003913
kruskal.test(WtCh~ Colony, data = TotWtRedx)
#boxplots just to check the data
ggplot(TotWtRedx, aes(Colony, WtCh)) + geom_boxplot()
ggplot(TotWtRedx, aes(Treatment, WtCh)) + geom_boxplot()
sum(TotWtRedx$Colony == 'WP9', na.rm=TRUE) #counting number of WP9 Colony indivs.
(NewTotWtRx <- TotWtRedx[TotWtRedx$Colony != "WP9", ]) #Removing WP9 samples because weights seem very low ##USE THESE DATA
#boxplots just to check the data again
ggdensity(NewTotWtRx$WtCh)
shapiro.test(NewTotWtRx$WtCh) #weight data ARE normally distributed. p=0.3677
ggplot(NewTotWtRx, aes(Colony, WtCh)) + geom_boxplot()
ggplot(NewTotWtRx, aes(Treatment, WtCh)) + geom_boxplot()

#Reconfigure DATA for Incremental Weight Change FIGURE
#using NewTotWtRedx, which excludes low-representative colonies
TFNTWRx<-NewTotWtRx
TFNTWRx$time0 = as.numeric(TFNTWRx$`StartWeight(g)`)
TFNTWRx$time3 = as.numeric(TFNTWRx$Weight_d3)
TFNTWRx$time6 = as.numeric(TFNTWRx$Weight_d6)
TFNTWRx$time9 = as.numeric(TFNTWRx$Weight_d9)
TFNTWRx$time12 = as.numeric(TFNTWRx$`EndWeight(g)_d12`)
WtIncNew <- TFNTWRx %>%
  pivot_longer(c('time0','time3','time6','time9','time12'), names_to = "timepoint", values_to = "beeweight")
WtIncNew2 <- select(WtIncNew,SampleID,Colony,Treatment,timepoint,beeweight)
WtIncNew2A<-WtIncNew2 %>% 
  separate(timepoint, into = c("text", "Timepoint"), sep = 4)
WtIncNew2A
WtIncNew2A$Timepoint=as.numeric(WtIncNew2A$Timepoint)
#consolidate the data into format needed for figure - low rep colonies & WP9 excluded.
WtIncNew4<-WtIncNew2A %>%
  group_by(Treatment, Timepoint) %>%
  summarise_each(funs(
    mean(beeweight, na.rm=T), n2=sum(!is.na(beeweight)), 
    se2=sd(beeweight, na.rm=T)/sqrt(sum(!is.na(beeweight)))))
WtIncNew4

#subsetting for CTL+NS figure
WtIncNS = filter(WtIncNew4, Treatment %in% c("control","NS3d","NS6d","NS9d","NS12d"))
WtIncNS
WtIncNS$mean<-as.numeric(WtIncNS$beeweight_mean)
WtIncNS$SE<-as.numeric(WtIncNS$beeweight_se2)
WtIncNS$Treatment <- factor(WtIncNS$Treatment,levels = c("control","NS3d","NS6d","NS9d","NS12d"), labels = c("Control","NS3d","NS6d","NS9d","NS12d"))
FigNS<-ggplot(WtIncNS, aes(x=Timepoint, y=WtIncNS$mean, color=Treatment, group=Treatment)) + 
  geom_line(size=0.5) + geom_point(aes(shape=Treatment),fill = "white", size=2)+ 
  scale_x_continuous(breaks = c(0,3,6,9,12)) +geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2) +
  theme_bw() + ylab("Mean ± s.e.m. bee weight (g)") + xlab("Days") + scale_color_manual(values = BluesPalette)+
  scale_shape_manual(values = c(16,17,15,8,7))+theme(text = element_text(size = 14))+
  scale_y_continuous(limits = c(0.35,0.65))
FigNS

#subsetting for CTL+PS figure
WtIncPS = filter(WtIncNew4, Treatment %in% c("control","PS3d","PS6d","PS9d","PS12d"))
WtIncPS
WtIncPS$Mean<-as.numeric(WtIncPS$beeweight_mean)
WtIncPS$SE<-as.numeric(WtIncPS$beeweight_se2)
WtIncPS$Treatment <- factor(WtIncPS$Treatment,levels = c("control","PS3d","PS6d","PS9d","PS12d"), labels = c("Control","PS3d","PS6d","PS9d","PS12d"))
FigPS<-ggplot(WtIncPS, aes(x=Timepoint, y=Mean, color=Treatment, group=Treatment)) + 
  geom_line(size=0.5) + geom_point(aes(shape=Treatment),fill = "white", size=3)+ 
  scale_x_continuous(breaks = c(0,3,6,9,12))+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2) +
  theme_bw() + ylab("Mean ± s.e.m. bee weight (g)") + xlab("Days") +scale_color_manual(values = OranPalette) +
  scale_shape_manual(values = c(16,17,15,8,7))+theme(text = element_text(size = 14))+
  scale_y_continuous(limits = c(0.35,0.65))
FigPS


#figure 2A and 2B to export 
fig2 <- ggarrange(
  FigNS, FigPS, ncol = 1, nrow =  2,
  labels = c("A", "B"),
  common.legend = FALSE, legend = "right",
  widths = c(3, 3)
)
fig2

#####Lipids Assasy#####
Lipid <- read.csv("QNSTLipidsConcCalc.csv")
head(Lipid)
dim(Lipid)
Lipid$SampleID = as.factor(Lipid$SampleID)
Lipid$Colony = as.factor(Lipid$Colony)
Lipid$Treatment = as.factor(Lipid$Treatment)
Lipid$Conc = as.numeric(Lipid$TotSampleConc.mg.mL.)
Lipid$AbdDryWt = as.numeric(Lipid$AbdDryWt)
Lipid<-mutate(Lipid, NormLipConc = Conc/(AbdDryWt*1000)) #making a new column that is Conc.of Lipid normalized by bee size (mass). This is NormLipConc
head(Lipid)
dim(Lipid)

ggplot(Lipid, aes(Treatment, Conc)) + geom_boxplot() + theme_bw() + ylab("Lipid concentration (mg/mL)") + theme(legend.position = "none")
library(ggpubr)
ggdensity(Lipid$Conc)
ggqqplot(Lipid$Conc)
shapiro.test(Lipid$Conc) #Lipid data are NOT normally distributed. p=0.01376

#identification of outliers
x <- Lipid$Conc
qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
caps <- quantile(x, probs=c(.05, .95), na.rm = T)
H <- 1.5 * IQR(x, na.rm = T)
x[x < (qnt[1] - H)] <- caps[1]
x[x > (qnt[2] + H)] <- caps[2]
scores(x)  #z-scores => (x-mean)/sd
scores(x, type="chisq")  #chi-sq scores => (x - mean(x))^2/var(x)
scores(Lipid$Conc, type="t", prob=0.95)
Lipid$Outliers<-scores(Lipid$Conc, type="z", prob=0.95) #these give which samples are above the 95% CI
#this removes the TRUE outliers from this dataset
LipRedX<-Lipid %>%
  filter(Outliers == "FALSE")
head(LipRedX) 
dim(LipRedX) #not removed
#with manual cutoffs for outliers (and drops NAs) now with NormLipConc
LipRedX<-Lipid %>% 
  filter(Treatment != "NS9d")
head(LipRedX)
dim(LipRedX) #removed 1 outlier

#check the normality Conc (mg/ml)
ggplot(LipRedX, aes(Treatment, Conc)) + geom_point()
ggdensity(LipRedX$Conc)
ggqqplot(LipRedX$Conc)
shapiro.test(LipRedX$Conc) #NOT normally distributed. p-value = 0.02156
#attempting a log-transformation of the data
LipRedX$Log.Con <-log(LipRedX$Conc)
ggdensity(LipRedX$Log.Con)
ggqqplot(LipRedX$Log.Con)
shapiro.test(LipRedX$Log.Con) #STILL NOT normally distributed. p-value = 5.665e-05

#lmer with Concetration (mg/ml) as response
lnull <- glmer(Conc ~1 + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l1 <- glmer(Conc ~ Treatment + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l2 <- glmer(Conc ~ WingAvgmm + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l3 <- glmer(Conc ~ AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l4 <- glmer(Conc ~ Treatment + WingAvgmm + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l5 <- glmer(Conc ~ Treatment * WingAvgmm + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l6 <- glmer(Conc ~ Treatment + AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l7 <- glmer(Conc ~ Treatment * AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l8 <- glmer(Conc ~ WingAvgmm + AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l9 <- glmer(Conc ~ WingAvgmm * AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l10 <- glmer(Conc ~ Treatment + WingAvgmm + AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l11 <- glmer(Conc ~ Treatment * WingAvgmm + AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l12 <- glmer(Conc ~ Treatment + WingAvgmm * AbdDryWt + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
l13 <- glmer(Conc ~ Treatment * AbdDryWt + WingAvgmm + (1|Colony),data=LipRedX, family=Gamma(link = "inverse"))
#assessing which model is the best fit
model.sel(lnull, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13)
summary(l1)
summary(glht(l1, linfct = mcp(Treatment = "Tukey")))
anova(l1,lnull)
lrtest(l1,lnull)


####FIGURE - Lipids####
LipRedX$Treatment <- factor(LipRedX$Treatment,levels = c("CTL","NS3d","NS6d","PS3d","PS6d","PS9d","PS12d"), labels = c("Control","NS3d","NS6d","PS3d","PS6d","PS9d","PS12d"))
LipFig<-ggplot(LipRedX, aes(Treatment, Conc)) + geom_boxplot(aes(fill=Treatment)) + theme_bw() + ylab("Lipid concentration (mg/mL)") + 
  theme(legend.position = "none",text = element_text(family = "Arial", size = 12)) + scale_fill_manual(values = BlOrReducedPal) +
  geom_bracket(xmin = "PS3d", xmax = c("PS9d"), y.position = c(0.025), label = c("*"), label.size = 5)
LipFig


#####Glycogen Assasy#####
Glyc <- read.csv("QNSTGlycConcCalc.csv")
head(Glyc)
dim(Glyc)
Glyc$SampleID = as.factor(Glyc$SampleID)
Glyc$Colony = as.factor(Glyc$Colony)
Glyc$Treatment = as.factor(Glyc$Treatment)
Glyc$Conc = as.numeric(Glyc$ReTotGlycConc.mg.mL.)
Glyc$WingAvgmm = as.numeric(Glyc$WingAvgmm)
Glyc$AbdDryWt = as.numeric(Glyc$AbdDryWt)
Glyc<-mutate(Glyc, NormGlycConc = Conc/(AbdDryWt*1000)) #making a new column that is Conc.of Glycogen normalized by bee size (mass). This is NormGlycConc

#identification of outliers
x <- Glyc$Conc
qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
caps <- quantile(x, probs=c(.05, .95), na.rm = T)
H <- 1.5 * IQR(x, na.rm = T)
x[x < (qnt[1] - H)] <- caps[1]
x[x > (qnt[2] + H)] <- caps[2]
scores(x)  #z-scores => (x-mean)/sd
scores(x, type="chisq")  #chi-sq scores => (x - mean(x))^2/var(x)
scores(Glyc$Conc, type="t", prob=0.95)
Glyc$Outliers<-scores(Glyc$Conc, type="z", prob=0.95) #these give which samples are above the 95% CI
#this removes the TRUE outliers from this dataset
GlycRedX<-Glyc %>%
  filter(Outliers == "FALSE")
head(GlycRedX) 
dim(GlycRedX) #removed 3 outliers
#with manual cutoffs for outliers (and drops NAs) now with NormGlycConc
GlycRedX<-Glyc %>% 
  filter(NormGlycConc <0.065)
GlycRedX
#if we want to cut out the other two outliers, use this cutoff 
GlycRedX<-Glyc %>%
  filter(NormGlycConc <0.036)  #this drops it down to 52 samples
GlycRedX
head(GlycRedX)
dim(GlycRedX)

#check the normality Conc (mg/ml)
ggplot(GlycRedX, aes(Treatment, Conc)) + geom_point()
ggdensity(GlycRedX$Conc)
ggqqplot(GlycRedX$Conc)
shapiro.test(GlycRedX$Conc) #NOT normally distributed. p-value = 0.0002238
#attempting a log-transformation of the data
GlycRedX$Log.Con <-log(GlycRedX$Conc)
ggdensity(GlycRedX$Log.Con)
ggqqplot(GlycRedX$Log.Con)
shapiro.test(GlycRedX$Log.Con) #normally distributed. p-value = 0.1401

#lmer with Concetration (mg/ml) as response
gnull <- lmer(Log.Con ~ 1 + (1|Colony),data=GlycRedX, REML=FALSE)
g1 <- lmer(Log.Con ~ Treatment + (1|Colony),data=GlycRedX, REML=FALSE)
g2 <- lmer(Log.Con ~ WingAvgmm + (1|Colony),data=GlycRedX, REML=FALSE)
g3 <- lmer(Log.Con ~ AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g4 <- lmer(Log.Con ~ Treatment + WingAvgmm + (1|Colony),data=GlycRedX, REML=FALSE)
g5 <- lmer(Log.Con ~ Treatment * WingAvgmm + (1|Colony),data=GlycRedX, REML=FALSE)
g6 <- lmer(Log.Con ~ Treatment + AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g7 <- lmer(Log.Con ~ Treatment * AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g8 <- lmer(Log.Con ~ WingAvgmm + AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g9 <- lmer(Log.Con ~ WingAvgmm * AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g10 <- lmer(Log.Con ~ Treatment + WingAvgmm + AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g11 <- lmer(Log.Con ~ Treatment * WingAvgmm + AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g12 <- lmer(Log.Con ~ Treatment + WingAvgmm * AbdDryWt + (1|Colony),data=GlycRedX, REML=FALSE)
g13 <- lmer(Log.Con ~ Treatment * AbdDryWt + WingAvgmm + (1|Colony),data=GlycRedX, REML=FALSE)
#assessing which model is the best fit
model.sel(gnull, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13)
summary(g3)
anova(g3,gnull)
summary(glht(g1, linfct = mcp(Treatment = "Tukey")), test = adjusted("holm"))
lrtest(g3,gnull)


####FIGURE - Glycogen####
GlycRedX$Treatment <- factor(GlycRedX$Treatment,levels = c("CTL","NS3d","NS6d","PS3d","PS6d","PS9d","PS12d"), labels = c("Control","NS3d","NS6d","PS3d","PS6d","PS9d","PS12d"))
GlycFig<-ggplot(GlycRedX, aes(Treatment, Conc)) + geom_boxplot(aes(fill=Treatment)) + theme_bw() + ylab("Glycogen concentration (mg/mL)") + 
  theme(legend.position = "none",text = element_text(family = "Arial", size = 12)) + scale_fill_manual(values = BlOrReducedPal)
GlycFig


#figure 3A and 3B to export 
fig3 <- ggarrange(
  LipFig, GlycFig, ncol = 1, nrow =  2,
  labels = c("A", "B"),
  widths = c(4, 4)
)
fig3

