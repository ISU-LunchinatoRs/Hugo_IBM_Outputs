######## LOAD PACKAGES ##########################################################################
library(ggplot2)
library(emmeans)
library(ggResidpanel)
library(lmerTest)
library(MASS)
library(rsq)

################ STATISTICAL MODEL #####################################################################################\
### LOAD AND WRANGLE DATA
data<-read.csv(file="esali_Data.csv", header=TRUE,sep=";")
data$name<-as.factor(data$name) ### Name of cells as factor
summary(data) ### LOOK AT THE DATA

### Change landcover naming
levels(data$cellCover) <- c(levels(data$cellCover),"Degraded","Refuge","Native")
data$cellCover[data$cellCover=='Landcover(1)']<-"Degraded"
data$cellCover[data$cellCover=='Landcover(0)']<-"Refuge"
data$cellCover[data$cellCover=='Landcover(2)']<-"Native"
data$cellCover<-droplevels(data$cellCover)
summary(data)

### Drop Refuge
dataFor<-data[which(data$cellCover!="Refuge"),]
dataFor$cellCover <- factor(dataFor$cellCover)
summary(dataFor)

### Plot Seeds
ggplot(dataFor, aes(x=NonSeeds)) + geom_histogram(binwidth=50)
ggplot(dataFor, aes(x=NatSeeds)) + geom_histogram(binwidth=100)
ggplot(dataFor, aes(x=TotSeeds)) + geom_histogram(binwidth=100)

### LETS DO GLMS #########################################
## Total Seeds
totGLM<-glm(TotSeeds~Disperser*cellCover*Treatment+DistNat+DistRef, dataFor, family="quasipoisson")
summary(totGLM)
resid_panel(totGLM)
rsq(totGLM,adj=TRUE)

testMeans <- emmeans(totGLM, ~ Disperser*cellCover*Treatment)
testMeans
contrast(testMeans,"consec",simple="each",combine=TRUE,adj="mvt")
emmip(testMeans,Disperser~Treatment|cellCover)

# CAN TEST WITH NEGATIVE BINOMIAL
totGLM<-glm.nb(TotSeeds~Disperser*cellCover*Treatment+DistNat+DistRef, dataFor)
summary(totGLM)
rsq(totGLM,adj=TRUE)
plot(totGLM)

# CAN TEST WITH RANDOM EFFECt TO ACCOUNT FOR OVERDISPERSION
totGLM<-glmer(TotSeeds~Disperser*cellCover*Treatment+DistNat+DistRef+(1|name/Simulation), dataFor, family="poisson")
summary(totGLM)
rsq(totGLM,adj=TRUE)
plot(totGLM)

## Native Seeds
natGLM<-glm(NatSeeds~Disperser*cellCover*Treatment+DistNat+DistRef, dataFor, family="quasipoisson")
summary(natGLM)
resid_panel(natGLM)
rsq(natGLM,adj=TRUE)

testMeans <- emmeans(natGLM, ~ Disperser*cellCover*Treatment)
testMeans
contrast(testMeans,"consec",simple="each",combine=TRUE,adj="mvt")
emmip(testMeans,Disperser~Treatment|cellCover)

## NonNative Seeds
nonGLM<-glm(NonSeeds~Disperser*cellCover*Treatment+DistNat+DistRef, dataFor, family="quasipoisson")
summary(nonGLM)
resid_panel(nonGLM)
rsq(nonGLM,adj=TRUE)

testMeans <- emmeans(nonGLM, ~ Disperser*cellCover*Treatment)
contrast(testMeans,"consec",simple="each",combine=TRUE,adj="mvt")
emmip(testMeans,Disperser~Treatment|cellCover)

# IF DONE WITH POISSON, SIGNIFICANT DIFFERENCE IN BORDER SCENARIO
nonGLM<-glm(NonSeeds~Disperser*cellCover*Treatment+DistNat+DistRef, dataFor, family="poisson")
testMeans <- emmeans(nonGLM, ~ Disperser*cellCover*Treatment)
contrast(testMeans,"consec",simple="each",combine=TRUE,adj="mvt")