library(tidyverse)
library(ordinal)
data <- read_delim("C:/Surfdrive/Projecten/CEVA/IBVexperiment/DailyClinicalSignsPerAnimal.csv",";")
names(data)
#transform data to long form
longform = data%>%reshape2::melt(id.vars= c("Group","Subgroup",  "Tag no", "Sample Type", "Vaccinated","Challenged" ), 
                                 variable.name ="dpch",
                                 value.name = "score")

#remove spaces in variable names
names(longform)<- gsub(" ","",names(longform))
longform$score <- factor(longform$score)
longform$Vaccinated <- factor(longform$Vaccinated)
longform$Challenged <- factor(longform$Challenged)
longform$Tagno <-factor(longform$Tagno)

#ordinal regression ####
#regression exclude dpch = 1 dpch because it causes singularities
fitclmm<- clmm(score ~dpch +Vaccinated+Challenged +(1|Tagno), data = longform%>%filter(dpch != "1 dpch"))
drop1(fitclmm,test = "Chisq")
fitclmm<- clmm(score ~dpch +Vaccinated +(1|Tagno), data = longform%>%filter(dpch != "1 dpch"))
drop1(fitclmm,test = "Chisq")
summary(fitclmm)
#too little data for ordinal regression


#Use continuity and dichotomize at each level ####
library(lme4)
for(dl in c(2:6)){
   print(paste("Cutt-off = score ", dl))
   longform$dich <- as.numeric(longform$score) >= dl
   fit.dich.vac<- glmer(dich ~ Vaccinated+ dpch  + (1|Tagno), data = longform%>%filter(dpch != "1 dpch"), family = binomial(link = "logit"))
   fit.dich.chal.vac<- glmer(dich ~ Vaccinated +dpch+ Challenged + (1|Tagno), data = longform%>%filter(dpch != "1 dpch"), family = binomial(link = "logit"))
   #print(summary(fit.dich))
   print(drop1(fit.dich.vac,test = c( "Chisq")))
   print(drop1(fit.dich.chal.vac,test = c( "Chisq")))
 }



