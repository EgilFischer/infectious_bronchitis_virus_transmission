## First specify the packages of interest
packages = c("tidyverse",
             "openxlsx")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

##load data 
dataIBV = read.xlsx("C://Surfdrive//Projecten//CEVA//IBVexperiment//U-P-014 Animal Data-Christophe_38-40 neg (updated30Sep & 07Oct).xlsx"
                    , sheet = "Vac (Group1)")
dataIBV = rbind(dataIBV,read.xlsx("C://Surfdrive//Projecten//CEVA//IBVexperiment//U-P-014 Animal Data-Christophe_38-40 neg (updated30Sep & 07Oct).xlsx"
                    , sheet = "Non Vac (Group 2)"))
write.csv(dataIBV, file = "C://Surfdrive//Projecten//CEVA//IBVexperiment//IBVdata.csv")

#reclassify groups
dataIBV$Group = paste0(dataIBV$Strain,"_",dataIBV$Group)


sir.data = reform.SIR.data(dataIBV,
                sampleday.vars = names(dataIBV)[7:19],
                cut.off = 39.99,
                SIR.state = F)


fit<- glm(cbind(C,S-C)~Group,
    offset = -log(I/N), 
    family = binomial(link ="cloglog"), 
    data =sir.data[[1]][sir.data[[1]]$I>0,])
summary(fit)
sir.data