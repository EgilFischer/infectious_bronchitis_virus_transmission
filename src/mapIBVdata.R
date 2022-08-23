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
data.file = "C://Surfdrive//Projecten//CEVA//IBVexperiment//IBVTransmission//data//054 Animal Data Christophe.xlsx"
dataIBV = read.xlsx(data.file
                    , sheet = "Vac (Group 1)")
dataIBV = rbind(dataIBV,read.xlsx(data.file
                    , sheet = "Non-Vac (Group 2)"))
#remove empty line
dataIBV = dataIBV[!is.na(dataIBV$Vaccinated),]

write.csv(dataIBV, file = "C://Surfdrive//Projecten//CEVA//IBVexperiment//IBVTransmission//data//IBVdata.csv")

#reclassify groups
dataIBV$Group = paste0(dataIBV$Strain,"_",dataIBV$Group,"_",dataIBV$REPL)

sir.data = reform.SIR.data(dataIBV,
                sampleday.vars = names(dataIBV)[7:20],
                cut.off = 36.0,
                min.positive = 2,
                model = "SIR",
                exclude.seeder = T,
                SIR.state = F)



substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
sir.data[[1]]$replication = substrRight(sir.data[[1]]$Group,1)

#check
dataIBV[dataIBV$bird.id == 983,]
sir.data[[2]][sir.data[[2]]$bird.id == 983,]

#fill in 40 for NA's
dataIBVnona <- dataIBV
dataIBVnona[is.na(dataIBVnona)]<-40.0000
sir.data.nona =  reform.SIR.data(dataIBVnona,
                                 sampleday.vars = names(dataIBVnona)[7:20],
                                 cut.off = 30.0,
                                 model = "SIR",
                                 exclude.seeder = T,
                                 SIR.state = F)

sir.data.nona[[1]]$replication = substrRight(sir.data.nona[[1]]$Group,1)
