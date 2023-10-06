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

##load data for Ct values
data.file = "./data/U-P-054_postChalL_DMV_CopyNo_Merge.csv"
dataIBV = read.csv(data.file, sep = ";")

#remove empty line
dataIBV = dataIBV[!is.na(dataIBV$group),]

#reclassify groups
dataIBV$Vaccinated = c("Yes","No")[str_detect(dataIBV$group, pattern = "no vac")+1]
dataIBV$Sample = "ON swab"
dataIBV$Strain = "DMV 1639"
dataIBV$Challenge = c("direct","Contact")[(str_detect(dataIBV$group, pattern = "con"))+1]
dataIBV$Group = paste0(dataIBV$Strain,"_",str_remove(dataIBV$group,"-con"),"_",dataIBV$Rep)
dataIBV$bird.id =   dataIBV$Tag.No
rm(sir.data)
sir.data = reform.SIR.data(dataIBV,
                sampleday.vars = names(dataIBV)[4:17],
                inf.rule = 2,
                rec.rule =2,
                 cut.off = 0.0,
                cut.off.dir = "g",
                min.positive = 2,
                model = "SIR",
                exclude.seeder = T,
                SIR.state = F)
sir.data[[2]][order(sir.data[[2]]$bird.id ),]


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
sir.data[[1]]$replication = substrRight(sir.data[[1]]$Group,1)

#check
dataIBV[dataIBV$bird.id == 983,]
sir.data[[2]][sir.data[[2]]$bird.id == 983,]

#fill in 0 for NA's
dataIBVnona <- dataIBV
dataIBVnona[is.na(dataIBVnona)]<-0.0000
rm(sir.data.nona)
sir.data.nona =  reform.SIR.data(dataIBVnona,
                                 sampleday.vars = names(dataIBVnona)[4:17],
                                 cut.off = 0.0,
                                 cut.off.dir = "g",
                                 model = "SIR",
                                 exclude.seeder = T,
                                 SIR.state = F)

sir.data.nona[[1]]$replication = substrRight(sir.data.nona[[1]]$Group,1)

