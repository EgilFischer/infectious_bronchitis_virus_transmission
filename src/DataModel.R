####################################################################
#       IBV data model                                             #
#       Casts a data set with measurements pre day into a data set #

#       that can be used for parameters estimation                 #
#       State 0 = susceptible, state 1 = infectious, state 2 = recovered
####################################################################

#fill in groups and types of birds when NA
fill.type.na.rename <- function(data){
  new.data <- data;
  records <- length(data[,1]);
  for(i in c(1:records))
  {
      types <- c(ifelse(!is.na(new.data[i,]$Group),new.data[i,]$Group, types[1]),#set group
                 ifelse(!is.na(new.data[i,]$Subgroup),new.data[i,]$Subgroup, types[2]),#set group
                 ifelse(!is.na(new.data[i,]$Role.in.the.study),new.data[i,]$Role.in.the.study, types[3]),#set group
                 ifelse(!is.na(new.data[i,]$bird.id),new.data[i,]$bird.id, types[4]))#set group
      new.data[i,c(1:4)]<- types
        
  }
  colnames(new.data)[c(1:3,5)]<- c("Vaccinated", "Group","ci","Sample")
  new.data$Vaccinated <- ifelse(new.data$Vaccinated =="Vaccinated", "Yes", "No")
  new.data$Challenge <- ifelse(new.data$ci =="seeder", "seeder", "contact")
  return(new.data[!is.na(data$Sample.type),])#remove empty records
}


#uses file with states(positive / negative) already determined manually from Ct value ####
reform.SIR.data <- function(data, 
                            sampleday.vars,
                            exclude.seeder = T,
                            SIR.state = T,
                            ...
){
  reform.data(data,
              sampleday.vars = sampleday.vars,
              exclude.seeder = exclude.seeder,
                  ...)
}

#count infections based on each type of sample separately 
reform.data<- function(data, 
                       model = c("SIS","SIR","SI")[1],#transmission model default SIS model
                       inf.rule =1, #rule with number of consecutive positive samples to be determined infectious
                       rec.rule =1, #rule with number of consecutive negative samples to be determined recovered
                       sampleday.vars,
                       group.by = NULL, 
                       positive.if = c("max","min")[1],#max will use at least 1 positive sample, min all samples positive.
                       cut.off = 0,
                       cut.off.dir = "leq",#direction of cut-off
                       min.positive = 1,
                       exclude.seeder = T,
                       SIR.state = F,#False means that SIR status needs to be determined
                       category = "ON swab"
)
{
  #remove manual classification which was added for the publication
  data <-data[(data$Sample==category&!is.na(data$Sample)),];
  #
  if(!SIR.state){
    #start with a recoding data to binary yes/no
    binary.data <- data; 
    compare.function = if(cut.off.dir == "leq"){function(x){ifelse(x <= cut.off,1,0)}}else
      if(cut.off.dir == "geq"){function(x){ifelse(x >= cut.off,1,0)}}else
        if(cut.off.dir == "l"){function(x){ifelse(x < cut.off,1,0)}}else
          if(cut.off.dir == "g"){function(x){ifelse(x > cut.off,1,0)}}
      
    binary.data[, sampleday.vars] <- lapply(binary.data[, sampleday.vars], 
                                             FUN = compare.function);
    
    #check for minimum of positive samples otherwise make all of them negative
    for(j in c(1:length(binary.data[,1])))
      {
    
       if(sum(binary.data[j, sampleday.vars],na.rm = T)<min.positive)
       {
         binary.data[j, sampleday.vars]<- binary.data[j, sampleday.vars]*0
       }
    }
    
      
    
    #group positive samples based on positive.if statement if samples of the same bird are taken.
    aggregate.data <- aggregate.data.frame(binary.data[,sampleday.vars], 
                                           by = list(Group = binary.data$Group, 
                                                     bird.id = binary.data$bird.id, 
                                                     Vaccinated = binary.data$Vaccinated,
                                                     Challenge = binary.data$Challenge,
                                                     Strain = binary.data$Strain
                                                     )
                                           ,
                                           FUN = positive.if)
    #recode based on transmission model
    for(ani in c(1:length(aggregate.data[,1])))
    {
   
      
      aggregate.data[ani,sampleday.vars]<-SIR.state(as.vector(aggregate.data[ani,sampleday.vars][1,]),
                                                    model = model,
                                                    inf.rule = inf.rule,
                                                    rec.rule = rec.rule, 
                                                    ani_nr = ani)
      
    }
    
    #change the colnames
colnames(aggregate.data)[1:5]<- c("Group","bird.id","Vaccinated","Challenge","Strain")
  }
  else{
    aggregate.data <- data #already as input presented.
  }
  #aggregate to number of infections prior to time step TODO is to put this in a separate function
  id<-NULL; inf <- NULL; rec<- NULL;n <- NULL; cases <- NULL;susceptibles <- NULL;
  for(day in c(1:length(sampleday.vars)))
  {
    
    #determine the number of infectious on this day
    inf <- rbind(inf, data.frame(dpch = sampleday.vars[day],
                                 ndpch = day,
                                 aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                      by = list(aggregate.data$Group, 
                                                                aggregate.data$Vaccinated,
                                                                aggregate.data$Strain), 
                                                      FUN = function(x){sum(x==1,na.rm = TRUE)})
    ))
    #determine the number of recovered on this day
    rec <- rbind(rec, data.frame(dpch = sampleday.vars[day],
                                 ndpch = day,
                                 aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                      by = list(aggregate.data$Group, 
                                                                aggregate.data$Vaccinated,
                                                                aggregate.data$Strain), 
                                                      FUN = function(x){sum(x==2, na.rm = TRUE)})
    ))
    #determine the number of birds on this day
    n <- rbind(n, data.frame(dpch = sampleday.vars[day],
                             ndpch = day,
                             aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                  by = list(aggregate.data$Group, 
                                                            aggregate.data$Vaccinated,
                                                            aggregate.data$Strain), 
                                                  FUN = function(x) sum( !is.na(x) )) 
                             
    ))
    
    #determine the number of susceptibles
    if(exclude.seeder){
      susceptibles <- rbind(susceptibles, data.frame(dpch = sampleday.vars[day],
                                                     ndpch = day,
                                                     aggregate.data.frame(aggregate.data[tolower(aggregate.data$Challenge) == "contact",sampleday.vars[day]], 
                                                                          by = list(aggregate.data$Group[tolower(aggregate.data$Challenge) == "contact"], 
                                                                                    aggregate.data$Vaccinated[tolower(aggregate.data$Challenge) == "contact"],
                                                                                    aggregate.data$Strain[tolower(aggregate.data$Challenge) == "contact"]), 
                                                                          FUN = function(x) {sum(x==0, na.rm = TRUE)}))) 
      
    }else
    {
      susceptibles <- rbind(susceptibles, data.frame(dpch = sampleday.vars[day],
                                                     ndpch = day,
                                                     aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                                          by = list(aggregate.data$Group, 
                                                                                    aggregate.data$Vaccinated,
                                                                                    aggregate.data$Strain), 
                                                                          FUN = function(x) {sum(x==0, na.rm = TRUE)})))
    }
    
    #determine the number of new cases at this day
    if(day < length(sampleday.vars)){
      
      if(exclude.seeder){
        cases <- rbind(cases, data.frame(dpch = sampleday.vars[day],
                                         ndpch = day,
                                         cases =  aggregate.data.frame(
                                           ((aggregate.data[tolower(aggregate.data$Challenge) == "contact",sampleday.vars[day+1]]==1)&(aggregate.data[tolower(aggregate.data$Challenge) == "contact",sampleday.vars[day]]==0)), 
                                           by = list(aggregate.data$Group[tolower(aggregate.data$Challenge) == "contact"], 
                                                     aggregate.data$Vaccinated[tolower(aggregate.data$Challenge) == "contact"],
                                                     aggregate.data$Strain[tolower(aggregate.data$Challenge) == "contact"]), 
                                           FUN = sum, na.rm = TRUE)   )      )
      }
      else
      {
        cases <- rbind(cases, data.frame(dpch = sampleday.vars[day],
                                         ndpch = day,
                                         cases =  aggregate.data.frame((
                                         aggregate.data[,sampleday.vars[day+1]]-aggregate.data[,sampleday.vars[day]])>0, 
                                         by = list(aggregate.data$Group, aggregate.data$Vaccinated,aggregate.data$Strain), 
                                         FUN = sum, na.rm = TRUE)))
      }
    }else {
      
      #determine cases and potential exclude the seeders
      if(exclude.seeder){
        cases <- rbind(cases, 
                       data.frame(dpch = sampleday.vars[day],
                                  ndpch = day,
                                  cases = aggregate.data.frame(aggregate.data[tolower(aggregate.data$Challenge) == "contact",sampleday.vars[day]], 
                                                               by = list(aggregate.data$Group[tolower(aggregate.data$Challenge) == "contact"],
                                                                         aggregate.data$Vaccinated[tolower(aggregate.data$Challenge) == "contact"],
                                                                         aggregate.data$Strain[tolower(aggregate.data$Challenge) == "contact"]), 
                                                               FUN = function(x){0}))
        )}
      else{
        cases <- rbind(cases, 
                       data.frame(dpch = sampleday.vars[day],
                                  ndpch = day,
                                  cases = aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 

                                                               by = list(aggregate.data$Group, aggregate.data$Vaccinated,aggregate.data$Strain), 
                                                               FUN = function(x){0}))
        )}
    }
  }
  
  #combine infectious, population and cases
  output = cbind(inf,
                 S = susceptibles$x, 
                 C = cases$cases.x,
                 R = rec$x,
                 N = n$x)
  

  colnames(output)[c(3:10)]<- c("Group","Vaccinated","Strain","I","S","C","R","N")

  return(list(output,aggregate.data))
  
  
}


###Determine state based on model ####
SIR.state<- function(in.data,#vector of consecutive samples
                     model,#type of transmission model to use
                     inf.rule =1,
                     rec.rule =1,
                     only.last = T, #if only last cell is positive remove
                     ani_nr = -999
){
  out.data <- in.data;
  #determine potential infection and recovery events
  change <- c(0,in.data[2:length(in.data)]-in.data[1:(length(in.data)-1)])
  #replace NA by change
  change[is.na(change)]<- 1
  #replace a positive first sample by a change
  if(!is.na(in.data[1])){if(in.data[1]==1){change[1]<-1}}
  
  #then based on the model recode it to being susceptible (0), infectious (1) or recovered (2)
  
  if(model == "SI")
  {
    #if the first sample was positive return a vector of 1's
    if(!is.na(in.data[1])){
    if(in.data[1] == 1)return(rep(1,length(out.data)))}
    #first positive samples and loop such that infection rule is fullfilled
    found = F; index = 1;index.first <- length(change)+1
    while(index <= length(change) & found == F & max(change==1,na.rm = T) == 1 )
    {
      index.first <- c(1:length(out.data))[change==1][index]
      index.first.na <- c(1:length(out.data))[is.na(change)][index]
      if(!is.na(index.first.na))
      {
        if(out.data[index.first.na]==1){index.first <- min(index.first,index.first.na)}}
      if(!is.na(index.first)){
        if(min(in.data[index.first:min(index.first+inf.rule-1, length(in.data))])==1)
        {
          found = T
          }
        else
        {
          #false positive removal
          out.data[index.first]<- 0;
          #add 
          index = index + 1;}
      }else{index = index + 1;}
    }
    #get the NA's
    index.nas <- c(1:length(out.data))[is.na(out.data)]
    #set all to one
    if(!is.na(index.first)){
      if(index.first <= length(out.data)){
        out.data[index.first:length(out.data)] <- 1}
      }
    out.data[index.nas] <- NA
  }
  
  if(model == "SIR")
  {
    #first positive samples and loop such that infection rule is full filled
    found = F; index = 1;index.first <- length(change)+1
    while(index <= length(change) & found == F & max(change==1,na.rm = T) == 1)
    {
      index.first <- c(1:length(out.data))[change==1][index]
      
      if(!is.na(index.first)){
       
        if(min(in.data[index.first:min(index.first+inf.rule-1, length(in.data))])==1){
          found = T
          if(index.first > 1){
          #remove changes before the first index 
          change[1:(index.first-1)] <- 0
          #set all to be susceptible
          out.data[1:(index.first-1)]<-0
          }
        }
        else{
          #false positive removal
          out.data[index.first]<- 0;
          #add 
          index = index + 1;}
      }else{
        #or actual value is false
        index = index + 1;
        }
    }
    #if the first index is larger than the length of the data set it should actually be the first positive sample
    if(!is.na(index.first)){
    if(found == F & index.first>length(change) & min(unlist(change[!is.na(change)]))==-1){
      index.first <- min(c(1:length(change))[!is.na(out.data)])
      }
    }
  
    
    #first negative samples and loop such that infection rule is full filled
    found = F; index = 1;index.last <- length(change)+1; 
    while(index <= length(change) & found == F & sum(out.data)>0)#there are changes and found is false and there is at least one moment that the animal became positive.
    {
      index.last <- c(1:length(out.data))[change== -1 & !is.na(change)][index]
      if(!is.na(index.last))
      {
        if(max(in.data[index.last:min(index.last+rec.rule-1, length(in.data))])==0)
        {
          found = T
        }
        else{
          #false negative removal
          out.data[index.last]<- 1;
          #add 
          index = index + 1;}
      } else{index = index + 1}
    }
    #get the NA's from the original data
    index.nas <- c(1:length(in.data))[is.na(in.data)]
    if(is.na(index.first)|index.first > length(out.data))
    {
      index.first <- NA #no infection
      index.last<- NA #no recovery without infection
    }
    
    #set all to one
    if(!is.na(index.first)){
      out.data[index.first:length(out.data)] <- 1 #infectious
      if(!is.na(index.last)){
        out.data[(index.last):length(out.data)] <- 2}} #recovered
    out.data[index.nas] <- NA
  }
  
  
  #if(model == "SIS") {stop("not yet implemented")}
  
  if(only.last){
    if(sum(as.numeric(out.data[1,]),na.rm = T)==1){
      if(tail(out.data[!is.na(out.data)], n = 1) == 1)
      {
        out.data[!is.na(out.data)]<- 0}#set all to be 0
    }
    
  }
  
  #replace na's in the data
  
  
  return(out.data)
}

