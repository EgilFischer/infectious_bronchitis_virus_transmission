##############################
##############################

library(ggplot2)
library(reshape2)
# library(rstan)
# library(cowplot)


############################
# plotting functions       #
############################
#plot the raw data when in an untidy format
raw.plot <-function(raw.data){
  raw.data.melt<- reshape2::melt(raw.data,
                                 id.vars = list("Group",
                                                "bird.id",
                                                "Vaccinated",
                                                "Challenge",
                                                "Strain"),
                                 variable.name = "day")
  raw.data.melt$ndpch <- reshape2::colsplit(raw.data.melt$day,"[.]",c("ndpch","X"))[,1]
  raw.data.melt$ndpch <-as.numeric(raw.data.melt$ndpch)
  raw.data.melt$delta <- 0;
  for(id in unique(raw.data.melt$bird.id)){
    raw.data.melt[raw.data.melt$bird.id== id,]$delta <-
      raw.data.melt[raw.data.melt$bird.id== id,]$ndpch-c(0,head(raw.data.melt[raw.data.melt$bird.id== id,]$ndpch,-1))  
  }
  
  
  return(raw.data.melt)
}

# head(raw.plot(sir.data[[2]]))
# 
# ggplot(raw.plot(sir.data[[2]]),
#        aes(as.factor(bird.id),fill = as.factor(value)))+
#   geom_bar(position = "stack")+coord_flip()+
#   facet_grid(Strain~Vaccinated )

################
# change data  #
################

transform.data <- function(sirdata){
  out <-data.frame(animal_id =c(),
                   inoculation=c(),
                   group=c(),
                   treatment=c(), 
                   lastneg = c(),
                   firstpos=c(),
                   lastpos=c(),
                   firstrec=c(),
                   lastObs=c())
  #iterate data
  for(i in c(1:length(sirdata[,1])))
  {
          out <- rbind(out,
                 data.frame(animal_id = sirdata[i,"bird.id"],
                            inoculation = sirdata[i,"Challenge"],
                            group = sirdata[i,"Group"],
                            treatment = sirdata[i,"Vaccinated"], 
                            lastneg = max(0,min(which(sirdata[i,] == 1)-5)-1),
                            firstpos = min(which(sirdata[i,] == 1)-5),
                            lastpos = max(which(sirdata[i,] == 1)-5),
                            firstrec = min(which(sirdata[i,] == 2)-5),
                            lastObs = max(which(!is.na(sirdata[i,]))-5))
    )
    
  }
  return(out)
}

#plot the input from a transmission experiment
#data require the following headings: group, animal_id, inoc,lastneg, firstpos, lastpos, firstrec, lastObs
input.plot <- function(data, color_scheme = c("#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), period.labels = c("Negative", 'Neg->Pos', 'Positive', 'Pos-> Neg',"Negative again")){
  plot.data <- data[,c(1:4)];#only select indicators
    #replace all negative numbers for lastObs
  data$lastneg[data$lastneg<0 ]<- data$lastObs[data$lastneg<0 ]
  data$firstpos[data$firstpos<0 ]<- data$lastObs[data$firstpos<0 ]
  data$lastpos[data$lastpos<0 ]<- data$lastObs[data$lastpos<0 ]
  data$firstrec[data$firstrec<0 ]<- data$lastObs[data$firstrec<0 ]
  
  #replace real time with interval
  plot.data$neg1 <- data$lastneg;
  plot.data$negpos <- data$firstpos - data$lastneg;
  plot.data$pos <-data$lastpos - data$firstpos;
  plot.data$posneg2 <- data$firstrec - data$lastpos;
  plot.data$neg2 <-  data$lastObs-data$firstrec;
  
  
  #melt down
  plot.data <- melt(plot.data,id = c("animal_id","inoculation","group","treatment"))
 #create a plot
  p1<-ggplot(data = plot.data[plot.data$treatment == "No",])+
    geom_bar(aes(y = value, x = as.factor(animal_id), fill = factor(variable,levels = c("neg1","negpos","pos","posneg2","neg2")[c(5,4,3,2,1)])),stat= 'identity')+
    coord_flip()+ 
    scale_fill_manual('State', 
                     values = color_scheme,
                     labels = period.labels[c(5,4,3,2,1)])+ylab("Days post challenge")+
    xlab("animal ID")+
    ggtitle("Unvaccinated")+
    scale_y_continuous(breaks = c(0:14) )+
    scale_x_discrete(limits = rev)+
          if(length(unique(plot.data$group))>1){facet_grid(group + inoculation ~ ., scales = "free_y")}
  
  p2<-ggplot(data = plot.data[plot.data$treatment == "Yes",])+
    geom_bar(aes(y = value, x = as.factor(animal_id), fill = factor(variable,levels = c("neg1","negpos","pos","posneg2","neg2")[c(5,4,3,2,1)])),stat= 'identity')+
    coord_flip()+ 
    scale_fill_manual('State', 
                      values = color_scheme,
                      labels = period.labels[c(5,4,3,2,1)])+ylab("Days post challenge")+
    xlab("animal ID")+
    ggtitle("Vaccinated")+
    scale_y_continuous(breaks = c(0:14) )+
    scale_x_discrete(limits = rev)+
    if(length(unique(plot.data$group))>1){facet_grid(group + inoculation ~ ., scales = "free_y")}
  return(list(plot1 = p1,plot2 = p2, data = plot.data))
}


#plot the input from a transmission experiment and the distribution of estimated unobserved data 
#data require the following headings: group, animal_id, inoc,lastneg, firstpos, lastpos, firstrec, lastObs
estimated.times.plot <- function(stan.fit, input.data){
  #extract estimated transitions
  x<-data.frame(extract(stan.fit)$infT)
  names(x) <- c(input.data$animal_id)
  ex.data.infT <- reshape2::melt(x)
  ex.data.infT$group <- rep(input.data$group, each = length(x[,1]))
  ex.data.infT$treatment <- rep(input.data$treatment, each = length(x[,1]))
  
  x<-data.frame(extract(stan.fit)$latT)
  names(x) <- c(input.data$animal_id)
  ex.data.latT <- reshape2::melt(x)
  
  ex.data.latT$group <- rep(input.data$group, each = length(x[,1]))
  ex.data.latT$treatment <- rep(input.data$treatment, each = length(x[,1]))
  
 
  
  #add extracted estimation times to the input plot
  p <- input.plot(input.data)$plot+
     geom_violin(data = ex.data.infT, aes(x =variable, y = value), colour = "gray", fill = "yellow")+
     geom_violin(data = ex.data.latT, aes(x =variable, y = value), colour = "gray", fill ="green")+
      coord_flip()+
     scale_x_discrete()+
     scale_y_discrete(limits = rev)+
     xlab("Animals")+
     ylab("time")+
     if(length(unique(input.data$group))>1){facet_grid(group ~ treatment, scales = "free_y")}
   return(list(plot = p))
}


#plot prior versus the posterior distributions 
#requires the following input: prior-distributions, prior-parameters, parameter names, stan.fit object
prior_posterior.plot <- function(prior.dist, prior.par, par.name, stan.fit){
  #extract estimated transitions
  d<-data.frame(k = extract(stan.fit)[par.name])
  dist <- get(prior.dist[par.name])
  pr.dist <- function(x){dist(x,  unlist(prior.par[par.name])[1],  unlist(prior.par[par.name])[2])}
  pp.plot <- ggplot()  +
    geom_density(aes(x = d[,par.name], y = ..density..), fill = "grey", alpha = .7, size = 0.5)   +
    geom_area(aes(x = x, y = y), 
              data = data.frame(x = seq(0,1.5 * max(d[,par.name]),max(d[,par.name])/30), 
                                y = sapply(X = seq(0,1.5 * max(d[,par.name]),max(d[,par.name])/30), FUN = pr.dist)),
              size =1,
              colour = "red",
              fill ="red",
              alpha = 0.1)  + xlab(par.name) 
  
  #return plot
  return(pp.plot)
}

#plot all
prior_posterior.plots <- function(prior.dist, prior.par,  stan.fit, print.plot = T){
  plots <- list();
  pars <- names(prior.dist)
  for(p.n in pars)
  {
    plots[[p.n]] <- prior_posterior.plot(prior.dist, prior.par, p.n, stan.fit);
  }
  if(print.plot){pp.grid <- plot_grid(plotlist = plots);
    return(pp.grid)}
  return(plots)
  
}
