######################################################################################
#                                                                                    #
#                   Final size estimation with a fast implementation                 #
#                                                                                    #
#                     Copyright R-code: Egil A.J.Fischer                             #
#                                       e.a.j.fischer@uu.nl/egil@egilfischer.nl      #
###################################################################################### 

require(tidyverse)
### this function will first create a generation based table with states and then gives the final size distribution as 
#output
distFSfast <- function(R,s0in,i0in ,r0in = NULL)
  {
    #set size of the population 
    nin = s0in+i0in;
    if(!is.null(r0in)) 
           nin = s0in+i0in+r0in;
    
    #infection probability
    pinf =function(s,n){s*R/(s*R + n)};
    
     #create final size distribution for each s0 and i0 combination
     max.cases <- max(s0in);
     gen.mat.out = matrix(0,ncol = (max.cases+1),nrow = length(s0in));
     for(m in 1:length(s0in))
       {
       s0= s0in[m];
       i0= i0in[m];
       n = nin[m];
       #set the state matrix
       gen.mat <- matrix(rep(0, (s0+i0+1)*(s0+1)),nrow= (s0+1));
       gen.mat[s0+1,i0+1]<- 1;
       #fill first row 
       for(j in c((i0+1):1)){
         gen.mat[s0+1,j-1]<- gen.mat[s0+1,j]*(1-pinf(s0,n)) 
         }
          
        #next rows added by infection
        for(i in c(1:(s0+1))){
          for(j in c((s0+i0+1):1)){
            #enter state by infection (only for j > 2, j = 1 no infections, j = 2 is 1 infection)
            if(j>2){
              gen.mat[s0+1-i,j]<-gen.mat[s0+1-i,j] + gen.mat[s0+2-i,j-1]*pinf(s0+1-i,n)
             };
          if(j < s0+i0+1){
            #enter state by recovery
            gen.mat[s0+1-i,j]<- gen.mat[s0+1-i,j] + gen.mat[s0+1-i,j+1]*(1-pinf(s0-i,n));
            }
           
          }
        }
    #store the final size distribution -> reversed so that it is 0 cases first
    gen.mat.out[m,1:(s0+1)] = rev(gen.mat[,1]);
     }
     #report it reversed 
     return(gen.mat.out)
}


pFS <- function(R,x,s0,i0,r0= NULL){
  #produce final size distribution for value r
  final.size.dist <- distFSfast(R,s0,i0,r0);
  #return the outcome for x 
  return(prod(mapply(function(i,j)final.size.dist[i,j+1],c(1:length(s0)),x)))
}


###Function to determine the probability of more extreme values given R####
# r = R, x = final number of cases, s0 = initial susceptibles, i0 = initial infectious, comp = the type of extreme 
#rm(pExtremes)

pExtremes<-  function(r,x,s0in,i0in,r0in = NULL, comp = `<`){
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0in),nrow = prod(s0in+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0in)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0in[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0in[k]),
                            FUN = function(x){rep(x,repetition)}), 
                     ncol = 1, 
                     nrow =prod(s0in+1))[,1]
  }
  #produce final size distribution for value r
  final.size.dist <- distFSfast(r,s0in,i0in, r0in);

  #define function for this distribution for the probability of a certain number of cases x in each of the trials
  pFSloc<- function(v){
    return(prod(mapply(function(i,j) {final.size.dist[i,j+1]},
                       c(1:length(s0in)),
                       v)))
  }
  
  #select the extremes by selecting those for which the total number of cases
  #and the sum of x are given by the comparison "comp"  thus either <,>,<= or >= 
  #calculate for each extreme the probability of the Final Size under the hypothesis R = r
  #and sum all probabilities of these extreme outcomes
  #if only one experiment
  if(length(s0in)==1)
  {
    return(sum(final.size.dist[out[comp(apply(out,1,sum),sum(x)),]+1]))
  }
  #If only one extreme outcome possible: 
   if(is.null(dim(out[comp(apply(out,1,sum),sum(x)),]))){
     #only requires the final size distribution
     return(pFSloc(out[comp(apply(out,1,sum),sum(x)),]))
     }
    #else
       return(c(sum(apply(out[comp(apply(out,1,sum),sum(x)),],
                          1,
                          function(ext){pFSloc(ext)}))))
}

##################################################################################################
#                                                                                                #
#                  Function to estimate R using the final size method                            #
#                  Includes confidence interval size 1- alpha and R >= 1 test                    #
#                                                                                                #
##################################################################################################

# Final Size function for x cases, s0 initially susceptible and i0 initially infectious animals.r0 for recovered
# Optional set sigficance level alpha, onesided testing 
# max.val is the maximum value used for optimization.  ####
FinalSize<- function(x,s0,i0,r0 = NULL, 
                     alpha = 0.05, 
                     onesided = FALSE, 
                     max.val = 250, 
                     decimals = 4){
  res <- data.frame(point.est = -999)
  #determine the point estimate by optimization of the log-likelihood function
  tmp <- if(sum(x)==0) {(list(minimum = c(0), objective = c(0)))}else
    {if(sum(x)==sum(s0)){(list(minimum = c(Inf),objective =  c(0)))}else
      {(optimize(interval = c(0.,max.val),
                       f = function(R){-log(pFS(R,x,s0,i0,r0))}))}}
                              
  res$point.est = tmp$minimum
  res$LL = tmp$objective
  #determine the confidence intervals
  #if one-sided is FALSE both sides, either only lower or upper limit of CI
  #lowerlimit is found for values of R for which the probability of extremes below the observations 
  res$ci.ll <- ifelse(sum(x)==0,0,
                      uniroot(interval = c(10^-10,max.val),extendInt = "yes",
                          f = function(R){(pExtremes(R,x,s0,i0,r0,comp = `>=`) - alpha / (2 - onesided))})$root)
  #upperlimit is found for values of R for which the probability of extremes above the observations
  res$ci.ul <- ifelse(sum(x)==sum(s0),Inf,
                      uniroot(interval = c(0,max.val),extendInt = "yes",
                              f = function(R){( pExtremes(R,x,s0,i0,r0,comp = `<=`) - alpha / (2 - onesided))})$root)
  
  #probability of R >= 1 is found be calculating the probability to find an equal or less positive under the assumption R0 = 1
  res$pval.above1 = pExtremes(1,x,s0,i0,r0,comp = `<=`)
  res$pval.below1 = pExtremes(1,x,s0,i0,r0,comp = `>=`)
  res$pval.equal1 = pExtremes(1,x,s0,i0,r0,comp = `>=`) * pExtremes(1,x,s0,i0,r0,comp = `<=`)
  return(signif(res, digits = decimals))
}



##################################################################################################
#                                                                                                #
#                  Function to compare R of two groups using the final size method               #
#                  ref: Velthuis et al2007                                                                                                #
##################################################################################################



Test.TwoPops<-  function(x,s0,i0,treat){
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0),nrow = prod(s0+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0[k]),FUN = function(x){rep(x,repetition)}), ncol = 1, nrow =prod(s0+1))[,1]
    
  }
  #select those for which the rowSum == sum(x)
  out<- out[rowSums(out)==sum(x),]
  
  #if needed change treatment to boolean
  treat = as.logical(str_replace(str_replace(c("Yes","No"),"Yes","TRUE"),"No","False"));
  
  #Select those possible outcomes for which the difference is >= sum(x[!treat]) - sum(x[treat])
  extremes <- out[sapply(X = apply(out[,!treat],1,sum)-apply(out[,treat],1,sum), FUN = function(z){abs(z) >=  sum(x[!treat])-sum(x[treat])}),]
  #and sum all probabilities of these same or extreme outcomes
  likeFun <- function(R){
    sum(apply(extremes,1,function(ext){pFS(R,ext,s0,i0)}))
  }
  #check extreme R-values if the difference is very large
  ext1 <- likeFun(0);
  ext2 <- likeFun(10);
  if(ext1 == ext2) return(ext1)
  #find R with the maximum probability of occurence, this is the p-value for R1 = R2
  return(optimize(interval = c(0,25),f = likeFun, maximum = TRUE)$objective)
}





