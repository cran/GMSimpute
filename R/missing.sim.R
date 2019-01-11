#' Missing peaks generating procedure for simulation study
#'
#' missing.sim generates various types of missing peaks based on specified missing proportion.  
#' @param complete.data The full abundance matrix without missing value, with features in rows and samples in columns.
#' @param total.missing A scalar or vector of proportions. It is the total percentage of missing peaks throughout the full matrix.  
#' @param random A scalar or vector of proportions. It is the percentage of random missing in all the missing peaks. 
#' @param pct.full A scalar for the percentage of alighned features (metabolites or peptides) without missing peaks.   
#' @param seednum The seed set for generating missing peaks index. Default seed is seednum=365.
#' @return 
#'    \item{simulated.data}{The list of all simulated scenarios}
#'    \item{Labels}{The description for each simulated scenario}
#' @examples 
#' data('tcga.bc.full')
#' # tcga.bc.full contains mass specturm abundance of 100 metabolites for 30 breast cancer 
#' # tumor and normal tissue samples without missing values.
#' 
#' 
#' simulated.data=missing.sim(tcga.bc.full,total.missing=c(0.2,0.4),random=c(0.3,0.5,0.7),pct.full=0.4)
#' # Generate missing (NA) values in full abundance matrix tcga.bc.full permuting all scenarios
#' 
#' @export 

missing.sim<-function(complete.data,total.missing,random,pct.full,seednum=365){
  if(sum(total.missing>1)>0)stop('Total missing percentage must be less than 1')
  if(sum(total.missing<0)>0)stop('Total missing percentage must be greater than 0')
  if(sum(random>1)>0)stop('Random missing percentage must be less than 1')
  if(sum(random<0)>0)stop('Random missing percentage must be greater than 0')
  
  p_miss=total.missing
  
  p_randOfMiss=random[random<1]
  
  n_probs1 <- length(p_miss)
  n_probs2 <- length(p_randOfMiss)
  if(1 %in% random){
  p_Allrand=total.missing
  n_probs3 <- length(p_Allrand)
  } else {n_probs3=0}
  
  tot_data   <- ncol(complete.data)#total number of samples
  
  set.seed(seednum)
  ct=1 # to loop through where datasets will be saved in list
  list_names <- rep(NA,(n_probs1*n_probs2+n_probs3))#we will now have a scenario that samples from low end, not just taking the entire low portion
  list_description <- rep(NA,n_probs1*n_probs2+n_probs3)
  data_miss   <- vector("list", n_probs1*n_probs2+n_probs3)
  
  
  for (i in 1:n_probs1){# go through each overall missing probability
    p_missi <- p_miss[i]
    for (k in 1:n_probs2) {#go through each probability of being random if missing (the remaining are missing at low end)
      p_randOfMissk <- p_randOfMiss[k]
      p_rand <- p_missi*p_randOfMissk  #probability of missing at random (overall)
      p_low  <- p_missi*(1-p_randOfMissk) #probability of missing at LOD (overall)
      
      n_rand_data   <- round(p_rand*tot_data)# total number of samples that are missing at random
      
      
      # find indices to treat as random missing
      if(n_rand_data!=0){
        miss_rand_data   <- t(apply(complete.data,1,function(x){x[sample(1:length(x),n_rand_data)]<-NA;return(x)}))#indices to delete (random)
         
      }else{
        miss_rand_data=complete.data
        
      }
      
      
      
      # get missing portion at low end (remove bottom portion)
      vec_data   <- as.vector(miss_rand_data)
      
      n_low_data    <- round(p_low*length(vec_data)) # total number of sample missing at low end
      
      indP <- order(vec_data)
      
      toDel_low_data   <- indP[1:n_low_data] #indices to delete (low end)- Complete truncation
      
      #replace those indices with NA
      vec_datamiss <- vec_data
      vec_datamiss[toDel_low_data] <- NA
      
      
      data_miss[[ct]] <- matrix(vec_datamiss,ncol=dim(complete.data)[2])
      
      rownames(data_miss[[ct]])=rownames(miss_rand_data)
      
      list_names[ct] <- paste('pRand',gsub("[[:punct:]]", "", p_rand),'_pLow',gsub("[[:punct:]]", "", p_low),sep='')
      list_description[ct] <- paste('P(missing at random)=',p_rand,', P(missing at LOD)=', p_low,sep='')      
      ct=ct+1
      
      
    }
  }
  
  # the last scenario- all missing is at random
  if(1 %in% random){
  n_rand_data=p_Allrand*tot_data
  
  for(i in 1:n_probs3){
    
    miss_rand_data   <- t(apply(complete.data,1,function(x){x[sample(1:length(x),n_rand_data[i])]<-NA;return(x)}))#indices to delete (random)
    
    data_miss[[ct]]   <- miss_rand_data
    
    list_names[ct]  <- paste('pRand',gsub("[[:punct:]]", "", p_Allrand[i]),'_pLow0',sep='')
    list_description[ct] <- paste('P(missing at random)=',p_Allrand[i],', P(missing at LOD)=', 0,sep='')      
    ct=ct+1
  }}
  
  names(data_miss) <- list_names
   
  miss=list(data=data_miss, Labels=list_description)
  
  mask=list()
  nonmiss=round(nrow(complete.data)*pct.full)
  for(i in 1:length(miss$data)){
    miss.n=rowSums(is.na(miss$data[[i]]))
    total.nonmiss=max(sum(miss.n<3),nonmiss)
    nonmiss.id=order(miss.n)[1:total.nonmiss]
    mask[[i]]=rbind(complete.data[nonmiss.id,],miss$data[[i]][-nonmiss.id,])
    mask[[i]]=mask[[i]][rownames(complete.data),]
    filter=rowSums(!is.na(mask[[i]]))>0.3*ncol(mask[[i]])
    mask[[i]]=mask[[i]][filter,]
  }
  
  names(mask)=list_names
  output=list(simulated.data=mask,Labels=list_description)
  
  return(output)
}