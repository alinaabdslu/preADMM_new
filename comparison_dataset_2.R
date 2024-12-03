# Use the same algorithm as for function preADMM()
# library(splatter)
# library(scPADGRN)
# library(tsibble)


comparison_dataset_2 = function(repetitions){
  
  reps = repetitions
  
  inside_function = function(){
    
    net <- preADMM(N,cluGeneExp,0.1,0.1,0.01)[[1]]
    for(i in c(1:(N-1))){
      diag(net[[i]])<-0
      net[[i]][net[[i]]>0]<- net[[i]][net[[i]]>0]/max(net[[i]][net[[i]]>0])
      net[[i]][net[[i]]<0]<- net[[i]][net[[i]]<0]/max(abs(net[[i]][net[[i]]<0]))
      colnames(net[[i]])<-tf[,1]
      rownames(net[[i]])<-tf[,1]
      net[[i]][abs(net[[i]])<0.2268]<-0
      net[[i]] <- abs(net[[i]])
      net[[i]][net[[i]]>0] <- 1
      net[[i]][net[[i]]<0] <- -1
    }
    net_original = net
    
    net <- preADMM_new(N,cluGeneExp,0.1,0.1,0.01)[[1]]
    for(i in c(1:(N-1))){
      diag(net[[i]])<-0
      net[[i]][net[[i]]>0]<- net[[i]][net[[i]]>0]/max(net[[i]][net[[i]]>0])
      net[[i]][net[[i]]<0]<- net[[i]][net[[i]]<0]/max(abs(net[[i]][net[[i]]<0]))
      colnames(net[[i]])<-tf[,1]
      rownames(net[[i]])<-tf[,1]
      net[[i]][abs(net[[i]])<0.2268]<-0
      net[[i]] <- abs(net[[i]])
      net[[i]][net[[i]]>0] <- 1
      net[[i]][net[[i]]<0] <- -1
    }
    net_new = net
    
    M_1_original = as.matrix(net_original[[1]])
    M_2_original = as.matrix(net_original[[2]])
    M_3_original = as.matrix(net_original[[3]])
    
    M_1_new = as.matrix(net_new[[1]])
    M_2_new = as.matrix(net_new[[2]])
    M_3_new = as.matrix(net_new[[3]])
    
    M_1 = M_1_new - M_1_original
    M_2 = M_2_new - M_2_original
    M_3 = M_3_new - M_3_original
    
    result = c(sum(abs(c(M_1))),
               sum(abs(c(M_2))),
               sum(abs(c(M_3))))
    
    result}
  
  result = replicate(reps, {inside_function()})
  
  difference_1 = mean(result[1,])
  difference_2 = mean(result[2,])
  difference_3 = mean(result[3,])
  
  data_tsib <- tsibble(
    value = c(difference_1, difference_2, difference_3),
    Index = c(1, 2, 3),
    index = Index
  )
  
  data_tsib$value
  
}
