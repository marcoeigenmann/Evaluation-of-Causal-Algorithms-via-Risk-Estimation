source("~/UtilityFunctions.R")
library(mvtnorm)
library(pcalg)
library(doParallel)
library(stringr)


p = NULL
IntMean=NULL
OracleNumberDe=F
no_cores = 4
BiDir <- T
i <- 1

prov <- NULL


### Select a subset of the simulations (if desired) (here shown for p and IntMean)

### p ###
filelist <- NULL
filelist <- list.files(path = "~/", pattern = "p=")#"/scratch/userdata/marcoei/May/Simtdistr"  ## "/scratch/userdata/marcoei/Apr2019"
filelistDiscard <- list.files(path = "~/", pattern = "DataSim")
filelist <- setdiff(filelist, filelistDiscard)





### IntMean ###
if(!is.null(IntMean)){
  prov <- paste(prov, "IntMean=",IntMean,"type",sep="")
  prov.filelist <- list.files(path = "~/", pattern = prov)
  prov <- NULL
  
  if(!is.null(filelist)){
    filelist <- intersect(filelist, prov.filelist)
  }else{
    filelist <- prov.filelist
  }
}




### Main function
if(length(filelist)>0){
  cat("Considering ", length(filelist), " datasets..")
  
  load(paste("~/",filelist[1],sep=""))
  
  q.value <- 4
  
  
  ### Parallelize
  registerDoParallel(cores = no_cores)
  NotNeeded <- foreach(i_par=1:length(filelist)) %dopar% {


    
    
    
    L <- Lstar <- LCV5 <- LCV5star <- LIota <- LIotastar <-  matrix(data = 0, nrow = 1, ncol = 7)

    
    load(paste("~/",filelist[i_par],sep=""))
    
    
    
    ### We stop and set to 1 Low and LowInt if we do not have enough descendants/interventions
    
    if(sum(tosave$DAG[[1]][unique(tosave$TStatistics[[1]][,2]),]!=0)==0){
      return()
    }
    if(sum(tosave$DAG[[1]][unique(tosave$TStatistics[[1]][,2]),]!=0)<=3){
      Low <- 1
    }else{
      Low <- 0
    }
    
    if(length(unique(tosave$TStatistics[[1]][,2]))<2){
      LowInt <- 1
    }else{
      LowInt <- 0
    }      
    
    p <- nrow(tosave$DAG[[1]])
    

    ### Get the interventions and the statistics
    interventions <- unique(tosave$TStatistics[[1]][,2])
    n.int <- length(interventions)
    Tstatistics <- tosave$TStatistics[[1]]
    order.Tstatistics <- order(abs(Tstatistics[,1]), decreasing = TRUE)
    ordered.Tstatistics <- Tstatistics[order.Tstatistics,]
    ordered.Tstatistics.extended <- ordered.Tstatistics
    if(length(interventions)!=p){
      for(i.allnodes in 1:length(setdiff(1:p,interventions))){
        toAdd <- matrix(data = 0, nrow = p-1, ncol = 3)
        toAdd[,2] <- setdiff(1:p,interventions)[i.allnodes]
        toAdd[,3] <- setdiff(1:p,toAdd[1,2])
        ordered.Tstatistics.extended <- rbind(ordered.Tstatistics.extended, toAdd)
      }
    }
    
    
    samplesize <- as.numeric(substr(str_match(string = filelist[[1]], pattern = "n=.*ENS"), 3,5))

    
    GIES <- tosave$GIES.tosave[[1]]
    GES <- tosave$GES.tosave[[1]]
    DAG <- tosave$DAG[[1]]
    DAG[DAG!=0] <- 1
    TransClosedDAG <- transClos(DAG)
    ### TC stands for transitive closure
    NrEdgesTCDAG_int <- sum(TransClosedDAG[interventions,])
    
    
    ### Compute the descendants either with oracle infromation or with a cut-off
    De <- matrix(data = 0, nrow = p,ncol = p)
    large.Tstatistics <- sum(abs(ordered.Tstatistics[,1])>q.value)
    if(length(interventions)>0){
      if(OracleNumberDe){
        q <- 1
        while(NrEdgesTCDAG_int > sum(transClos(De))){
          De[ordered.Tstatistics[q,2], ordered.Tstatistics[q,3]] <- 1
          q <- q+1
        }
      }else{
        col.q.value <- which.min((abs(n.int/p - c(0.1,0.2,0.3,0.4,0.5,1))))
        row.q.value <- which.min((abs(p - c(10,25,50,100,200,400))))
        ### Use the value you consider more appropriate
        q.value <- 4
        large.Tstatistics <- sum(abs(ordered.Tstatistics[,1])>q.value)
        for(q in 1:large.Tstatistics){
          De[ordered.Tstatistics[q,2], ordered.Tstatistics[q,3]] <- 1
        }
      }
    }

    
    ### Compute the correaltion based graph estimates
    Cov <- tosave$CovMatAll[tosave$invPermutation, tosave$invPermutation]
    AbsCov <- abs(Cov)
    Cor <- cov2cor(Cov)
    AbsCor <- abs(Cor)
    Cor02 <- Cor
    Cor02[AbsCor>=0.2] <- 1
    Cor02[AbsCor<0.2] <- 0
    diag(Cor02) <- 0
    Cor05 <- Cor02
    Cor05[AbsCor<0.5] <- 0
    diag(Cor05) <- 0
    Cor08 <- Cor05
    Cor08[AbsCor<1] <- 0
    diag(Cor08) <- 0
    MaxCorEntry <- which(max(((AbsCor - diag(x = 1, nrow = p, ncol = p))[interventions,]))==(AbsCor - diag(x = 1, nrow = p, ncol = p)), arr.ind=TRUE)
    MaxCorEntry.indices <- matrix(data = interventions[MaxCorEntry], nrow = 2, ncol = 2)
    Cor02 <- Cor02*0
    Cor08 <- Cor
    Cor08[AbsCor < sort(AbsCor-diag(x = 1, nrow = p, ncol = p), decreasing = T)[2*sum(tosave$DAG[[1]]!=0)]] <- 0
    Cor08[Cor08 != 0] <- 1
    diag(Cor08) <- 0
    Cor02[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
    Cor05[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
    Cor08[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
    Cor02[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
    Cor05[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
    Cor08[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
    Cor02 <- transClos(Cor02)
    Cor05 <- transClos(Cor05)
    Cor08 <- transClos(Cor08)
    
    
    
    ### Compute the extended naive risk estimator
    GIES <- GIES[tosave$invPermutation, tosave$invPermutation]
    GIESBiDir <- GIES
    colnames(GIESBiDir) <- colnames(Cor05)
    rownames(GIESBiDir) <- rownames(Cor05)

    if(!BiDir){
      GIES[GIES==t(GIES)] <- 0
    }
    

    
    GES <- GES[tosave$invPermutation, tosave$invPermutation]
    GESBiDir <- GES
    colnames(GESBiDir) <- colnames(Cor05)
    rownames(GESBiDir) <- rownames(Cor05)
    if(!BiDir){
      GES[GES==t(GES)] <- 0
    }
    

    

    DAGWithoutNames <- tosave$DAG[[1]]
    dimnames(DAGWithoutNames) <- NULL
    IDAGStar <- DAGWithoutNames*0
    
    


    if(length(interventions)>0){
      for(n.int in 1:p){
        
        L[i,1] <- L[i,1] + RiskEstimation(G1 = De, G2 = transClos(GES), i = n.int, type = "Jaccard")[[2]]/p
        L[i,3] <- L[i,3] + RiskEstimation(G1 = De, G2 = transClos(GIES), i = n.int, type = "Jaccard")[[2]]/p
        L[i,5] <- L[i,5] + RiskEstimation(G1 = De, G2 = Cor02, i = n.int, type = "Jaccard")[[2]]/p
        L[i,6] <- L[i,6] + RiskEstimation(G1 = De, G2 = Cor05, i = n.int, type = "Jaccard")[[2]]/p
        L[i,7] <- L[i,7] + RiskEstimation(G1 = De, G2 = Cor08, i = n.int, type = "Jaccard")[[2]]/p
        
        
        Lstar[i,1] <- Lstar[i,1] + RiskEstimation(G1 = TransClosedDAG, G2 = transClos(GES), i = n.int, type = "Jaccard")[[2]]/p
        Lstar[i,3] <- Lstar[i,3] + RiskEstimation(G1 = TransClosedDAG, G2 = transClos(GIES), i = n.int, type = "Jaccard")[[2]]/p
        Lstar[i,5] <- Lstar[i,5] + RiskEstimation(G1 = TransClosedDAG, G2 = Cor02, i = n.int, type = "Jaccard")[[2]]/p
        Lstar[i,6] <- Lstar[i,6] + RiskEstimation(G1 = TransClosedDAG, G2 = Cor05, i = n.int, type = "Jaccard")[[2]]/p
        Lstar[i,7] <- Lstar[i,7] + RiskEstimation(G1 = TransClosedDAG, G2 = Cor08, i = n.int, type = "Jaccard")[[2]]/p  
        
        
      }
      
      
      
      
    }
    ############ End
    
    
    ### Compute the naive risk estimator
    
    fPermutation <- invPerm(tosave$invPermutation)
    
    
    Cov <- tosave$CovMatAll[tosave$invPermutation, tosave$invPermutation]
    AbsCov <- abs(Cov)
    Cor <- cov2cor(Cov)
    AbsCor <- abs(Cor)
    Cor02 <- Cor
    Cor02[AbsCor>=0.2] <- 1
    Cor02[AbsCor<0.2] <- 0
    diag(Cor02) <- 0
    Cor05 <- Cor02
    Cor05[AbsCor<0.5] <- 0
    diag(Cor05) <- 0
    Cor08 <- Cor05
    Cor08[AbsCor<1] <- 0
    diag(Cor08) <- 0
    MaxCorEntry <- which(max(((AbsCor - diag(x = 1, nrow = p, ncol = p))[interventions,]))==(AbsCor - diag(x = 1, nrow = p, ncol = p)), arr.ind=TRUE)
    MaxCorEntry.indices <- matrix(data = interventions[MaxCorEntry], nrow = 2, ncol = 2)
    Cor02 <- Cor02*0
    Cor08 <- Cor
    Cor08[AbsCor < sort(AbsCor-diag(x = 1, nrow = p, ncol = p), decreasing = T)[2*sum(tosave$DAG[[1]]!=0)]] <- 0
    Cor08[Cor08 != 0] <- 1
    diag(Cor08) <- 0
    Cor02[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
    Cor05[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
    Cor08[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
    Cor02[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
    Cor05[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
    Cor08[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
    Cor02 <- transClos(Cor02)
    Cor05 <- transClos(Cor05)
    Cor08 <- transClos(Cor08)
    
    
    if(length(interventions)>0){
      
      p <- length(interventions)
      
      for(n.int2 in 1:length(interventions)){
        
        n.int <- interventions[n.int2]
        
        LIota[i,1] <- LIota[i,1] + RiskEstimation(G1 = De, G2 = transClos(GES), i = n.int, type = "Jaccard")[[2]]/p
        LIota[i,3] <- LIota[i,3] + RiskEstimation(G1 = De, G2 = transClos(GIES), i = n.int, type = "Jaccard")[[2]]/p
        LIota[i,5] <- LIota[i,5] + RiskEstimation(G1 = De, G2 = Cor02, i = n.int, type = "Jaccard")[[2]]/p
        LIota[i,6] <- LIota[i,6] + RiskEstimation(G1 = De, G2 = Cor05, i = n.int, type = "Jaccard")[[2]]/p
        LIota[i,7] <- LIota[i,7] + RiskEstimation(G1 = De, G2 = Cor08, i = n.int, type = "Jaccard")[[2]]/p
        
        
        LIotastar[i,1] <- LIotastar[i,1] + RiskEstimation(G1 = TransClosedDAG, G2 = transClos(GES), i = n.int, type = "Jaccard")[[2]]/p
        LIotastar[i,3] <- LIotastar[i,3] + RiskEstimation(G1 = TransClosedDAG, G2 = transClos(GIES), i = n.int, type = "Jaccard")[[2]]/p
        LIotastar[i,5] <- LIotastar[i,5] + RiskEstimation(G1 = TransClosedDAG, G2 = Cor02, i = n.int, type = "Jaccard")[[2]]/p
        LIotastar[i,6] <- LIotastar[i,6] + RiskEstimation(G1 = TransClosedDAG, G2 = Cor05, i = n.int, type = "Jaccard")[[2]]/p
        LIotastar[i,7] <- LIotastar[i,7] + RiskEstimation(G1 = TransClosedDAG, G2 = Cor08, i = n.int, type = "Jaccard")[[2]]/p  
        
        
        
      }
      p <- ncol(tosave$DAG[[1]])
    }
    ############ End
    
    
    
    ### Compute the CV based risk estimator
    ### We use 5-fold CV
    if(length(interventions)>0){
    
      
      
      

      fPermutation <- invPerm(tosave$invPermutation)
      
      ### the 5 estimates GIES
      GIESAll <- tosave$CV5ToSave[[1]]
      
      ### the 5 estimates GES
      GESAll <- tosave$CV5ToSave[[2]]
      
      ### the 5 set of interventions
      intervs <- tosave$CV5ToSave[[3]]
      
      ### the 5 covariance matrices
      CovAll <- tosave$CV5ToSave[[4]]
      
      
      
      length.u <- 0
      for(u in 1:length(GIESAll)){
        length.u <- length.u + length(intervs[[u]])
      }
      
      if(!is.null(GIESAll)){
        if(length(interventions)>0){
          for(u in 1:length(GIESAll)){
            
            Cov <- CovAll[[u]][tosave$invPermutation, tosave$invPermutation]
            AbsCov <- abs(Cov)
            Cor <- cov2cor(Cov)
            AbsCor <- abs(Cor)
            Cor02 <- Cor
            Cor02[AbsCor>=0.2] <- 1
            Cor02[AbsCor<0.2] <- 0
            diag(Cor02) <- 0
            Cor05 <- Cor02
            Cor05[AbsCor<0.5] <- 0
            diag(Cor05) <- 0
            Cor08 <- Cor05
            Cor08[AbsCor<1] <- 0
            diag(Cor08) <- 0
            MaxCorEntry <- which(max(((AbsCor - diag(x = 1, nrow = p, ncol = p))[interventions,]))==(AbsCor - diag(x = 1, nrow = p, ncol = p)), arr.ind=TRUE)
            MaxCorEntry.indices <- matrix(data = interventions[MaxCorEntry], nrow = 2, ncol = 2)
            Cor02 <- Cor02*0
            Cor08 <- Cor
            Cor08[AbsCor < sort(AbsCor-diag(x = 1, nrow = p, ncol = p), decreasing = T)[2*sum(tosave$DAG[[1]]!=0)]] <- 0
            Cor08[Cor08 != 0] <- 1
            diag(Cor08) <- 0
            Cor02[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
            Cor05[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
            Cor08[MaxCorEntry[1,1],MaxCorEntry[1,2]] <- 1
            Cor02[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
            Cor05[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
            Cor08[MaxCorEntry[1,2],MaxCorEntry[1,1]] <- 1
            Cor02 <- transClos(Cor02)
            Cor05 <- transClos(Cor05)
            Cor08 <- transClos(Cor08)
            
            GIES <- GIESAll[[u]]
            GES <- GESAll[[u]]

            GIES <- GIES[tosave$invPermutation, tosave$invPermutation]
            GES <- GES[tosave$invPermutation, tosave$invPermutation]

            if(!BiDir){
              GIES[GIES==t(GIES)] <- 0
              GES[GES==t(GES)] <- 0
            }

            
            for(n.int2 in 1:length(intervs[[u]])){
              
              n.int <- fPermutation[intervs[[u]]][n.int2]
              
              LCV5[i,1] <- LCV5[i,1] + RiskEstimation(G1 = De, G2 = transClos(GES), i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5[i,3] <- LCV5[i,3] + RiskEstimation(G1 = De, G2 = transClos(GIES), i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5[i,5] <- LCV5[i,5] + RiskEstimation(G1 = De, G2 = Cor02, i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5[i,6] <- LCV5[i,6] + RiskEstimation(G1 = De, G2 = Cor05, i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5[i,7] <- LCV5[i,7] + RiskEstimation(G1 = De, G2 = Cor08, i = n.int, type = "Jaccard")[[2]]/length.u
              
              
              LCV5star[i,1] <- LCV5star[i,1] + RiskEstimation(G1 = DAG, G2 = transClos(GES), i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5star[i,3] <- LCV5star[i,3] + RiskEstimation(G1 = DAG, G2 = transClos(GIES), i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5star[i,5] <- LCV5star[i,5] + RiskEstimation(G1 = DAG, G2 = Cor02, i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5star[i,6] <- LCV5star[i,6] + RiskEstimation(G1 = DAG, G2 = Cor05, i = n.int, type = "Jaccard")[[2]]/length.u
              LCV5star[i,7] <- LCV5star[i,7] + RiskEstimation(G1 = DAG, G2 = Cor08, i = n.int, type = "Jaccard")[[2]]/length.u
              
            }
          }
        }
        
        ############ End of CV5
        
      }
    }
    

    
    setwd("~/")
    DataToSave <- list(L, Lstar, LIota, LIotastar, LCV5, LCV5star, Low, LowInt)
    save(DataToSave,file = paste("Data",filelist[[i_par]],"OracleNumberDe=",OracleNumberDe,"func=",tosave$fn,"DoInt=",tosave$DoInt,"Low=",Low,"LowInt=",LowInt,sep=""), version = 2)
  }
}


