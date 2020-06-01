### Compute the graph estimators to be used to evaluate the risk estimation framework
  
library(pcalg)
library(ggm)
source("~/UtilityFunctions.R")


CompGraphs <- function(p, n, ENS, nrep, interventions, distr = "Normal", method="er", par1=NULL,
                        par2=NULL, SNR.input=NULL, IntMean, TS_input = 10, lambda = 0.5, Simindex = 1,
                        max.kk = 10, normalize = FALSE, seed = NULL, df = 4, DoInt, fn=NULL){
  
  
  ## Initialize
  
  if(!is.null(seed)){
    set.seed(seed)
  }

  p.edge <- ENS/(p-1)
  
    
  if(!is.null(SNR.input)){
    if(SNR.input<1){
      SNR <-  (1-SNR.input)/SNR.input
    }else{
      SNR <- SNR.input
    }
  }else{
    SNR <- NULL
  }
    


  ## Generate all true DAGs
  DAG <- GenMatDAGs(nrep = nrep, p = p, prob = p.edge, method = method, par1 = par1, par2 = par2, fn = fn)
  FunctionType <- DAG[[2]]
  DAG <- DAG[[1]]
  
  tosaveTstats <- NULL
  DesEst <- NULL
  TStatistics <- NULL
  GIES.tosave <- NULL
  GES.tosave <- NULL

  
  GIES.score.tosave <- NULL
  GIES.score.tosave.All <- GES.score.tosave.All <- NULL

    
    
    
  ## Main for loop for the nrep iterations
  for(i in 1:nrep){
    
    ## Pick the i-th DAG and generate data
    MatDAG <- DAG[[i]]
    d <- GenData(DAG = DAG[[i]], n = n, interventions = interventions,permute = T,
                 IntMean = IntMean, IntVariance = NULL, NTot = NULL, distr = distr, SNR = SNR,
                 normalize = normalize, df = df, DoInt = DoInt, FunctionType = FunctionType)
    data <- d[[1]]
    invPermutation <- d[[6]]
    n <- d[[2]]
    interv <- d[[8]]
    n.sum <- n
    if(length(n)>1){
      for(q in 2:length(n)){
        n.sum[q] <- n.sum[q] + n.sum[q-1]
      }
    }
    
    CovMatAll <- cov(data)

    
    E1 <- NULL
    E2 <- NULL
    
    ## Compute the GIES estimates
    targets <- list(integer(0))
    targets <- append(targets,as.list(invPermutation[interv]))
    
    
    target.index <- rep(1,n[1])
    if(length(targets)>1){
      for(ii in 2:(length(targets))){
        target.index <- c(target.index,rep(ii,n[ii]))
      }
    }
    

  score.E1 <- new("GaussL0penIntScore", data, targets, target.index, intercept=T, lambda = lambda*log(nrow(data)))
  E1 <- gies(score.E1,verbose = F, phase = c("forward", "backward"), iterate = FALSE)
  GIES.tosave <- append(GIES.tosave,list(as(E1$essgraph,"matrix")))
    

    
  score.ges <- new("GaussL0penObsScore", data, intercept=T, lambda = lambda*log(nrow(data)))
  prov.ges <- ges(score = score.ges, verbose = FALSE, iterate = FALSE, phase = c("forward", "backward"))
  GES.tosave <- append(GES.tosave,list(as(prov.ges$essgraph,"matrix")))
    

  
  
  TS <- TS_input

  ### Compute the test statistics
  EstTstats <- TTestDesGraph(data = data[,invPermutation], n = n, interventions = interv, alpha = TS,
                               OptAlphaEqBeta = F, TStatBool = TRUE)
  
  
  DesEst <- append(DesEst,list(EstTstats[[1]]))
  TStatistics <- append(TStatistics,list(EstTstats[[3]]))

  
  
  
  interv2 <- invPermutation[interv]
  CV5FoldEstGIES <- CV5FoldEstGES <- CV5FoldInt <- CV5FoldCovMat <- NULL

  
  ## Compute the 5 outputs for the different algorithm used in the CV based risk estimator
  CV5ToSave <- NULL
  if(length(interv)>1){

    ### Below k-fold (5) CV
    
    int.length <- length(interv)
    fold.length <- max(1,int.length/5)
    cv5partition <- sample(x = 1:int.length, size = int.length, replace = FALSE)
    for(kk in 1:min(int.length,5)){
      if(kk==min(int.length,5)){
        k <- invPermutation[interv][sort(cv5partition[round(fold.length*(kk-1)+1):length(interv)])]
      }else{
        k <- invPermutation[interv][sort(cv5partition[round(fold.length*(kk-1)+1):round(fold.length*kk)])]
      }
      num.int <- which(invPermutation[interv]%in%k)
      num.int <- sort(num.int, decreasing = TRUE)
      
      dataMinus <- data
      for(ijk in 1:length(num.int)){
        kk.index <- num.int[ijk]
        dataMinus <- dataMinus[-((n.sum[kk.index]+1):n.sum[kk.index+1]),]
      }
      
      targets2 <- list(integer(0))
      targets2 <- append(targets2,as.list(setdiff(interv2,k)))
      
      target.index2 <- rep(1,n[1])
      if(length(targets2)>1){
        
        for(ijk in 1:length(interv)){
          if(length(intersect(invPermutation[interv][ijk],k))==0){
            
            target.index2 <- c(target.index2,rep(max(target.index2)+1,n[ijk+1]))
            
          }
        }
        
      }
      stopifnot(nrow(dataMinus) == length(target.index2))
      
      score.gies <- new("GaussL0penIntScore", dataMinus, targets2, target.index2, intercept=T, lambda = lambda*log(nrow(data)))


        E2 <- gies(score.gies,verbose = F, phase = c("forward", "backward"), iterate = FALSE)
      
      
      
      CovMat <- cov(dataMinus)
      CV5FoldEstGIES <- append(CV5FoldEstGIES, list(as(E2$essgraph,"matrix")))
      CV5FoldInt <- append(CV5FoldInt, list(k))
      CV5FoldCovMat <- append(CV5FoldCovMat, list(CovMat))
      
      score.ges <- new("GaussL0penObsScore", dataMinus, intercept=T, lambda = lambda*log(nrow(data)))

        E2 <- ges(score.ges,verbose = F, phase = c("forward", "backward"), iterate = FALSE)
      
      
      
      CV5FoldEstGES <- append(CV5FoldEstGES, list(as(E2$essgraph,"matrix")))

    }
    CV5ToSave <- list(CV5FoldEstGIES, CV5FoldEstGES, CV5FoldInt, CV5FoldCovMat)
    ### Start eta(iota)
    
  }
    
  }
    
  
  
  sysInfo <- sessionInfo()
  
  if (length(n)==1){
    n <- c(n,n)
  }
  
  name <- paste("SimROCData,",Simindex, "p=",p, "n=",n[2], "ENS=", ENS,"interventions=",interventions,
                "nrep=", nrep,"distr=",distr,"lambda=",lambda, "SNR=",SNR.input,
                "IntMean=",IntMean, "df=",df,"seed=",seed,sep = "")
  tosave <- list(DAG=DAG, TStatistics=TStatistics, DesEst=DesEst, GIES.tosave=GIES.tosave, GES.tosave = GES.tosave, GIES.score.tosave.All=GIES.score.tosave.All, GES.score.tosave.All=GES.score.tosave.All, 
                 name=name, invPermutation = invPermutation,
                 CovMatAll = CovMatAll, CV5ToSave = CV5ToSave, seed = seed,
                 DoInt = DoInt, fn = fn, FunctionType = FunctionType, sysInfo = sysInfo)

  setwd("~/")

  save(tosave, file = name, version = 2)
}



