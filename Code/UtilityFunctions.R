
library(Matrix)
library(ggm)
library(Rfast)
library(pcalg)


### Generate nrep pxp random matrices (random DAGs). All DAGs will be in the topological order, that is, matrices are upper triangular.
GenMatDAGs <- function(nrep=100, p=20, prob=NULL, method = "er", par1=NULL, par2=NULL,
                       lowerbound = 1, upperbound = 3, negValues=TRUE, fn){
  
  if(is.null(prob)){
    prob <- p/(p*(p-1)/2)
  }
  
  ENS <- prob * (p-1)
  
  lB = lowerbound
  uB = upperbound
  
  
  
  ## Generate random DAG
  wFun <- function(m, lB, uB){runif(m,lB,uB)*(-1)^(rbinom(n = 1, size = 1, prob = 0.5)*negValues)}
  DAG <- randDAG(n = p, d = ENS, method=method, weighted = TRUE, wFUN = list(wFun,lB,uB), par1=NULL, par2=NULL)
  
  ### Transform to upper triangular
  DAG <- as(DAG,"matrix")
  p_prov <- PermuteToUpperTri(G = DAG)
  DAG <- DAG[p_prov,p_prov]
  
  FunctionType <- matrix(data = "N", nrow = p, ncol = p)
  FunctionType[which(DAG!=0, arr.ind = T)] <- sample(x = c("linear", "sigmoid"), size = sum(DAG!=0), replace = T)
  if(fn=="linear"){
    FunctionType[FunctionType=="sigmoid"] <- "linear"
  }
  if(fn=="sigmoid"){
    FunctionType[FunctionType=="linear"] <- "sigmoid"
  }
  
  
  ## Transform random DAG into adj. matrix
  # MatDAG <- list(as(DAG,"matrix"))
  MatDAG <- list(DAG)
  
  
  if(nrep>1){
    for(i_nrep in 2:nrep){
      
      ## Generate random DAG
      DAG <- randDAG(n = p, d = ENS, method=method, weighted = TRUE, wFUN = list(wFun,lB,uB), par1=NULL, par2=NULL)
      
      DAG <- as(DAG,"matrix")
      p_prov <- PermuteToUpperTri(G = DAG)
      DAG <- DAG[p_prov,p_prov]
      
      ## Transform random DAG into adj. matrix
      # MatDAG <- append(MatDAG,list(as(DAG,"matrix")))
      MatDAG <- append(MatDAG,list(DAG))
      
      
      
    }
  }
  
  return(list(MatDAG,FunctionType))
  
}






### Given a matrix (DAG) generates n i.i.d. observation from it
GenDataFromMat <- function(Mat,n=1000,means=NULL,var.bounds=NULL,scale=F, distr = "Normal", df=4,
                           SNR=NULL, normalize = FALSE, FunctionType){
  
  ### Assume Mat is upper triangular
  if(sum(diag(abs(Mat)))>0){
    stop("Mat has non-zero entries on the diagonal")
  }
  
  if(sum(abs(Mat[lower.tri(Mat)]))>0){
    stop("Mat is not upper triangular")
  }
  
  ### Initialize
  p<-ncol(Mat)
  X<-matrix(data = 0,nrow = n,ncol = p)
  err<-matrix(data = 0,nrow = n,ncol = p)
  
  
  ### Set the variance and mean of the observations (distribution)
  if(is.null(var.bounds)){
    var.bounds <- cbind(rep(0.5,p),rep(1.5,p))### 0.5 1
    if(!is.null(SNR)){
      if(normalize){
        w <- 1/(SNR+1)
        var.bounds <- cbind(rep(w,p),rep(w,p))
      }
      
    }
  }
  
  
  
  if(is.null(means)){
    means <- numeric(p)
  }
  
  var <- numeric(p)
  for(i in 1:p){
    var[i] <- runif(n = 1,min = var.bounds[i,1],max = var.bounds[i,2])
  }
  
  
  ## Can control the SNR ratio.
  ## Normalize is used to set the variance at each node of the graph to 1.
  ## This menas that the sum of the variance of the signal and those of the noise term has to be 1
  ## Hence, the variance of the noise cannot be larger than 1.
  if(!is.null(SNR)){
    if(normalize){
      if(max(var)>=1){
        stop("Variance too large for SNR = TRUE")
      }
    }
  }
  
  
  ## Generate observations from a distribution
  for(i in 1:p){
    
    if(distr == "Cauchy"){
      X <- matrix(data = rcauchy(n = n*p, location = means[i], scale = 2), nrow = n, ncol = p, byrow = FALSE)
    }else if(distr == "Uniform"){
      X <- matrix(data = runif(n = n*p, min = (-sqrt(3*var[i]) + means[i]), max = (sqrt(3*var[i]) + means[i])), nrow = n, ncol = p, byrow = FALSE)
    }else if(distr == "t"){
      X[,i] <- rt(n = n, df = df) + means[i]
      stopifnot(df>2)
      X[,i] <- X[,i]/sqrt(df/(df-2))*sqrt(var[i])
    }else if(distr == "lognormal"){
      X[,i] <- rlnorm(n = n) - exp(0.5) + means[i]
      X[,i] <- X[,i]/sqrt((exp(1)-1)*exp(1))*sqrt(var[i])
    }else {
      X[,i] <- rnorm(n = n,mean = means[i], sd = sqrt(var[i]))
    }
  }
  ## Noises is probably not used anymore
  Noises <- X
  
  ## For the control the the SNR we use C. Heinze-Deml, M.H. Maathuis and N. Meinshausen (2018). Causal structure learning. Annual Review of Statistics and Its Applications
  ## This paper is for the normalized version
  if(!is.null(SNR)){
    
    if(!normalize){
      
      for(i in 1:p){
        ## We sum over the incoming variables but without the noise. So this is the signal
        prov <- 0
        if(sum(Mat[,i]!=0)==0){
          prov <- X[,1]*Mat[1,i]
          varprov <- var(prov)
        }else{
          for(j in 1:(i)){
            if(FunctionType[j,i]=="sigmoid"){
              prov <- prov + Sigmoid(X[,j])*Mat[j,i]
            }else{
              prov <- prov + X[,j]*Mat[j,i]
            }
          }
          varprov <- var(prov)
        }
        ## varprov = 0 only if there are no incoming variables (constant vector 0)
        if(varprov!=0){
          
          Mat[,i] <- Mat[,i]*sqrt(SNR*var[i]/varprov)
          
        }
        
        for(j in 1:(i)){
          
          if(FunctionType[j,i]=="sigmoid"){
            X[,i]<-X[,i]+Sigmoid(X[,j])*Mat[j,i]
          }else{
            X[,i]<-X[,i]+X[,j]*Mat[j,i]
          }
        }
      }
      return(list(X,var,Noises,Mat))
      
    }else{
      ## Same but with normalization
      for(i in 1:p){
        
        prov <- 0
        if(sum(Mat[,i]!=0)==0){
          prov <- X[,1]*Mat[1,i]
          varprov <- var(prov)
        }else{
          for(j in 1:(i)){
            if(FunctionType[j,i]=="sigmoid"){
              prov <- prov + Sigmoid(X[,j])*Mat[j,i]
            }else{
              prov <- prov + X[,j]*Mat[j,i]
            }
          }
          varprov <- var(prov)
        }
        
        
        varprov <- var(prov)
        if(varprov!=0){
          
          w2 <- sqrt((1- var[i])/varprov)
          Mat[,i] <- Mat[,i]*w2
          
          for(j in 1:(i)){
            if(FunctionType[j,i]=="sigmoid"){
              X[,i]<-X[,i]+Sigmoid(X[,j])*Mat[j,i]
            }else{
              X[,i]<-X[,i]+X[,j]*Mat[j,i]
            }
          }
          
          
        }else{ # varprov is 0
          for(j in 1:(i)){
            if(FunctionType[j,i]=="sigmoid"){
              X[,i]<-X[,i]+Sigmoid(X[,j])*Mat[j,i]
            }else{
              X[,i]<-X[,i]+X[,j]*Mat[j,i]
            }
            if(Mat[j,i]!=0){
              stop("Coefficient should be zero")
            }
          }
          ## Since we normalize it is important to know the variance of the noise term which in the case of
          ## no parents it has to be set to 1.
          X[,i] <- X[,i]/sqrt(var[i])
          var[i] <- 1
        }
        
      }
      return(list(X,var,Noises,Mat))
    }
    
    
  }
  
  
  inv <- solve(diag(p) - Mat)
  X <- X %*% inv

  
  return(list(X,var,Noises,Mat))
}





## Generate obsetvational and interventional data using GenDataFromMat
GenData <- function(DAG, n=1000, oracle=F, hidden=0, interventions=0, permute = T,
                    IntMean = NULL, IntVariance = NULL, NTot=NULL, distr = "Normal",
                    SNR=NULL, normalize = FALSE, df, DoInt, FunctionType){
  
  ## Hidden: probability of a variable to be unobserved (not used at the moment)
  ## interventions: probability of a variable to be intervened on among the observed variables
  ## For now: Observational data are always present

  
  if(!is.null(SNR)){
    if(DoInt){
      IntVariance <- 1
    }else{
      IntVariance <- 1/(SNR+1)
    }
  }
  
  p <- ncol(DAG)
  data <- NULL
  
  permutation <- sample(x = p)
  if(!permute){
    permutation <- 1:p
  }
  
  
  ## Generate the observational data
  dataprov <- GenDataFromMat(Mat = DAG, n = n[1], distr = distr, SNR = SNR, normalize = normalize, df=df, FunctionType = FunctionType)
  if(!is.null(SNR)){
    DAG <- dataprov[[4]]
  }
  variance <- dataprov[[2]]
  
  
  data <- rbind(data,dataprov[[1]])
  if(oracle){
    CovDataArray.oracle <- solve(diag(p)-t(DAG))%*%diag(variance)%*%solve(diag(p)-DAG)
    
    empiricalCovArray <- cov(data)
    e <- eigen(empiricalCovArray)
    sqrt.cov.mat.old <- e$vectors%*%sqrt(diag(e$values))
    e <- eigen(CovDataArray.oracle)
    sqrt.cov.mat.new <- e$vectors%*%sqrt(diag(e$values))
    
    
    
    data <- t(sqrt.cov.mat.new%*%solve(sqrt.cov.mat.old,t(data)))
  }
  ### We permute the comumns in order to avoid possible use of the ordered variables
  data <- data[,permutation]
  
  interventions.index <- hidden.index <- NULL
  
  hidden.bool <- rbinom(n = p, size = 1, prob = hidden)
  hidden.index <- which(hidden.bool==1)
  if(length(hidden.index)>0){
    data <- data[,-hidden.index]
  }
  
  ## Allocate the interventions
  interventions.bool <- rbinom(n = p, size = 1, prob = interventions) * (1 - hidden.bool)
  
  
  for(i in 1:p){
    
    ## For each node where we have an intervention we compute the new DAG and generate data for it
    ## according to the input variance and mean
    if(interventions.bool[i]){
      
      
      if(is.na(n[sum(interventions.bool[1:i])+1])){
        n <- c(n,n[1])
      }
      
      MatDAGint <- DAG
      varianceInt <- variance
      if(DoInt){
        MatDAGint[,i] <- 0
      }
      
      
      
      
      means <- numeric(p)
      
      if(is.null(IntMean)){
        means[i] <- runif(n = 1,min = 1,max = 3)
      }else{
        means[i] <- IntMean
      }
      
      if(is.null(IntVariance)){
        varianceInt[i] <- runif(n = 1,min = 0.5,max = 1)
      }else{
        varianceInt[i] <- IntVariance      
      }
      
      
      dataprov <- GenDataFromMat(Mat = MatDAGint, n = n[sum(interventions.bool[1:i])+1],means = means, var.bounds = cbind(varianceInt,varianceInt), distr = distr, SNR = NULL, normalize = normalize, df=df, FunctionType = FunctionType)
      
      
      dataToTransform <- dataprov[[1]]
      if(oracle){
        
        
        CovDataArray.oracle <- solve(diag(p)-t(MatDAGint))%*%diag(varianceInt)%*%solve(diag(p)-MatDAGint)
        
        empiricalCovArray <- cov(dataToTransform)
        e <- eigen(empiricalCovArray)
        sqrt.cov.mat.old <- e$vectors%*%sqrt(diag(e$values))
        e <- eigen(CovDataArray.oracle)
        sqrt.cov.mat.new <- e$vectors%*%sqrt(diag(e$values))
        
        
        
        dataToTransform <- t(sqrt.cov.mat.new%*%solve(sqrt.cov.mat.old,t(dataToTransform)))
      }
      
      interventions.index <- c(interventions.index, i)
      
      dataToTransform <- dataToTransform[,permutation]
      
      if(length(hidden.index)>0){
        dataToTransform <- dataToTransform[,-hidden.index]
      }
      
      data <- rbind(data,dataToTransform)
      
    }  
  }
  
  
  ## Now we have the observational and interventional data
  
  if(length(hidden.index)>0){
    for(qq in seq_along(hidden.index)){
      qqq <- length(hidden.index) + 1 - qq
      
      a <- hidden.index[qqq]
      
      permutation <- permutation[-which(permutation==a)]
      permutation[which(permutation>a)] <- permutation[which(permutation>a)] - 1
      
      interventions.index[which(interventions.index>a)] <- interventions.index[which(interventions.index>a)] - 1
    }
  }
  
  interventions.bool <- interventions.bool[which(hidden.bool==0)][permutation]
  
  
  ## Compute the place of the interventions for the permutated DAG
  perm_p <- ncol(data)
  inv_perm <- invPerm(permutation)
  interventions.index <- interventions.bool[inv_perm]*(1:perm_p)#[inv_perm]
  interventions.index <-interventions.index[interventions.index!=0]
  
  
  n <- n[1:(length(interventions.index) + 1)]
  
  
  # ## NTot is used if instead to decide on the number of observations per interventions
  # ## we fix the total number of observations (obs + int)
  # ## in this case we need to consider only part of the data
  # 
  # ## Not used anymore
  # if(!is.null(NTot)){
  #   if(NTot < sum(n)){
  #     
  #     rat <- NTot/sum(n)
  #     
  #     n <- round(n * rat)
  #     
  #     n.sum <- n
  #     if(length(n)>1){
  #       for(q in 2:length(n)){
  #         
  #         n.sum[q] <- n.sum[q] + n.sum[q-1]
  #         
  #       }
  #     }
  #     
  #     dat_prov <- data[1:n[1],]
  #     
  #     
  #     if(length(n)>1){
  #       for(i in 2:length(n)){
  #         
  #         dat_prov <- rbind(dat_prov, data[(n.sum[i-1]+1):n.sum[i],])
  #         
  #       }
  #     }
  #     
  #     data <- dat_prov
  #   }
  # }
  
  
  return(list(data,n,variance,hidden.bool,interventions.bool,inv_perm,permutation,interventions.index, DAG))
  
}




RiskEstimation <- function(G1,G2,i,AlsoObs=T, TCG2=F, TCG1=F, type = "Jaccard"){
  
  ## This plays a role when we have to compute the distance between nodes.
  ## The function is therefore not symmetric w.r.t G1 and G2
  
  
  p <- ncol(G1)
  
  
  
  ## TCG1 and 2 are always FALSE now, we input directly the "right" graph
  if(TCG1){
    T1 <- transClos(G1) 
  }else{
    T1 <- G1
  }
  
  if(TCG2){
    T2 <- transClos(G2)
  }else{
    T2 <- G2
  }
  
  
  ## Compute the jaccard index
  if(type == "Jaccard"){
    Jaccard.index <- 0
    
    h <- i
    T1 <- transClos(T1)
    DescSetA <- which(T1[h,]!=0)
    DescSetB <- which(T2[h,]!=0)
    
    if(length(DescSetA)+length(DescSetB)!=0){
      Jaccard.index <-Jaccard.index + 1 - length(intersect(DescSetA,DescSetB))/(length(DescSetA)+length(DescSetB) - length(intersect(DescSetA,DescSetB)))
    }else{
      Jaccard.index <- Jaccard.index + 0
    }
    
    
    
    ## The Frobenius norm is not used anymore, but we keep it here
    T1[-i,] <- 0
    T2[-i,] <- 0
    
    
    
    FNorm <- norm(T1-T2,type = "F")
    return(list(FNorm, Jaccard.index))
  }
  
  ## Compute the loss with power decay
  if(type=="Power"){
    
    h <- i
    
    ## We set the distance between nodes that cannot be reached to infinity and then compute the shortest path
    ## between all nodes
    T1floyd <- T1
    T1floyd[T1floyd==0] <- Inf
    a <- floyd(T1floyd)
    a[a==Inf] <- 0
    
    
    ## Now we compute the transitive closure of T1 so that we can compare the descendants more easily
    T1 <- transClos(T1)
    DescSetAprov <- which(T1[h,]!=0)
    DescSetB <- which(T2[h,]!=0)
    
    UnionAB <- union(DescSetAprov, DescSetB)
    
    DescSetA <- setdiff(DescSetAprov, DescSetB)
    DescSetB <- setdiff(DescSetB, DescSetAprov)
    
    
    
    stopifnot(sum(a[h,DescSetA]) != Inf)
    stopifnot(sum(a[h,DescSetB]) != Inf)
    
    powerloss <- 0
    
    if(length(DescSetA)!=0){
      for(ides in seq_along(DescSetA)){
        powerloss <- powerloss + 2/2^(a[h,DescSetA[ides]])
      }
    }
    if(length(DescSetB)!=0){
      for(ides in seq_along(DescSetB)){
        powerloss <- powerloss + 1### Maybe something better later!
      }
    }
    if(length(UnionAB)>0){
      powerloss <- powerloss / length(UnionAB)
    }
    return(list(powerloss, powerloss))
    
  }
  
  
  ## Similar but with polynomial decay
  if(type=="Poly"){
    
    h <- i
    
    T1floyd <- T1
    T1floyd[T1floyd==0] <- Inf
    a <- floyd(T1floyd)
    a[a==Inf] <- 0
    
    T1 <- transClos(T1)
    DescSetAprov <- which(T1[h,]!=0)
    DescSetB <- which(T2[h,]!=0)
    
    UnionAB <- union(DescSetAprov, DescSetB)
    
    
    DescSetA <- setdiff(DescSetAprov, DescSetB)
    DescSetB <- setdiff(DescSetB, DescSetAprov)
    
    stopifnot(sum(a[h,DescSetA]) != Inf)
    stopifnot(sum(a[h,DescSetB]) != Inf)
    
    polyloss <- 0
    
    if(length(DescSetA)!=0){
      for(ides in seq_along(DescSetA)){
        polyloss <- polyloss + 1/(a[h,DescSetA[ides]])
      }
    }
    if(length(DescSetB)!=0){
      for(ides in seq_along(DescSetB)){
        polyloss <- polyloss + 1### Maybe something better later!
      }
    }
    
    if(length(UnionAB)>0){
      polyloss <- polyloss / length(UnionAB)
    }
    return(list(polyloss, polyloss))
    
  }
  

}






PermuteToUpperTri <- function(G){
  
  p <- ncol(G)
  
  if(p==1){
    return(1)
  }
  G <- (G!=0)*1
  
  colsum <- apply(X = G,MARGIN = 2, FUN = sum)
  colsum2 <- colsum
  
  perm <- order(colsum)
  
  
  ## This is now the candidate permutation to transform G
  perm <- (1:p)[perm]
  
  k=0
  ## If the candidate is good, the sum of its lower trinagular part should be zero
  while(sum(lower.tri(G[perm,perm]) * G[perm,perm])!=0){
    
    ## If the iteration ran over all nodes, it is not possible to transform the matrix
    if(k == (p-1)){
      stop("Cannot transform to upper triangular")
    }
    
    ## The idea is that if the sum over the columns is different for each column, then we are done in one step
    ## If this is not the case, this menas that at least two columns have the same sum
    
    ## We start from the first column and add 1, if this column had a lower sum, then it will now have again a lower (or equal) sum
    ## In any case, this does not change its position when using order()
    ## If it had the same sum as the second (and maybe third..) column, they will flip
    k <- k+1
    perm2 <- perm
    
    G2 <- G[perm2,perm2]
    
    colsum <- apply(X = G2[k:p,],MARGIN = 2, FUN = sum)
    
    colsum[k] <- colsum[k]+1
    
    perm2 <- order(colsum)
    perm <- perm[perm2]
    
  }
  
  return(perm)
}



## Compute the T statistics
TTestDesGraph <- function(data, n, interventions, alpha, OptAlphaEqBeta=F, TStatBool=F){
  
  stopifnot(alpha>=0)
  p <- ncol(data)
  Des <- matrix(data = 0, nrow = p, ncol = p)
  TTests <- array(data = 0, dim = c((p-1)*length(interventions), 3)) ## Last two colums are parents and children
  ### Create cumulative sum of index n, if n is NULL, then all same size
  
  if(length(interventions)>0){
    
    
    n.sum <- n
    
    
    if(length(n)>1){
      for(i in 2:length(n)){
        
        n.sum[i] <- n.sum[i] + n.sum[i-1]
        
      }
    }
    
    ## Observatonal data
    dataObs <- data[1:n.sum[1],]
    
    
    
    for(i in 1:length(interventions)){
      
      ## Interventional data (for specific intervention)
      dataInt <- matrix(data = data[(n.sum[i] + 1) : n.sum[i+1],], ncol = p, byrow = FALSE)
      
      PossDes <- setdiff(1:p, interventions[i])
      
      if(length(PossDes)>1){
        
        for(k in 1:length(PossDes)){
          
          target <- PossDes[k]
          
          # ## Not used, always FALSE
          # if(VarAdaptive){
          #   
          #   n.alphabetaOpt <- n[1]
          #   sp <- (sd(dataInt[,target])+ sd(dataObs[,target]))/sqrt(2)
          #   meanshift <- 2 ## Can be adapted to first check for the minimal meanshift! The interventions are known
          #   
          #   alpha <- exp(-seq(from=1, to = 20, length.out = 1000))
          #   
          #   alpha2 <- alpha
          #   
          #   
          #   val <- qt(p = (1-alpha2/2), df = 2*(n.alphabetaOpt-1))
          #   pow <- 1-pt(q = val, df = 2*(n.alphabetaOpt-1), ncp = meanshift*sqrt(n.alphabetaOpt)/(sp*sqrt(2)))
          #   # alpha[min(which(pow >0.99))]
          #   
          #   if(OptAlphaEqBeta){
          #     alpha <- alpha[which.min(abs(1-pow-alpha))]
          #   }else{
          #     alpha <- alpha[which.max(pow-alpha)]
          #   }
          #   
          # }
          

          
          ## Until here very similar to the first version
          ## Now we save the values of the tstatistics with the corresponding parent and child
          
          TTests[(i-1)*length(PossDes)+k,1] <- t.test(x = dataObs[,target], y = dataInt[,target], alternative = "two.sided", paired = F, var.equal = T)$statistic
          TTests[(i-1)*length(PossDes)+k,2] <- interventions[i]
          TTests[(i-1)*length(PossDes)+k,3] <- target
          

        }
        
      }
      
    }
  }
  
  
  ## We construct a matrix where all edges are present in each direction and the number in the (adj.) matrix
  ## correspond to the order of entrance of that edge in the final matrix
  orderTstats <- order(TTests[,1])
  for(ee in 1:nrow(TTests)){
    
    Des[TTests[orderTstats[ee],2], TTests[orderTstats[ee],3]] <- ee
    
  }  
  
  
  Des2 <- Des
  Des3 <- matrix(data = 0, nrow = p, ncol = p)
  ee <- 1
  while((abs(TTests[orderTstats[ee],1]) > alpha) && (ee<=(p-1)*length(interventions)) ){
    
    ## Insert the edges step-by-step in order to be able to check for cycles
    Des3[Des2==ee] <- 1
    ee <- ee+1
  }
  
  return(list(Des3,alpha,TTests))
  
  
}


Sigmoid <- function(x){
  res <- 10/(1+exp(-0.65*x)) - 5
  return(res)
}