### Function used to run the simulation and sample the parameter space


library(reshape2)
library(ggplot2)
library(ggm)
library(doParallel)


source("~/UtilityFunctions.R")
source("~/CompGraphs.R")



ToAdd <- 0
p <- c(25,25)
n <- c(100,100)
nrep <- 1
interventions2 <- c(0.5,0.5)
lambda <- c(0.5)
ENS <- c(1.5,2.5)
SNR.input <- c(1,3,5)
IntMean <- c(5,5)
distr.input <- c("Normal", "t", "lognormal")
DoInt.input <- c(TRUE, FALSE)
fn.input <- c("linear", "sigmoid", "linsig")
method="er"
par1 <- NULL
par2 <- NULL
df <- c(5,5)
TS_input = log(n)
normalize <- TRUE


no_cores <- 4

registerDoParallel(cores = no_cores)
NotNeeded <- foreach(i=1:50) %dopar%
  {
    

    seed = (i+ToAdd)
    set.seed(seed)
    
    p.arg <- sample(x = p, size = 1)
    n.arg <- sample(x = n, size = 1)
    ENS.arg <- sample(x = ENS, size = 1)
    interventions.arg <- sample(x = interventions2, size = 1)
    SNR.input.arg <- sample(x = SNR.input, size = 1)
    IntMean.arg <- sample(x = IntMean, size = 1)
    Simindex <- (i+ToAdd)
    df.arg <- sample(x = df, size = 1)
    DoInt.arg <- sample(x = DoInt.input, size = 1)
    distr.arg <- sample(x = distr.input, size = 1)
    fn.arg <- sample(x = fn.input, size = 1)
    n.arg2 <- c(max(100,n.arg),rep(n.arg,p.arg))
    
    

    
    CompGraphs(p = p.arg,
                n = n.arg2,
                ENS = ENS.arg,
                nrep = nrep,
                interventions = interventions.arg,
                distr = distr.arg,
                df=df.arg,
                method = method, par1 = par1, par2 = par2,
                SNR.input = SNR.input.arg,
                IntMean = IntMean.arg,
                TS_input = TS_input,
                lambda = lambda,
                Simindex = Simindex,
                max.kk = max.kk, normalize = normalize,
                seed = (i+ToAdd),
                DoInt = DoInt.arg,
                fn = fn.arg)
  }


