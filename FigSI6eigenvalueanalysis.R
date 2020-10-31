############# SI FIG.6  ###############

rm(list=ls())

#packages
library("expm", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

######## GLOBAL VARIABLES ########################

# memory range
mem.range <- 1:31

# OFF birth rate for untreated
b.r.untreat <- 1

# OFF birth rate for treated
b.r.treat <- 1

# OFF death rate for untreated
d.r.untreat <- 0.98

# OFF death rate for treated
d.r.treat <- 1.02

# mu for untreated
mu.untreat <- 0

# mu for treated
mu.treat <- 0.2

# epsilon 
eps <- 0.25

# length of environment
env.duration <- 15

# environmental sequences, 1 means is no treatment, 0 is treatment.
env.scheme <- sapply(0:6, function(x) {y<-rep(1,6); y[(c(0:x))]<-0; return(y)})
env.scheme <- rbind(env.scheme,sapply(6:0, function(x) {y<-rep(0,6); y[(c(0:x))]<-1; return(y)}))
env.scheme <- env.scheme[-6,]
env.scheme <- env.scheme[,(ncol(env.scheme):1)]


# # random environmental sequences, 1 means is no treatment, 0 is treatment.
# env.scheme <- sample(c(0,1), size = 144, replace = TRUE)
# env.scheme <- matrix(env.scheme,ncol=12,nrow=12)




####### MATRIX MODELS ##########

# function to create matrix model
ppm <- function(OFFbirth, # a scalar birth rate equal for all compartments
                OFFdeath, # a scalar death rate equal for all compartments
                ONbirth,  # a scalar birth rate equal for all compartments
                ONdeath,  # a scalar death rate equal for all compartments
                my.mu,       # the scalar mu
                my.eps,      # a scalar leaching rate through 'on' compartments
                mem       # number of 'on' compartments
) {
  nr.compartments <- mem+1
  b.r <- c(OFFbirth, rep(ONbirth, mem))
  d.r <- c(OFFdeath, rep(ONdeath, mem))
  t.r <- c(my.mu,rep(my.eps,mem))
  
  if (nr.compartments > 2) {
    M <- diag(b.r-d.r-t.r)
    M <- M + rbind(0,cbind(diag(t.r[-(nr.compartments)]),0))
    M[1,nr.compartments] <- my.eps}
  
  else {
    M <- diag(b.r-d.r-t.r)
    M[2,1] <- my.mu
    M[1,2] <- my.eps
  }
  return(M)
}

# Eigenvalues of the projection matrix ordered in decreasing magnitude of the real part
lambdas <- function(ppm) (sort(Re(eigen(ppm, only.values = TRUE)[[1]]), decreasing = TRUE))[1]

####### FITNESSES ##########


fitnesses <- sapply(mem.range, function(x) {
                    untreated <- ppm(OFFbirth = b.r.untreat, OFFdeath = d.r.untreat, ONbirth = 0, ONdeath = 0, my.mu = mu.untreat, my.eps = eps, mem = x)
                    treated <- ppm(OFFbirth = b.r.treat, OFFdeath = d.r.treat, ONbirth = 0, ONdeath = 0, my.mu = mu.treat, my.eps = eps, mem =x)
                      sapply(1:ncol(env.scheme), function(y) {
                      envs <- env.scheme[ ,y]
                      time <- 1
                      n0 <- t(t(c(1000,rep(0,x))))
                      my.mat <-diag(x+1)
                      while (time <= length(envs)) {
                       my.mat <- expm(env.duration * if(envs[time] == 1) {treated} else {untreated}) %*%  my.mat
                       time <- time+1
                      }
                      eig.lambda <- lambdas(my.mat)
                      nsemifinal1 <- (my.mat %^% (1-1)) %*% n0
                      nfinal1 <- my.mat %*% nsemifinal1
                      est.lambda1 <- sum(nfinal1) / sum(nsemifinal1)
                      nsemifinal3 <- (my.mat %^% (3-1)) %*% n0
                      nfinal3 <- my.mat %*% nsemifinal3
                      est.lambda3 <- sum(nfinal3) / sum(nsemifinal3)
                      nsemifinal10 <- (my.mat %^% (10-1)) %*% n0
                      nfinal10 <- my.mat %*% nsemifinal10
                      est.lambda10 <- sum(nfinal10) / sum(nsemifinal10)
                      return(c(est.lambda1=est.lambda1, est.lambda3=est.lambda3,est.lambda10=est.lambda10,eig.lambda=eig.lambda))
                    }, simplify = TRUE)
                    }
                    )
fitnesses

par(mfrow=c(7,4),mar=(c(1.9, 2, 1, 2.1)))
sapply(1:7, function(x) {
plot(mem.range,fitnesses[ c(TRUE,FALSE,FALSE,FALSE), ][x,], xlab = "", ylab = "",pch=16)
plot(mem.range,fitnesses[ c(FALSE,TRUE,FALSE,FALSE), ][x,], xlab = "", ylab = "",pch=16)
plot(mem.range,fitnesses[ c(FALSE,FALSE,TRUE,FALSE), ][x,], xlab = "", ylab = "",pch=16)
plot(mem.range,fitnesses[ c(FALSE,FALSE,FALSE,TRUE), ][x,], xlab = "", ylab = "",pch=16)
},simplify = TRUE)