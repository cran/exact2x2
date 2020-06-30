# Calculation to show that for all n=1,2,...,100
# the coverage of the 95% confidence interval for the 
# Delta = unconditional difference =
# proportion positive minus proportion negative
# where prop pos + prop neg + prop zero = 1
# Also Delta = E(S), where S=sign(x)
# 
# We can write Delta = theta (2*beta-1)
#   where theta = Pr[S=1 or S=-1]
#         beta = Pr[S=1 | S !=0 ]
# We calculated over the grid
#   betas <- seq(0, 1, by = .01)
#   thetas <- seq(0, 1, by = .01)
# For each n, and each grid point we calculate 
#     lowerError = Pr[ L(X) > Delta]
#     upperError = Pr[ U(X) < Delta]
# where [L(X),U(X)] is the 95% confidence for X
# Then we find the maximum of the 2 errors over the grid for each n
# We want the maximum of each set of errors to be less than 2.5% for the 95% central CI
#
# Tests show that for all n=1,2,..100 the errors are less than 2.5% (within rounding error)
# The only n values with greater than 2.5% error are  n=71 and n=72,
#  but the errors are likely due to rounding error in the 
# numeric integration. Those errors  are  0.02503121 for n=72 (beta=0, theta=0.96) or (beta=1, theta=0.96)
# and 0.025000563 for n=71 (beta=0, theta=0.21) or (beta=1, theta=0.21)
#
# use exacgt2x2 version 1.6.4 or greater
library(exact2x2)
library(ggplot2)
library(grid)
library(gridExtra)

mcnemarSim <- function(n, nmc = 0, plots = TRUE) {
  
  enum <- expand.grid(0:n, 0:n)
  enum <- enum[enum$Var1 >= enum$Var2, ]
  names(enum) <- c("m", "x")
  enum <- enum[order(enum$m), ]
  enum$Lower <- numeric(nrow(enum))
  enum$Upper <- numeric(nrow(enum))
  
  for (i in 1:nrow(enum)) {
    if(nmc == 0){
      CI <- mcnemarExactDP(n = n, m = enum$m[i], x = enum$x[i])$conf.int
    } else {
      CI <- mcnemarExactDP(n = n, m = enum$m[i], x = enum$x[i], nmc = nmc)$conf.int
    }
    
    enum$Lower[i] <- CI[1]
    enum$Upper[i] <- CI[2]
    
  }
  
  betas <- seq(0, 1, by = .01)
  thetas <- seq(0, 1, by = .01)
  
  results <- expand.grid(betas, thetas)
  names(results) <- c("Beta", "Theta")
  results$LowErrorProb <- numeric(nrow(results))
  results$HighErrorProb <- numeric(nrow(results))
  
  count <- 1
  for (t in thetas) {
    for (b in betas) {
      
      delta <- t * (2 * b - 1)
      
      loopframe <- enum
      loopframe$LowError <- ifelse(delta < loopframe$Lower, 1, 0)
      loopframe$HighError <- ifelse(delta > loopframe$Upper, 1, 0)
      loopframe$Prob <-
        dbinom(loopframe$m, n, t) * dbinom(loopframe$x, loopframe$m, b)
      
      results$LowErrorProb[count] <-
        sum(loopframe$LowError * loopframe$Prob)
      
      results$HighErrorProb[count] <-
        sum(loopframe$HighError * loopframe$Prob)
      
      count <- count + 1
      
    }
  }
  
  if(plots){
    lowheatmap <-
      ggplot(results, aes(x = Beta, y = Theta, color = LowErrorProb, fill = LowErrorProb)) +
      geom_tile() +
      scale_fill_gradient2(low = "white", high = "black") +
      scale_color_gradient2(low = "white", high = "black") +
      ggtitle("Lower Error Probabilities")
    
    highheatmap <-
      ggplot(results, aes(x = Beta, y = Theta, color = HighErrorProb, fill = HighErrorProb)) +
      geom_tile() +
      scale_fill_gradient2(low = "white", high = "black") +
      scale_color_gradient2(low = "white", high = "black") +
      ggtitle("Upper Error Probabilities")
    
    grid.arrange(lowheatmap, highheatmap, ncol = 2)
  }
  
  return(results)
  
}

# Example: predict how much time the simulation will take
# run several times....
#t0<-proc.time()
#mc <- mcnemarSim(11, plots = TRUE)
#mc <- mcnemarSim(21, plots = TRUE)
# n=26 gives graph in paper
#mc <- mcnemarSim(26, plots = TRUE)
#dev.print(pdf,file="ConfIntErrors26.pdf")
#mc <- mcnemarSim(40, plots = TRUE)
#mc <- mcnemarSim(60, plots = TRUE)
#mc <- mcnemarSim(100, plots = TRUE)
#t1<-proc.time()
#sec<-(t1-t0)[1]

# Results of the times
#n<-c(11,21,26,40,60,100)
#sec<-c(9.02,32.03,49.36,120.86,278.6,789.12)
#n2<- n^2
#lout<-lm(sec~n+n2)
#coef(lout)
#     (Intercept)           n          n2 
#       2.25685660 -0.30005320  0.08169258 
#N<-1:100
#sec100<-predict(lout,newdata=data.frame(n=N,n2=N^2))
#plot(1:100,sec100,type="l",lwd=2,log="")
#points(n,sec,cex=2,pch=16)

# Predict all n=1:100, in hours
# sum(sec100)/(60*60)
##         7.31975  hours

# Run below to test values of n from 1 to 100
# for all valid values of m and x.

mcList <- vector("list", length = 100)
for(i in 1:100){
  
  print(i)
  mcList[[i]] <- mcnemarSim(i, plots = FALSE)
  
}

loerrorVec <- hierrorVec<-numeric(100)
for(i in 1:100){
  loerrorVec[i] <- max(mcList[[i]]$LowErrorProb)
  hierrorVec[i] <- max(mcList[[i]]$HighErrorProb)
} 
max(loerrorVec)
max(hierrorVec)
names(loerrorVec)<-names(hierrorVec)<-1:100
loerrorVec
hierrorVec

plot(1:100,loerrorVec,type="l",ylim=c(0.02,0.0255))

m72<-mcList[[72]]
m72[m72$LowErrorProb>0.025 | m72$HighErrorProb>0.025,]

m71<-mcList[[71]]
m71[m71$LowErrorProb>0.025 | m71$HighErrorProb>0.025,]




#all(errorVec < 0.025) # Hopefully evaluates to TRUE