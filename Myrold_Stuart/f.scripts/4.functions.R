# Table of Contents ----
# 1. Diameter of the sphere representing population distance from the multivariate optimum (Fisher's 'd')
# 2. Standardize the Measure of Mutation Size (Fisher's 'x') 
# 5. qqplots for test against exponential

# end table of contents

#==============================================================================#


# 1. Distance d from the multivariate optimum ----
# Following Orr, 1998, Evolution: The Population Genetics of Adaptation: The distribution of factors fixed during adaptive evolution
# assume the optimum sits at the origin of a multivariate set of coordinates. A population on an adaptive walk sits some distance from the origin, on the surface of an imaginary sphere.
# d is the diameter of the sphere. This value is used to calculate the standardized measure of a mutations effect size.





# 2. Standardize the Measure of Mutation Size ----
# Following Orr, 1998, Evolution: The Population Genetics of Adaptation: The distribution of factors fixed during adaptive evolution
#r is the raw step size in the original trait units. i.e., the euclidean distance in means from one step to the next in the adaptive walk.
#n is the number of traits considered in the multivariate optimum
#d is the diameter of the sphere that the adapting population is sitting on, which is 2 * the distance to the optimum

std.mutation.size <- function(r, n, dist.to.opt) {
  (r * sqrt(n)) / (2 * dist.to.opt)
}

# 3. C.I. calculation ----
c.i.calc <- function(std, n) {
  (1.960)(std/sqrt(n))
}

# 4. Sample Size ----
f.sample.size.noNAs <- function(x) length(x) - sum(is.na(x))

# 5. qqplot test ----
#y = r.mv
f.qqexp <-  function(y, mainlab, line=FALSE, step.type) { 
  y <- y[!is.na(y)]
  order.y <- y[order(y)]
  n <- length(y)
  x <- qexp(c(1:n)/(n+1))
  m <- mean(y)
  if (any(range(y)<0)) stop("Data contains negative values")
  ylim <- c(0, 3.3) #c(0,max(x)) YS Set this arbitrarily so that axes would be comparable across plots. 
  xlim <- c(0, 3.3)
  #par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))
  qqplot(x, y, xlab="Expected quantiles",xlim = xlim, ylim=ylim,ylab= "Sample Quantiles", cex = 0.6, pch = 16, las = 1, bty = "n") #use to be ylab=paste("Ordered step sizes (", step.type, ")", sep = "")
  mtext(side = 3, line = -2, text = mainlab, cex = 0.8)
  if (line) abline(0,m,lty=1)
  invisible()
  cor.qqplot <- cor(order.y, x)
  mtext(side = 3, line = -4, text = paste("qqplot correlation = ", round(cor.qqplot, 2), sep = ""), adj = 0.5, cex = 0.7)
  return(cor.qqplot)
}

# 6. qqexp ----
qqexp <-  function(y, line=FALSE, ...) { # Q-Qplot function found from: https://stats.stackexchange.com/questions/76994/how-do-i-check-if-my-data-fits-an-exponential-distribution
  y <- y[!is.na(y)]
  n <- length(y)
  x <- qexp(c(1:n)/(n+1))
  m <- mean(y)
  if (any(range(y)<0)) stop("Data contains negative values")
  ylim <- c(0,max(y))
  qqplot(x, y, xlab="Exponential plotting position",ylim=ylim,ylab="Ordered sample", ...)
  if (line) abline(0,m,lty=2)
  invisible()
}

#test.plot <- qqexp(x, line=TRUE) #resulting graph. seems like good fit?

#install.packages('KScorrect')
#library('KScorrect')
#ks.corrected <- LcKS(x,'pexp') #p=0.288
#hist(ks.corrected$D.sim)
