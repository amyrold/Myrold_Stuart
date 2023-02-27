# * 0.3 Data frame of r, x, and dist.to.opt at the start of each step during the adaptive walk. Using t1/2 method to define walk length ----
# r is the change in value from the start of the step to the end. Because the optimum end point is smaller in trait value than the lineage start, a positive value means a step toward the optimum. A negative value is a step away.
# distance is calculated at the start of a given step, so that the relative size of r can be calculated with respect to the distance from the optimum, giving x.
#make a new data frame that sets time 0 at the replacement event. Generate means, and 95% confidence intervals about them.
data.sbf.f <- as.data.frame(unique(data.l.sample.sbf.final[,6]) - 5000, row.names = as.character(unique(data.l.sample.sbf.final[,6])))
data.sbf.f <- na.omit(data.sbf.f) #remove last row

# find 95% C.I. for each trait and create column where year 5000 is set to 0
for(i in 1:3) {
  conf.int <- 1.96*sample.sd[,i]/sqrt(sample.N[,i])
  data.sbf.f <- cbind(data.sbf.f, sample.means[,i], conf.int)
};colnames(data.sbf.f) <- c('linII.time','ips','ips.CI','mds','mds.CI','ptt','ptt.CI')

#generate step size and distance-to-optimum statistics for downstream analysis, using the t1/2 method for length of adaptive walk
#ips
ips.r.x.d <- data.frame(matrix(NA, nrow = length(seq(ips.walk.start.end[1], ips.walk.start.end[2], by = 250)), ncol = 6))
colnames(ips.r.x.d) <- c("linII.time", "t.step.start", "dist.to.opt.step.start", "sample.mean.step.start" ,"r.step.size", "x.step.size")
ips.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(ips.walk.start.end[1], ips.walk.start.end[2], by = 250)] #time of step start
ips.r.x.d$t.step.start <- seq(ips.walk.start.end[1], ips.walk.start.end[2], by = 250) #time of step start duplicated as a check
ips.r.x.d$dist.to.opt.step.start <- data.sbf.f$ips[data.sbf.f$linII.time %in% ips.r.x.d$t.step.start] - optimum[1] #distance to optimum at the start of the step
ips.r.x.d$sample.mean.step.start <- data.sbf.f$ips[data.sbf.f$linII.time %in% ips.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:(nrow(ips.r.x.d)-1)) {
  ips.r.x.d$r.step.size[i] <- ips.r.x.d$sample.mean.step.start[i] - ips.r.x.d$sample.mean.step.start[i+1] #step made from the step start. positive implies step toward optimum.
}
ips.r.x.d$x.step.size <- std.mutation.size(r = ips.r.x.d$r.step.size, n = 1, dist.to.opt = ips.r.x.d$dist.to.opt.step.start)

#mds
mds.r.x.d <- data.frame(matrix(NA, nrow = length(seq(mds.walk.start.end[1], mds.walk.start.end[2], by = 250)), ncol = 6))
colnames(mds.r.x.d) <- c("linII.time", "t.step.start", "dist.to.opt.step.start", "sample.mean.step.start" ,"r.step.size", "x.step.size")
mds.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(mds.walk.start.end[1], mds.walk.start.end[2], by = 250)] #time of step start
mds.r.x.d$t.step.start <- seq(mds.walk.start.end[1], mds.walk.start.end[2], by = 250) #time of step start duplicated as a check
mds.r.x.d$dist.to.opt.step.start <- data.sbf.f$mds[data.sbf.f$linII.time %in% mds.r.x.d$t.step.start] - optimum[2] #distance to optimum at the start of the step
mds.r.x.d$sample.mean.step.start <- data.sbf.f$mds[data.sbf.f$linII.time %in% mds.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:(nrow(mds.r.x.d)-1)) {
  mds.r.x.d$r.step.size[i] <- mds.r.x.d$sample.mean.step.start[i] - mds.r.x.d$sample.mean.step.start[i+1] #step made from the step start. positive implies step toward optimum.
}
mds.r.x.d$x.step.size <- std.mutation.size(r = mds.r.x.d$r.step.size, n = 1, dist.to.opt = mds.r.x.d$dist.to.opt.step.start)

#ptt
ptt.r.x.d <- data.frame(matrix(NA, nrow = length(seq(ptt.walk.start.end[1], ptt.walk.start.end[2], by = 250)), ncol = 6))
colnames(ptt.r.x.d) <- c("linII.time", "t.step.start", "dist.to.opt.step.start", "sample.mean.step.start" ,"r.step.size", "x.step.size")
ptt.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(ptt.walk.start.end[1], ptt.walk.start.end[2], by = 250)] #time of step start
ptt.r.x.d$t.step.start <- seq(ptt.walk.start.end[1], ptt.walk.start.end[2], by = 250) #time of step start duplicated as a check
ptt.r.x.d$dist.to.opt.step.start <- data.sbf.f$ptt[data.sbf.f$linII.time %in% ptt.r.x.d$t.step.start] - optimum[3] #distance to optimum at the start of the step
ptt.r.x.d$sample.mean.step.start <- data.sbf.f$ptt[data.sbf.f$linII.time %in% ptt.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:(nrow(ptt.r.x.d)-1)) {
  ptt.r.x.d$r.step.size[i] <- ptt.r.x.d$sample.mean.step.start[i] - ptt.r.x.d$sample.mean.step.start[i+1] #step made from the step start. positive implies step toward optimum.
}
ptt.r.x.d$x.step.size <- std.mutation.size(r = ptt.r.x.d$r.step.size, n = 1, dist.to.opt = ptt.r.x.d$dist.to.opt.step.start)



# scrap -----
#fitting distributions  
#https://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf
x.norm<-rnorm(n=200,m=10,sd=2)  
hist(x.norm,main="Histogram of observed data")  
plot(density(x.norm),main="Density estimate of data")
plot(ecdf(x.norm),main= "Empirical cumulative distribution function")

z.norm<-(x.norm-mean(x.norm))/sd(x.norm) ## standardized data 
qqnorm(z.norm) ## drawing the QQplot 
abline(0,1) ## drawing a 45-degree reference line

x.wei<-rweibull(n=200,shape=2.1,scale=1.1)  ##  sampling  from  a  Weibull distribution with parameters shape=2.1 and scale=1.1 
x.teo<-rweibull(n=200,shape=2,  scale=1)  ##  theorical  quantiles  from  a Weibull population with known paramters shape=2 e scale=1 
qqplot(x.teo,x.wei,main="QQ-plot distr. Weibull") ## QQ-plot 
abline(0,1) ## a 45-degree reference line is plotted

fitdistr(x.gam,"gamma")
ks.test(x.wei,"pweibull", shape=2,scale=1) 

# QUESTION: do we need a logged optimum? ---- YS: Not sure yet (4 June 2021)
# once the optimum is within the 95% C.I., a '0' in the .adp column signifies end of adaptive walk
trait.cols <- c(2,4,6)
adp.names <- c('ptt.adp','mds.apd','ptt.adp')
for(i in 1:3){ # normal data
  data.sbf.f[which(data.sbf.f[,trait.cols[i]]-data.sbf.f[,trait.cols[i]+1] <= optimum[i]), adp.names[i]] <-0
  data.sbf.f[which(data.sbf.f[,trait.cols[i]]-data.sbf.f[,trait.cols[i]+1] > optimum[i]), adp.names[i]] <-1
}
for(i in 1:3){ # logged data
  trait.cols <- c(2,4,6)
  log.adp.names <- c('log.ptt.adp','log.mds.apd','log.ptt.adp')
  log.data.sbf.f[which(log.data.sbf.f[,trait.cols[i]]-log.data.sbf.f[,trait.cols[i]+1] <= optimum[i]), log.adp.names[i]] <-0
  log.data.sbf.f[which(log.data.sbf.f[,trait.cols[i]]-log.data.sbf.f[,trait.cols[i]+1] > optimum[i]), log.adp.names[i]] <-1
} 

# * 0.4 Create Cleaned Data Frame 
data.sbf.cleaned <- cbind(data.sbf.f[,c(1,2,3,8)],log.data.sbf.f[,c(2,3,8)],data.sbf.f[,c(4,5,9)],log.data.sbf.f[,c(4,5,9)],data.sbf.f[,c(6,7,10)],log.data.sbf.f[,c(6,7,10)])
write.csv(data.sbf.cleaned, file = paste(p.data.clean,'data.sbf.means.CI.csv', sep = ''))

# * 0.5 generate a new data frame of r, x, and dist.to.opt values for each step in the adaptive walk 
data.sbf.walk

#starting at the first arrival of lineage II, how does it evolve toward the optimum
movement <- vector()
movement.toward.optimum <- vector()
distance.from.optimum <- vector()
angle.from.optimum <- vector()

for(i in 1:(nrow(sample.means)-1)){
  sm.i <- sample.means[i, ] #the sample means at ti, where i starts at 5000 and proceeds to 22750. This is the starting point of each evolutionary trajectory.
  sm.iplus1 <- sample.means[i + 1, ] #this is the end point of each evolutionary trajectory.
  
  vector.timei.to.timeiplus1 <- sm.i - sm.iplus1 #this is the evolutionary vector from ti to ti+1
  vector.timei.to.opt <- sm.i - optimum #this is the vector toward the optimum at ti
  
  #angle between two vectors
  #theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  angle.between.movement.and.optimum <- angle.calc(x = vector.timei.to.timeiplus1, y = vector.timei.to.opt)
  
  #
  distance.in.direction.of.optimum <- cos(angle.between.movement.and.optimum) * sqrt(sum(vector.timei.to.timeiplus1^2))
  
  movement <- c(movement, round(sqrt(sum(vector.timei.to.timeiplus1^2)), 3))
  distance.from.optimum <- c(distance.from.optimum, round(sqrt(sum(vector.timei.to.opt^2)), 3))
  movement.toward.optimum <- c(movement.toward.optimum, round(distance.in.direction.of.optimum, 3))
  angle.from.optimum <- c(angle.from.optimum, round(angle.between.movement.and.optimum*(180/pi), 2))
}

#add vertical line when ps starts evolving

hist(movement.toward.optimum)
mean(movement.toward.optimum)
sd(movement.toward.optimum)
std.movement.toward.optimum <- (movement.toward.optimum - mean(movement.toward.optimum))/sd(movement.toward.optimum)
hist(std.movement.toward.optimum)
std.movement.toward.optimum

#data frame
movement.df <- data.frame(sort(unique(data.l.sample.sbf.final$interval)[-73])-5000, movement, movement.toward.optimum, distance.from.optimum, angle.from.optimum); colnames(movement.df) <- c("time.step.start", "movement.during.time.step", "movement.toward.optimum", "start.distance.from.optimum", "angle.of.movement.from.optimum")

#generate table to recreate ips.ci.r.x.d
boot.t <- data.frame(matrix(NA, nrow = length(unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))), ncol=7)) ; colnames(boot.t) <- c("linII.time",'boot.mean.start','boot.dist.opt','boot.step.size','is.real.step','boot.sderr','boot.ci.interval') ; boot.t$linII.time <- unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))

previous.bs <- NULL 
row = 0
#begin re-sampling data by time period
for(i in unique(seq(ips.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))){
  current.bs <- c() #used to store the 100 means of the current time period 'i'
  for(j in 2:101){ # re-sample and calculate means for time time period 'i'
    ips.boot.temp <- ips.boot.data$ips[ips.boot.data$linII.time == i] #generate observed data in first section
    ips.boot.temp.resample <- sample(x = ips.boot.temp, size = length(ips.boot.temp), replace = TRUE)
    current.bs <- append(current.bs, mean(ips.boot.temp.resample, na.rm = TRUE)) #store the mean from the current re-sampling
  }
  if(is.null(previous.bs) == FALSE){ # once we have 2 sets of means, begin calculating step sizes and confidence intervals
    boot.t.test <- t.test(x = unlist(previous.bs - current.bs)) #boot.t.test
    boot.ci <- c(unname(boot.t.test$conf.int)) # boot.ci
    boot.t[row,2] <- mean(previous.bs, na.rm = TRUE) # boot.mean
    boot.t[row,3] <- mean(previous.bs, na.rm = TRUE) - optimum[1] # distance from optimum 
    boot.t[row,4] <- mean(previous.bs - current.bs, na.rm = TRUE) #boot.step.size
    if(boot.ci[1]*boot.ci[2] < 0){boot.t[row,5] <- 0}else{boot.t[row,5] <- 1} #if it is a real step, assign value of 1
    boot.t[row,6] <- unname(boot.t.test$stderr) # boot.sderr
    boot.t[row,7] <- boot.ci[2] - unname(boot.t.test$estimate)
  }
  if( i == 9750){ #insert data for final row. there will be no step size data.
    boot.t[row+1,2] <- mean(current.bs, na.rm = TRUE) # boot.mean
    boot.t[row+1,3] <- mean(current.bs, na.rm = TRUE) - optimum[1] # distance from optimum
  }
  previous.bs <- current.bs # used to store 2 sets of means concurrently
  row = row+1
}



#AM:NB I believe this is some of our old tests. Up to you whether we dig this back up or stick to simply estimating the rate parameter?
norm.t <- rnorm(n = 200000, mean = 0, sd = 1)
fitdistr(x = norm.t, densfun = "normal")

exp.t <- rexp(n = 200, rate = 3)
hist(exp.t)
fitdistr(x = exp.t, densfun = "exponential")

# Not done
#fit exponential
mv.exp.rate <- fitdistr(x = na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0]) , densfun = "exponential")$estimate

ks.test(x = na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0]), y = "pexp", rate = mv.exp.rate) #D = 0.161, p = 0.514. Not different from exponential

#fit normal
mv.exp.mean.sd <- fitdistr(x = na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0]) , densfun = "normal")$estimate
ks.test(x = na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0]), y = "pnorm", mean = mv.exp.mean.sd[1], sd = mv.exp.mean.sd[2]) #D = 0.28, p = 0.037. Significantly different from normal at 0.05 level

#^^is there some circularity in fitting the distribution to get the parameters (fitdistr), then doing ks.test?^^ 

#Yup. Can't do that. https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best. See discussion of gamlss though. Promising.
x <- as.vector(na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0]))
hist(x, breaks = seq(0, 0.65, by = 0.025))
fit <- fitDist(x, k = 2, type = "realAll", trace = TRUE, try.gamlss = TRUE, )
summary(fit) #that was more complicated than I was hoping. maybe it's just a bunch of qqplots or something and we make the argument visually?
#picked inverse gaussian. Is that consistent with Kimura's intuition that the smallest effect size mutations will be lost to drift?


#YS: 10 June 2022. We haven't reported this anywhere, but let's keep the code. It would go back in analysis, since it's not plotting ----
# ** 3.5.1: use the function fitDist, which fits parametric gamlss.family distributions to a single data vector. The best distribution is picked using generalized AIC.
# *** 3.5.1.1 raw mv step size (r)
r.mv <- as.vector(na.omit(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size>0])) #multivariate step sizes toward the optimum
fit.r.mv <- fitDist(r.mv, k = 2, type = "realAll", trace = TRUE, try.gamlss = TRUE, )
summary(fit.r.mv) #SN2 (skew normal type 2 is best fit)
fit.r.mv$fits #gAIC for SN2 = -35.30. gAIC for Exponential = -27.26. dAIC = 8.04 units. 32 models better than exponential.

# *** 3.5.2.2 standardized mv step size (x)
x.mv <- as.vector(na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0])) #multivariate step sizes toward the optimum
fit.x.mv <- fitDist(x.mv, k = 2, type = "realAll", trace = TRUE, try.gamlss = TRUE, )
summary(fit.x.mv) #IG inverse gaussian is best fit
fit.x.mv$fits #gAIC for inverse gaussian = -48.50. gAIC for Exponential = -44.20. dAIC = 4.3 units. 19 models better than exponential.



# 4. using qqplot methods described in the comments below ----
#https://stats.stackexchange.com/questions/76994/how-do-i-check-if-my-data-fits-an-exponential-distribution
# * 4.1 raw mv step size (r) ----
r.mv <- as.vector(na.omit(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size>0])) #multivariate step sizes toward the optimum
r.mv.cor <- f.qqexp(r.mv, line = TRUE, step.type = "r") #correlation = 0.9875


# * 4.2 standardized mv step size (x) ----
x.mv <- as.vector(na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0])) #multivariate step sizes toward the optimum
x.mv.cor <- f.qqexp(x.mv, line = TRUE, step.type = "x") #correlation = 0.9830


# QUESTION: do we use the code below? ----
#YS response to question: I think it's used in slightly different form for the mv r and x histograms. Could delete this.
# * 4.1 Estimating rate parameters from the data, assuming and exponential
par(mfrow = c(2,1),mar = c(4,4,1,0), oma = c(1, 0, 0, 1))
# * 2.4 estimating rate parameter for mv step size (r) assuming exponential dist
eexp(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size > 0], ci = TRUE)
mv.r.rate <- unname(eexp(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size > 0])$parameters)
hist(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size > 0], xlab = "mv step size (r)", main = "", prob = TRUE)
curve(dexp(x, rate = mv.r.rate), col = 2, lty = 2, lwd = 2, add = TRUE)

# * 4.2 estimating rate parameter for mv standardized step size (x) assuming exponential dist 
eexp(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size > 0], ci = TRUE)
mv.r.rate <- unname(eexp(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size > 0])$parameters)
hist(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size > 0], xlab = "mv step size (r)", main = "", prob = TRUE)
curve(dexp(x, rate = mv.r.rate), col = 2, lty = 2, lwd = 2, add = TRUE)



# * 1.4 Length of adaptive walk in years - global variable ----
#NB: 5 June 2022. Not used. Chose to define end of adaptive walk by when trait mean CIs overlap with new optimum
#based on the t-1/2 estimates from Hunt et al. 2008
#t-1/2 was the amount of time (in generations) Hunt et al. 2008 estimated it took to get halfway to the optimum.
#assumed 2 years per generation
#adaptive walk time = t-1/2 * 2 * 2

#adaptive walk time dorsal spines
awt.mds <- 853 * 2 * 2

#adaptive walk time number of touching pterygiophores

awt.ptt <- 580 * 2 * 2

#adaptive walk time pelvic score
awt.ips <- 635 * 2 * 2

adaptive.walk.duration <- c(awt.ips, awt.mds, awt.ptt)

ips.walk.start.end <- c(3500, 3500 + adaptive.walk.duration[1]) #starts at linII.time-3500, because that is the last time mean(ips) was 3.0.
mds.walk.start.end <- c(0, 0 + adaptive.walk.duration[2])
ptt.walk.start.end <- c(0, 0 + adaptive.walk.duration[3])