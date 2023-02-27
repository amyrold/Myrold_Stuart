# Table of Contents ----
# 0. Read in Cleaned Data
      #0.1. Subset to only lineage II fish starting at year 5000 of temporal sequence L and going to end. #done
      #0.2. Create means by sample. raw and ln. add columns with 95% confidence intervals about the mean (single column with the interval size 1.96*s/sqrt(n)). add an adjusted year column where year 5000 is set to 0 and count up from there. choose t0-t1 or t1-t0 and we will be consistent from here out.
      #0.3. Define the adaptive walk temporal intervals for each of the three traits. 
            #e.g., for pelvic score, adaptive walk starts at adjusted year 3000. It will end the first time a the 95% confidence interval of a population sample contains the pelvic score optimum. Create an indicator variable column with 1s for being in adaptive walk and 0s for being outside the walk.
      #0.4. Generate data frames that contain r, x, and dist.to.opt values for each step in the adaptive walk, using t1/2 interval length
      #0.5. Generate data frames that contain r, x, and dist.to.opt values for each step in the adaptive walk, using confidence intervals to define adaptive walk

#Altered ToC:

# 0. Read Clean Data
    # 0.1 Subset to only Lineage II fish
    # 0.2 Create means by sample, mimicking bell et al

# 1. Analysis
    # 1.1 Univariate Analysis
    # 1.2 Multivariate Analysis

# 1. Bootstrapping
    # 1.1 Read Clean Data and re-sample
    # 1.2 ips bootstrapping and analysis
    # 1.3 mds bootstrapping and analysis
    # 1.4 ptt bootstrapping and analysis
    # 1.5 Multivariate Bootstrapping and Analysis
        # 1.5.1 Create variables and lists for use in MV analysis
        # 1.5.2 Bootstrapping and analysis
        # 1.5.3 t-tests

# 3. Trait correlations.  

# end Table of Contents
#==============================================================================#

# 0. Read Clean Data ----
data.l.sample.sbf #Data frame created in 2.data.manip.R

# * 0.1 Subset to only Lineage II fish ----
#Hunt et al. (2008) fit models to Lineage II fish only, starting at the replacement event.
tapply(X = data.l.sample.sbf$ips, INDEX = data.l.sample.sbf$interval, FUN = mean, na.rm = T) # pelvic score jumps at 5000 years, suggesting a mix of lineage I and lineage II fish; i.e., the one sample representing the replacement event.
data.l.sample.sbf$ips[data.l.sample.sbf$interval == 5000] #mix of 1.0s and 3.0s. Keep the 3.0s as lineage II fish.

#subsample to only the lineage II fish.
data.l.sample.sbf.final <- rbind(data.l.sample.sbf[which(data.l.sample.sbf$ips == 3 & data.l.sample.sbf$interval == 5000), ], data.l.sample.sbf[data.l.sample.sbf$interval > 5000, ])
data.l.sample.sbf.final <- data.l.sample.sbf.final[-which(data.l.sample.sbf.final$interval == 23250) , ] # removing a section without data

# * 0.2 Create means by sample, mimicking Bell et al. 2006 and Hunt et al. 2008 ----
#means, sum, and sd by sample. Useful for univariate and multivariate analysis.
sample.means <- as.data.frame(apply(data.l.sample.sbf.final[ , c(2, 3, 5)], MAR = 2, FUN = tapply, data.l.sample.sbf.final$interval, mean, na.rm = TRUE))
sample.N <- apply(X = data.l.sample.sbf.final[, 2:ncol(data.l.sample.sbf.final)], MARGIN = 2, FUN = tapply, data.l.sample.sbf.final$interval, f.sample.size.noNAs)
sample.sd <- as.data.frame(apply(data.l.sample.sbf.final[ , c(2, 3, 5)], MAR = 2, FUN = tapply, data.l.sample.sbf.final$interval, sd, na.rm = TRUE))

#make a new data frame that sets time 0 at the replacement event. Generate means, and 95% confidence intervals about them.
data.sbf.f <- as.data.frame(unique(data.l.sample.sbf.final[,6]) - 5000, row.names = as.character(unique(data.l.sample.sbf.final[,6])))
data.sbf.f <- na.omit(data.sbf.f) #remove last row

# find 95% C.I. for each trait and create column where year 5000 is set to 0
  for(i in 1:3) {
    conf.int <- 1.96*sample.sd[,i]/sqrt(sample.N[,i])
    data.sbf.f <- cbind(data.sbf.f, sample.means[,i], conf.int)
  };colnames(data.sbf.f) <- c('linII.time','ips','ips.CI','mds','mds.CI','ptt','ptt.CI')

# 1. Analysis ----
# * 1.1 Univariate Analysis ----
#Data frame of r, x, and dist.to.opt at the start of each step during the adaptive walk. Using confidence interval method to define walk length

#use confidence intervals about the mean and ask when mean estimates first overlap the optimum value from Hunt et al. 2008 (defined in 2.data.manip.R) to set adaptive walk interval

#ips
ips.lowerbound <- data.sbf.f$ips-data.sbf.f$ips.CI
data.sbf.f[which(ips.lowerbound < optimum[1]) , ] #sample at 9750 years is the first time it overlaps optimum. That year is the end of the walk, linII time. Walk starts at year 3500 linII time.
#mds
mds.lowerbound <- data.sbf.f$mds-data.sbf.f$mds.CI
data.sbf.f[which(mds.lowerbound < optimum[2]) , ] #sample at 3000 years is the first time it overlaps optimum. That year is the end of the walk. Walk starts at year 0 linIItime
#ptt
ptt.lowerbound <- data.sbf.f$ptt-data.sbf.f$ptt.CI
data.sbf.f[which(ptt.lowerbound < optimum[3]) , ] #sample at 2250 years is the first time it overlaps optimum. That year is the end of the walk. Walk starts at year 0 linIItime.

#adaptive walks, based on confidence intervals.
ips.ci.walk.start.end <- c(3500, 9750) #15 intervals longer than Hunt et al. t1/2 method. Did they make a mistake? Something to do with logged values?
mds.ci.walk.start.end <- c(0, 3000) #one interval shorter relative to Hunt et al. t1/2 method.
ptt.ci.walk.start.end <- c(0, 2250) #doesn't change from Hunt et al. t1/2 method.

#ips.ci
ips.ci.r.x.d <- data.frame(matrix(NA, nrow = length(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250)), ncol = 6))
colnames(ips.ci.r.x.d) <- c("linII.time", "t.step.start", "dist.to.opt.step.start", "sample.mean.step.start" ,"r.step.size", "x.step.size")
ips.ci.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250)] #time of step start
ips.ci.r.x.d$t.step.start <- seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250) #time of step start duplicated as a check
ips.ci.r.x.d$dist.to.opt.step.start <- data.sbf.f$ips[data.sbf.f$linII.time %in% ips.ci.r.x.d$t.step.start] - optimum[1] #distance to optimum at the start of the step
ips.ci.r.x.d$sample.mean.step.start <- data.sbf.f$ips[data.sbf.f$linII.time %in% ips.ci.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:(nrow(ips.ci.r.x.d)-1)) {
  ips.ci.r.x.d$r.step.size[i] <- ips.ci.r.x.d$sample.mean.step.start[i] - ips.ci.r.x.d$sample.mean.step.start[i+1] #step made from the step start. positive implies step toward optimum.
}
ips.ci.r.x.d$x.step.size <- std.mutation.size(r = ips.ci.r.x.d$r.step.size, n = 1, dist.to.opt = ips.ci.r.x.d$dist.to.opt.step.start)

#mds.ci
mds.ci.r.x.d <- data.frame(matrix(NA, nrow = length(seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250)), ncol = 6))
colnames(mds.ci.r.x.d) <- c("linII.time", "t.step.start", "dist.to.opt.step.start", "sample.mean.step.start" ,"r.step.size", "x.step.size")
mds.ci.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250)] #time of step start
mds.ci.r.x.d$t.step.start <- seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250) #time of step start duplicated as a check
mds.ci.r.x.d$dist.to.opt.step.start <- data.sbf.f$mds[data.sbf.f$linII.time %in% mds.ci.r.x.d$t.step.start] - optimum[2] #distance to optimum at the start of the step
mds.ci.r.x.d$sample.mean.step.start <- data.sbf.f$mds[data.sbf.f$linII.time %in% mds.ci.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:(nrow(mds.ci.r.x.d)-1)) {
  mds.ci.r.x.d$r.step.size[i] <- mds.ci.r.x.d$sample.mean.step.start[i] - mds.ci.r.x.d$sample.mean.step.start[i+1] #step made from the step start. positive implies step toward optimum.
}
mds.ci.r.x.d$x.step.size <- std.mutation.size(r = mds.ci.r.x.d$r.step.size, n = 1, dist.to.opt = mds.ci.r.x.d$dist.to.opt.step.start)

#ptt.ci
ptt.ci.r.x.d <- data.frame(matrix(NA, nrow = length(seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250)), ncol = 6))
colnames(ptt.ci.r.x.d) <- c("linII.time", "t.step.start", "dist.to.opt.step.start", "sample.mean.step.start" ,"r.step.size", "x.step.size")
ptt.ci.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250)] #time of step start
ptt.ci.r.x.d$t.step.start <- seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250) #time of step start duplicated as a check
ptt.ci.r.x.d$dist.to.opt.step.start <- data.sbf.f$ptt[data.sbf.f$linII.time %in% ptt.ci.r.x.d$t.step.start] - optimum[3] #distance to optimum at the start of the step
ptt.ci.r.x.d$sample.mean.step.start <- data.sbf.f$ptt[data.sbf.f$linII.time %in% ptt.ci.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:(nrow(ptt.ci.r.x.d)-1)) {
  ptt.ci.r.x.d$r.step.size[i] <- ptt.ci.r.x.d$sample.mean.step.start[i] - ptt.ci.r.x.d$sample.mean.step.start[i+1] #step made from the step start. positive implies step toward optimum.
}
ptt.ci.r.x.d$x.step.size <- std.mutation.size(r = ptt.ci.r.x.d$r.step.size, n = 1, dist.to.opt = ptt.ci.r.x.d$dist.to.opt.step.start)

#==============================================================================#

# * 1.2 Multivariate Analysis ----
# make r,x,d table in multivariate space

#xlabel = "Relative time of deposition (in years since Lineage II appeared)"
xlabel = "Relative time of deposition in years (since start of adaptive walk)"

# Adaptive walk interval
#adaptive walk starts at time 0 and ends at 9750 when ips reaches the optimum
#9750 based on confidence interval approach from above
mv.adaptive.walk <- c(0, 9750)

mv.ci.r.x.d <- data.frame(matrix(NA, nrow = length(seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250)), ncol = 8))
colnames(mv.ci.r.x.d) <- c("linII.time", "t.step.start", "ips.sample.mean.step.start", "mds.sample.mean.step.start", "ptt.sample.mean.step.start", "euclid.dist.to.opt.step.start","mv.r.step.size", "mv.x.step.size")
mv.ci.r.x.d$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250)]
mv.ci.r.x.d$t.step.start <- seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250) #time of step start duplicated as a check
mv.ci.r.x.d$ips.sample.mean.step.start <- data.sbf.f$ips[data.sbf.f$linII.time %in% mv.ci.r.x.d$t.step.start] #mean of the sample at start of step.
mv.ci.r.x.d$mds.sample.mean.step.start <- data.sbf.f$mds[data.sbf.f$linII.time %in% mv.ci.r.x.d$t.step.start] #mean of the sample at start of step.
mv.ci.r.x.d$ptt.sample.mean.step.start <- data.sbf.f$ptt[data.sbf.f$linII.time %in% mv.ci.r.x.d$t.step.start] #mean of the sample at start of step.
for(i in 1:length(mv.ci.r.x.d$euclid.dist.to.opt.step.start)){
  x.t <- mv.ci.r.x.d[i , 3:5] #sample means
  y.t <- optimum
  mv.ci.r.x.d$euclid.dist.to.opt.step.start[i] <- dist(rbind(x.t, y.t), method = "euclidean")
}
for(i in 1:(length(mv.ci.r.x.d$mv.r.step.size)-1)){
  mv.ci.r.x.d$mv.r.step.size[i] <- mv.ci.r.x.d$euclid.dist.to.opt.step.start[i] - mv.ci.r.x.d$euclid.dist.to.opt.step.start[i + 1] #step made from step start. positive implies step toward optimum
}
mv.ci.r.x.d$mv.x.step.size <- std.mutation.size(r = mv.ci.r.x.d$mv.r.step.size, n = 3, dist.to.opt = mv.ci.r.x.d$euclid.dist.to.opt.step.start)
#end section 1)

#==============================================================================#

# 2. Bootstrapping ----
#Goal: to estimate sampling error on the step sizes. 
#To do. 
# 1) For a given trait, within each temporal sample, re-sample the data with replacement
# 2) Then, with that new data set, calculate the step size of each step.
# 3) Repeat 100 times. 
# 4) Get a mean and standard deviation across bootstrap replicates for each step.
# 5) Plot the mean and sd on the plots generated in the data-visualization sections.
# 6) use a one-sample t-test to ask which step size distributions are significantly different from 0. These could be considered real steps, distinguished from sampling error. Those real steps are what are analyzed in 1.1.4, 1.2, 2.3.3, and 2.4

# * 2.1 Read Clean Data and Re-sample ----
bootstrap.data <- data.l.sample.sbf.final
bootstrap.data$linII.time <- bootstrap.data$interval - 5000

ips.boot.data <- bootstrap.data[bootstrap.data$linII.time %in% seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250) , c(2,7)] #data just from adaptive walk

mds.boot.data <- bootstrap.data[bootstrap.data$linII.time %in% seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250) , c(3,7)] #data just from adaptive walk

ptt.boot.data <- bootstrap.data[bootstrap.data$linII.time %in% seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250) , c(5,7)] #data just from adaptive walk

#multivariate adaptive walk starts at time 0 and ends at 9750 when ips reaches the optimum
#9750 based on confidence interval approach from above
mv.adaptive.walk <- c(0, 9750)
mv.boot.data <- bootstrap.data[bootstrap.data$linII.time %in% seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250) , c(2,3,5,7)] #data just from adaptive walk

#==============================================================================#

# * 2.2 ips Bootstrapping and Analysis ----
#generate table to gather means from bootstrap sampling
bootstrapped.ips.data.t <- data.frame(matrix(NA, nrow = length(unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))), ncol = 101)) ; colnames(bootstrapped.ips.data.t) <- c("linII.time", paste("bs.means.", seq(1:100), sep = "")) ; bootstrapped.ips.data.t$linII.time <- unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))

#add to table above
for(j in 2:101){
  #resample ips by section
  for(i in unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))){ #iterate through sections in adaptive walk
    ips.boot.temp <- ips.boot.data$ips[ips.boot.data$linII.time == i] #generate observed data in first section
    ips.boot.temp.resample <- sample(x = ips.boot.temp, size = length(ips.boot.temp), replace = TRUE)
    bootstrapped.ips.data.t[bootstrapped.ips.data.t$linII.time == i, j] <- mean(ips.boot.temp.resample, na.rm = TRUE)
  }
}

#generate table to recreate ips.ci.r.x.d
ips.boot.t <- data.frame(matrix(NA, nrow = length(unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))), ncol=12)) ; colnames(ips.boot.t) <- c("linII.time",'sample.mean.step.start.bs','sample.mean.step.start.sd.bs','dist.to.opt.step.start.bs','r.step.size.bs','r.step.size.ci.bs','r.step.size.t-stat.bs', 'r.step.size.p-value.bs','x.step.size.bs','x.step.size.ci.bs', 'x.step.size.t-stat.bs','x.step.size.p-value.bs') ; ips.boot.t$linII.time <- unique(seq(ips.ci.walk.start.end[1], ips.ci.walk.start.end[2], by = 250))
for(i in 1:length(ips.boot.t$linII.time)){
  ips.previous <- unlist(bootstrapped.ips.data.t[i,-1])
  ips.current <- unlist(bootstrapped.ips.data.t[i+1,-1]) # save data to calculate step size.
  
  if(i == 26){ # This will help terminate before errors are thrown in the last row
    ips.boot.t[i,2] <- mean(ips.previous, na.rm = TRUE) # boot.mean
    ips.boot.t[i,4] <- mean(ips.previous, na.rm = TRUE) - optimum[1] # distance from optimum 
    break;
  }
  ips.boot.t[i,2] <- mean(ips.previous, na.rm = TRUE) # sample.mean.step.start.bs
  ips.boot.t[i,3] <- sd(ips.previous, na.rm = TRUE) # sample.mean.step.start.sd.bs
  ips.boot.t[i,4] <- mean(ips.previous, na.rm = TRUE) - optimum[1] # dist.to.opt.step 
  ips.boot.t[i,5] <- mean(ips.previous - ips.current, na.rm = TRUE) # mean r.step.size.bs
  r.step.size.bs.t.test <- t.test(x = ips.previous - ips.current, mu=0) #t-test to see if raw step size distribution includes 0.
  r.step.size.bs.ci <- c(unname(r.step.size.bs.t.test$conf.int)) # boot.ci
  ips.boot.t[i,6] <- r.step.size.bs.ci[2] - unname(r.step.size.bs.t.test$estimate) # r.step.size.ci.bs. This is half of the CI.
  ips.boot.t[i,7] <- round(r.step.size.bs.t.test$statistic, 2) # r.step.size.t-stat.bs #t-test to see if raw step size distribution includes 0.
  ips.boot.t[i,8] <- round(r.step.size.bs.t.test$p.value, 3) # r.step.size.p-value.bs #pvalue to see if raw step size distribution includes 0.
  ips.boot.t[i,9] <- std.mutation.size(ips.boot.t[i,5], 1, ips.boot.t[i,4]) # x.step.size.bs 
  x.step.size.bs.t.test <- t.test(x = std.mutation.size((ips.previous - ips.current), 1, ips.previous), mu=0) #t-test to see if standardized step size distribution includes 0.
  x.step.size.bs.ci <- c(unname(x.step.size.bs.t.test$conf.int)) # boot.ci
  ips.boot.t[i,10] <- x.step.size.bs.ci[2] - unname(x.step.size.bs.t.test$estimate) # x.step.size.ci.bs
  ips.boot.t[i,11] <- round(unname(x.step.size.bs.t.test$statistic), 2) # x.step.size.t-stat.bs #t-test to see if standardized step size distribution includes 0.
  ips.boot.t[i,12] <- round(unname(x.step.size.bs.t.test$p.value), 3) # x.step.size.p-value #pvalue to see if standardized step size distribution includes 0.
  #if(ips.boot.ci[1]*ips.boot.ci[2] < 0){ips.boot.t[i,6] <- 0}else{ips.boot.t[i,6] <- 1} #if it is a real step, assign value of 1
}

#==============================================================================#

# * 2.3 mds Bootstrapping and Analysis ----
#generate table to gather means from bootstrap sampling
bootstrapped.mds.data.t <- data.frame(matrix(NA, nrow = length(unique(seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250))), ncol = 101)) ; colnames(bootstrapped.mds.data.t) <- c("linII.time", paste("bs.means.", seq(1:100), sep = "")) ; bootstrapped.mds.data.t$linII.time <- unique(seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250))

#add to table above
for(j in 2:101){
  #resample mds by section
  for(i in unique(seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250))){ #iterate through sections in adaptive walk
    mds.boot.temp <- mds.boot.data$mds[mds.boot.data$linII.time == i] #generate observed data in first section
    mds.boot.temp.resample <- sample(x = mds.boot.temp, size = length(mds.boot.temp), replace = TRUE)
    bootstrapped.mds.data.t[bootstrapped.mds.data.t$linII.time == i, j] <- mean(mds.boot.temp.resample, na.rm = TRUE)
  }
}

#mds.r.x.d$x.step.size <- std.mutation.size(r = mds.r.x.d$r.step.size, n = 1, dist.to.opt = mds.r.x.d$dist.to.opt.step.start)


#generate table to recreate mds.ci.r.x.d
mds.boot.t <- data.frame(matrix(NA, nrow = length(unique(seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250))), ncol=12)) ; colnames(mds.boot.t) <- c("linII.time",'sample.mean.step.start.bs','sample.mean.step.start.sd.bs','dist.to.opt.step.start.bs','r.step.size.bs','r.step.size.ci.bs','r.step.size.t-stat.bs', 'r.step.size.p-value.bs','x.step.size.bs','x.step.size.ci.bs', 'x.step.size.t-stat.bs','x.step.size.p-value.bs') ; mds.boot.t$linII.time <- unique(seq(mds.ci.walk.start.end[1], mds.ci.walk.start.end[2], by = 250))
for(i in 1:length(mds.boot.t$linII.time)){
  mds.previous <- unlist(bootstrapped.mds.data.t[i,-1])
  mds.current <- unlist(bootstrapped.mds.data.t[i+1,-1]) # save data to calculate step size.
  
  if(i == length(mds.boot.t$linII.time)){ # This will help terminate before errors are thrown in the last row
    mds.boot.t[i,2] <- mean(mds.previous, na.rm = TRUE) # boot.mean
    mds.boot.t[i,4] <- mean(mds.previous, na.rm = TRUE) - optimum[1] # distance from optimum 
    break;
  }
  mds.boot.t[i,2] <- mean(mds.previous, na.rm = TRUE) # sample.mean.step.start.bs
  mds.boot.t[i,3] <- sd(mds.previous, na.rm = TRUE) # sample.mean.step.start.sd.bs
  mds.boot.t[i,4] <- mean(mds.previous, na.rm = TRUE) - optimum[1] # dist.to.opt.step 
  mds.boot.t[i,5] <- mean(mds.previous - mds.current, na.rm = TRUE) # r.step.size.bs
  r.step.size.bs.t.test <- t.test(x = mds.previous - mds.current, mu=0)
  r.step.size.bs.ci <- c(unname(r.step.size.bs.t.test$conf.int)) # boot.ci
  mds.boot.t[i,6] <- r.step.size.bs.ci[2] - unname(r.step.size.bs.t.test$estimate) # r.step.size.ci.bs
  mds.boot.t[i,7] <- round(r.step.size.bs.t.test$statistic, 2) # r.step.size.t-stat.bs
  mds.boot.t[i,8] <- round(r.step.size.bs.t.test$p.value, 3) # r.step.size.p-value.bs 
  mds.boot.t[i,9] <- std.mutation.size(mds.boot.t[i,5], 1, mds.boot.t[i,4]) # x.step.size.bs 
  x.step.size.bs.t.test <- t.test(x = std.mutation.size((mds.previous - mds.current), 1, mds.previous), mu=0)
  x.step.size.bs.ci <- c(unname(x.step.size.bs.t.test$conf.int)) # boot.ci
  mds.boot.t[i,10] <- x.step.size.bs.ci[2] - unname(x.step.size.bs.t.test$estimate) # x.step.size.ci.bs
  mds.boot.t[i,11] <- round(unname(x.step.size.bs.t.test$statistic), 2) # x.step.size.t-stat.bs
  mds.boot.t[i,12] <- round(unname(x.step.size.bs.t.test$p.value), 3) # x.step.size.p-value
  #if(mds.boot.ci[1]*mds.boot.ci[2] < 0){mds.boot.t[i,6] <- 0}else{mds.boot.t[i,6] <- 1} #if it is a real step, assign value of 1
}

#==============================================================================#

# * 2.4 ptt Bootstrapping and Analysis ----
#generate table to gather means from bootstrap sampling
bootstrapped.ptt.data.t <- data.frame(matrix(NA, nrow = length(unique(seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250))), ncol = 101)) ; colnames(bootstrapped.ptt.data.t) <- c("linII.time", paste("bs.means.", seq(1:100), sep = "")) ; bootstrapped.ptt.data.t$linII.time <- unique(seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250))

#add to table above
for(j in 2:101){
  #resample ptt by section
  for(i in unique(seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250))){ #iterate through sections in adaptive walk
    ptt.boot.temp <- ptt.boot.data$ptt[ptt.boot.data$linII.time == i] #generate observed data in first section
    ptt.boot.temp.resample <- sample(x = ptt.boot.temp, size = length(ptt.boot.temp), replace = TRUE)
    bootstrapped.ptt.data.t[bootstrapped.ptt.data.t$linII.time == i, j] <- mean(ptt.boot.temp.resample, na.rm = TRUE)
  }
}

#ptt.r.x.d$x.step.size <- std.mutation.size(r = ptt.r.x.d$r.step.size, n = 1, dist.to.opt = ptt.r.x.d$dist.to.opt.step.start)

#generate table to recreate ptt.ci.r.x.d
ptt.boot.t <- data.frame(matrix(NA, nrow = length(unique(seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250))), ncol=12)) ; colnames(ptt.boot.t) <- c("linII.time",'sample.mean.step.start.bs','sample.mean.step.start.sd.bs','dist.to.opt.step.start.bs','r.step.size.bs','r.step.size.ci.bs','r.step.size.t-stat.bs', 'r.step.size.p-value.bs','x.step.size.bs','x.step.size.ci.bs', 'x.step.size.t-stat.bs','x.step.size.p-value.bs') ; ptt.boot.t$linII.time <- unique(seq(ptt.ci.walk.start.end[1], ptt.ci.walk.start.end[2], by = 250))
for(i in 1:length(ptt.boot.t$linII.time)){
  ptt.previous <- unlist(bootstrapped.ptt.data.t[i,-1])
  ptt.current <- unlist(bootstrapped.ptt.data.t[i+1,-1]) # save data to calculate step size.
  
  if(i == length(ptt.boot.t$linII.time)){ # This will help terminate before errors are thrown in the last row
    ptt.boot.t[i,2] <- mean(ptt.previous, na.rm = TRUE) # boot.mean
    ptt.boot.t[i,4] <- mean(ptt.previous, na.rm = TRUE) - optimum[1] # distance from optimum 
    break;
  }
  ptt.boot.t[i,2] <- mean(ptt.previous, na.rm = TRUE) # sample.mean.step.start.bs
  ptt.boot.t[i,3] <- sd(ptt.previous, na.rm = TRUE) # sample.mean.step.start.sd.bs
  ptt.boot.t[i,4] <- mean(ptt.previous, na.rm = TRUE) - optimum[1] # dist.to.opt.step 
  ptt.boot.t[i,5] <- mean(ptt.previous - ptt.current, na.rm = TRUE) # r.step.size.bs
  r.step.size.bs.t.test <- t.test(x = ptt.previous - ptt.current, mu=0)
  r.step.size.bs.ci <- c(unname(r.step.size.bs.t.test$conf.int)) # boot.ci
  ptt.boot.t[i,6] <- r.step.size.bs.ci[2] - unname(r.step.size.bs.t.test$estimate) # r.step.size.ci.bs
  ptt.boot.t[i,7] <- round(r.step.size.bs.t.test$statistic, 2) # r.step.size.t-stat.bs
  ptt.boot.t[i,8] <- round(r.step.size.bs.t.test$p.value, 3) # r.step.size.p-value.bs 
  ptt.boot.t[i,9] <- std.mutation.size(ptt.boot.t[i,5], 1, ptt.boot.t[i,4]) # x.step.size.bs 
  x.step.size.bs.t.test <- t.test(x = std.mutation.size((ptt.previous - ptt.current), 1, ptt.previous), mu=0)
  x.step.size.bs.ci <- c(unname(x.step.size.bs.t.test$conf.int)) # boot.ci
  ptt.boot.t[i,10] <- x.step.size.bs.ci[2] - unname(x.step.size.bs.t.test$estimate) # x.step.size.ci.bs
  ptt.boot.t[i,11] <- round(unname(x.step.size.bs.t.test$statistic), 2) # x.step.size.t-stat.bs
  ptt.boot.t[i,12] <- round(unname(x.step.size.bs.t.test$p.value), 3) # x.step.size.p-value
  #if(ptt.boot.ci[1]*ptt.boot.ci[2] < 0){ptt.boot.t[i,6] <- 0}else{ptt.boot.t[i,6] <- 1} #if it is a real step, assign value of 1
}

#==============================================================================# 

# * 2.5 multivariate Bootstrapping and Analysis ----
#Goal: bootstrap what we did in section 2.2. 
#re-sample the three traits at each sample, without mixing
#get one set of means, then statistics, then one set of means, then statistics, etc. which is different than the strategy above

# ** 2.5.1 Create variables and lists for use in mv analysis ----
mv.years <- seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250) #store the sequence of years in the adaptive walk
# create a list that stores all the fish as individuals
mv.boot.ls <- vector(mode = 'list', length = length(mv.years)) 
for(i in 1:length(mv.years)){
  mv.boot.t <- mv.boot.data[mv.boot.data$linII.time == mv.years[i],]
  mv.boot.ls.t <- vector(mode = 'list', length = length(mv.boot.t[,1]))
  for(j in 1:length(mv.boot.t[,1])){
    mv.boot.ls.t[[j]] <- c(unname(unlist(mv.boot.t[j,])))
  }
  mv.boot.ls[[i]] <- mv.boot.ls.t
}
#Create the list we will use to store the final bootstrapped data
mv.ci.r.x.d.bs <- data.frame(matrix(NA, nrow = length(seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250)), ncol = 8))
colnames(mv.ci.r.x.d.bs) <- c("linII.time", "t.step.start.bs", "ips.sample.mean.step.start.bs", "mds.sample.mean.step.start.bs", "ptt.sample.mean.step.start.bs", "euclid.dist.to.opt.step.start.bs","mv.r.step.size.bs", "mv.x.step.size.bs")
mv.ci.r.x.d.bs$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250)]
mv.ci.r.x.d.bs$t.step.start.bs <- seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250) #time of step start duplicated as a check

# ** 2.5.2 bootstrapping and analysis ----
for(i in 1:100){ # determine number of bootstraps to perform
  for(j in 1:length(mv.boot.ls)){ #re-sample data by time-series and append to a new data frame
    current.ls <- mv.boot.ls[[j]]
    mv.bs.t <- sample(x = current.ls, size = length(current.ls), replace = TRUE)
    if(j==1){
      mv.bs.df.t <- data.frame(matrix(unlist(mv.bs.t), nrow=length(mv.bs.t), byrow=TRUE))
    } else{
      mv.bs.df.t <- rbind(mv.bs.df.t, data.frame(matrix(unlist(mv.bs.t), nrow=length(mv.bs.t), byrow=TRUE)))
    }
  }
  #calculate the means by time-series
  mv.means.bs <- as.data.frame(apply(mv.bs.df.t[ , c(1, 2, 3, 4)], MAR = 2, FUN = tapply, mv.bs.df.t$X4, mean, na.rm = TRUE))
  if(i == 1){ # for our first bootstrap, we store our analyses into the final dataframe
    mv.ci.r.x.d.bs$ips.sample.mean.step.start.bs <- mv.means.bs$X1[mv.means.bs$X4 %in% mv.ci.r.x.d.bs$t.step.start.bs] #mean of the sample at start of step.
    mv.ci.r.x.d.bs$mds.sample.mean.step.start.bs <- mv.means.bs$X2[mv.means.bs$X4 %in% mv.ci.r.x.d.bs$t.step.start.bs] #mean of the sample at start of step.
    mv.ci.r.x.d.bs$ptt.sample.mean.step.start.bs <- mv.means.bs$X3[mv.means.bs$X4 %in% mv.ci.r.x.d.bs$t.step.start.bs] #mean of the sample at start of step.
    for(m in 1:length(mv.ci.r.x.d.bs$euclid.dist.to.opt.step.start.bs)){
      x.t <- mv.ci.r.x.d.bs[m , 3:5] #sample means
      y.t <- optimum
      mv.ci.r.x.d.bs$euclid.dist.to.opt.step.start.bs[m] <- dist(rbind(x.t, y.t), method = "euclidean")
    }
    for(m in 1:(length(mv.ci.r.x.d.bs$mv.r.step.size.bs)-1)){
      mv.ci.r.x.d.bs$mv.r.step.size.bs[m] <- mv.ci.r.x.d.bs$euclid.dist.to.opt.step.start.bs[m] - mv.ci.r.x.d.bs$euclid.dist.to.opt.step.start.bs[m + 1] #step made from step start. positive implies step toward optimum
    }
    mv.ci.r.x.d.bs$mv.x.step.size.bs <- std.mutation.size(r = mv.ci.r.x.d.bs$mv.r.step.size.bs, n = 3, dist.to.opt = mv.ci.r.x.d.bs$euclid.dist.to.opt.step.start.bs)
  } 
  else{ #for remaining bootstraps, we create a temporary data frame and store the mean between the temp and final into the final.
    mv.ci.r.x.d.bs.t <- data.frame(matrix(NA, nrow = length(seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250)), ncol = 8))
    colnames(mv.ci.r.x.d.bs.t) <- c("linII.time", "t.step.start.bs", "ips.sample.mean.step.start.bs", "mds.sample.mean.step.start.bs", "ptt.sample.mean.step.start.bs", "euclid.dist.to.opt.step.start.bs","mv.r.step.size.bs", "mv.x.step.size.bs")
    mv.ci.r.x.d.bs.t$linII.time <- data.sbf.f$linII.time[data.sbf.f$linII.time %in% seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250)]
    mv.ci.r.x.d.bs.t$t.step.start.bs <- seq(mv.adaptive.walk[1], mv.adaptive.walk[2], by = 250) #time of step start duplicated as a check
    
    mv.ci.r.x.d.bs.t$ips.sample.mean.step.start.bs <- mv.means.bs$X1[mv.means.bs$X4 %in% mv.ci.r.x.d.bs.t$t.step.start.bs] #mean of the sample at start of step.
    mv.ci.r.x.d.bs.t$mds.sample.mean.step.start.bs <- mv.means.bs$X2[mv.means.bs$X4 %in% mv.ci.r.x.d.bs.t$t.step.start.bs] #mean of the sample at start of step.
    mv.ci.r.x.d.bs.t$ptt.sample.mean.step.start.bs <- mv.means.bs$X3[mv.means.bs$X4 %in% mv.ci.r.x.d.bs.t$t.step.start.bs] #mean of the sample at start of step.
    for(m in 1:length(mv.ci.r.x.d.bs.t$euclid.dist.to.opt.step.start.bs)){
      x.t <- mv.ci.r.x.d.bs.t[m , 3:5] #sample means
      y.t <- optimum
      mv.ci.r.x.d.bs.t$euclid.dist.to.opt.step.start.bs[m] <- dist(rbind(x.t, y.t), method = "euclidean")
    }
    for(m in 1:(length(mv.ci.r.x.d.bs.t$mv.r.step.size.bs)-1)){
      mv.ci.r.x.d.bs.t$mv.r.step.size.bs[m] <- mv.ci.r.x.d.bs.t$euclid.dist.to.opt.step.start.bs[m] - mv.ci.r.x.d.bs.t$euclid.dist.to.opt.step.start.bs[m + 1] #step made from step start. positive implies step toward optimum
    }
    mv.ci.r.x.d.bs.t$mv.x.step.size.bs <- std.mutation.size(r = mv.ci.r.x.d.bs.t$mv.r.step.size.bs, n = 3, dist.to.opt = mv.ci.r.x.d.bs.t$euclid.dist.to.opt.step.start.bs)
    
    #we use the abind function to create a 3D array of data frames (mv.ci.r.x.d.bs)
    mv.ci.r.x.d.bs = abind(mv.ci.r.x.d.bs, mv.ci.r.x.d.bs.t, along = 3)
  }
}

# ** 2.5.3 t-tests ----
mv.ci.r.x.d.bs.f <- as.data.frame(rowMeans(mv.ci.r.x.d.bs, dims = 2))
for(i in 1:(length(mv.years)-1)){ #use
  if(i == 1){
    r.step.size.tstat.bs <- c()
    r.step.size.pvalue.bs <- c()
    r.step.size.ci.bs <- c()
    r.step.size.sd.bs <- c()
    x.step.size.ci.bs <- c()
    x.step.size.tstat.bs <- c()
    x.step.size.pvalue.bs <- c()
    x.step.size.sd.bs <- c()
  }
  r.t.test <- t.test(mv.ci.r.x.d.bs[i,7,], mu=0)
  x.t.test <- t.test(mv.ci.r.x.d.bs[i,8,], mu=0)
  r.ci <- r.t.test$conf.int
  x.ci <- x.t.test$conf.int
  r.step.size.sd.bs <- append(r.step.size.sd.bs, sd(mv.ci.r.x.d.bs[i,7,]))
  r.step.size.ci.bs <- append(r.step.size.ci.bs, (r.ci[2] - unname(r.t.test$estimate)) )
  r.step.size.tstat.bs <- append(r.step.size.tstat.bs, unname(r.t.test$statistic))
  r.step.size.pvalue.bs <- append(r.step.size.pvalue.bs, r.t.test$p.value)
  x.step.size.sd.bs <- append(x.step.size.sd.bs, sd(mv.ci.r.x.d.bs[i,8,]))
  x.step.size.ci.bs <- append(x.step.size.ci.bs, (x.ci[2] - unname(x.t.test$estimate)) )
  x.step.size.tstat.bs <- append(x.step.size.tstat.bs, unname(x.t.test$statistic))
  x.step.size.pvalue.bs <- append(x.step.size.pvalue.bs, x.t.test$p.value)
  if(i == (length(mv.years)-1)){
    mv.ci.r.x.d.bs.f$r.step.size.sd.bs  <- c(r.step.size.sd.bs, NA)
    mv.ci.r.x.d.bs.f$r.step.size.ci.bs  <- c(r.step.size.ci.bs, NA)
    mv.ci.r.x.d.bs.f$r.step.size.tstat.bs  <- c(r.step.size.tstat.bs, NA)
    mv.ci.r.x.d.bs.f$r.step.size.pvalue.bs  <- c(r.step.size.pvalue.bs, NA)
    mv.ci.r.x.d.bs.f$x.step.size.sd.bs  <- c(x.step.size.sd.bs, NA)
    mv.ci.r.x.d.bs.f$x.step.size.ci.bs  <- c(x.step.size.ci.bs, NA)
    mv.ci.r.x.d.bs.f$x.step.size.tstat.bs  <- c(x.step.size.tstat.bs, NA)
    mv.ci.r.x.d.bs.f$x.step.size.pvalue.bs  <- c(x.step.size.pvalue.bs, NA)
  }
}

# ** 2.5.4 euclidian dist C.I. ----
#calculate the SD at each time slice
mv.d.CI <- as.data.frame(mv.ci.r.x.d.bs[1:40,6,1:100])
mv.d.CI$interval <- NA
for(i in 1:length(mv.d.CI[,1])){
  mv.d.CI[i,'interval'] <- sd(mv.d.CI[i,], na.rm=TRUE)
}
mv.ci.r.x.d.bs.f$mv.CI <- mv.d.CI$interval


# 3. qq-plot tests for exponential ----
#https://stats.stackexchange.com/questions/76994/how-do-i-check-if-my-data-fits-an-exponential-distribution

# * 3.1 step size qq-plot test for exponential ----
par(mfrow = c(4,2), mar = c(4,4,1,1), oma = c(2,3,0,0))#, mfg=c(1,1))
# ips step size (r)
r.ips <- as.vector(na.omit(ips.ci.r.x.d$r.step.size[ips.ci.r.x.d$r.step.size>0])) #multivariate step sizes toward the optimum
r.ips.cor <- f.qqexp(r.ips, mainlab = "ips (r)", line = TRUE) #, line = TRUE, step.type = "r") #correlation = 0.9757

# mds step size (r)
par(mfg=c(2,1))
r.mds <- as.vector(na.omit(mds.ci.r.x.d$r.step.size[mds.ci.r.x.d$r.step.size>0])) #multivariate step sizes toward the optimum
r.mds.cor <- f.qqexp(r.mds, mainlab = "mds (r)", line = TRUE, step.type = "r") #correlation = 0.9682

# ptt step size (r)
par(mfg=c(3,1))
r.ptt <- as.vector(na.omit(ptt.ci.r.x.d$r.step.size[ptt.ci.r.x.d$r.step.size>0])) #multivariate step sizes toward the optimum
r.ptt.cor <- f.qqexp(r.ptt, mainlab = "ptt (r)", line = TRUE, step.type = "r") #correlation = 0.9543

# mv step size (r)
par(mfg=c(4,1))
r.mv <- as.vector(na.omit(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size>0])) #multivariate step sizes toward the optimum
r.mv.cor <- f.qqexp(r.mv, mainlab = "mv (r)", line = TRUE, step.type = "r") #correlation = 0.9875

# * 3.2 standardized step size qq-plot test for exponential ----
# ips standardized step size(x)
par(mfg=c(1,2))
x.ips <- as.vector(na.omit(ips.ci.r.x.d$x.step.size[ips.ci.r.x.d$x.step.size>0])) #multivariate step sizes toward the optimum
x.ips.cor <- f.qqexp(x.ips, mainlab = "ips (x)", line = TRUE, step.type = "x") #correlation = 0.9868

# mds standardized step size (x)
par(mfg=c(2,2))
x.mds <- as.vector(na.omit(mds.ci.r.x.d$x.step.size[mds.ci.r.x.d$x.step.size>0])) #multivariate step sizes toward the optimum
x.mds.cor <- f.qqexp(x.mds, mainlab = "mds (x)", line = TRUE, step.type = "x") #correlation = 0.9632

# ptt standardized step size (x)
par(mfg=c(3,2))
x.ptt <- as.vector(na.omit(ptt.ci.r.x.d$x.step.size[ptt.ci.r.x.d$x.step.size>0])) #multivariate step sizes toward the optimum
x.ptt.cor <- f.qqexp(x.ptt, mainlab = "ptt (x)", line = TRUE, step.type = "x") #correlation = 0.9921

# mv standardized step size (x)
par(mfg=c(4,2))
x.mv <- as.vector(na.omit(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size>0])) #multivariate step sizes toward the optimum
x.mv.cor <- f.qqexp(x.mv, mainlab = "mv (x)", line = TRUE, step.type = "x") #correlation = 0.9830
mtext(side=2, text=c('Pelvic Score','Dorsal Spine Number','Touching \n Pterygiophore Num.','Multivariate'),at=c(.90,.65,.40,.13),bty='n', line=0,outer=TRUE)


# 4. S1 table ----
s.table <- as.data.frame(sample.N)
s.table <- cbind(row.names(s.table), s.table)
s.table <- s.table[,c(1,2,3,5)]
colnames(s.table) <- c('Temporal_Sample','ips','mds','ptt')
s.table$Temporal_Sample <- as.integer(s.table$Temporal_Sample)-5000
# * 4.1 FIGURE S1 ----
#Plot the individual traits using different points for each line
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))
#WALK STOPS AT 9750 YRS, row 40
plot(s.table$Temporal_Sample[1:40], s.table[c(1:40),2], type = "b", lty = 3, pch=15, col='#66c2a5', xlab='',ylab='',xaxt='n', yaxt = 'n', las=1, bty='n', ylim = c(0, 125))
lines(s.table$Temporal_Sample[1:40], s.table[1:40,3], type = "b", lty = 3,pch=17, col='#8da0cb')
lines(s.table$Temporal_Sample[1:40], s.table[1:40,4], type = "b", lty = 3,pch=20, col='#fc8d62')
legend(5000,120, legend=c("Pelvic Score", "Dorsal Spine Number", 'Touching Pterygiophore Num.'),
       col=c("#66c2a5", "#8da0cb", '#fc8d62'), lty = 3, pch=c(15,17,20), cex=0.8, bty='n')
#add in the labels
axis(side = 1, at=s.table$Temporal_Sample[1:40], labels = TRUE)
axis(side = 2, at= seq(from = 0, to = max(s.table[1:40, 2:4]), by = 10), labels = TRUE)
mtext('Sample Size',side=2, line = 3)
mtext('Relative time of deposition in years (since start of adaptive walk)',side=1, line=3)

# 5. Confirming log-linear Alpha assumption ----
a <- seq(0,1,by=0.05)
a.rse <- as.data.frame(matrix(data=NA,nrow=length(a),ncol=5));colnames(a.rse)<-c('alpha','ips.rse','mds.rse','ptt.rse','mv.rse');a.rse$alpha<-a

# * 5.1 Prepping the data ----
ips.ci.r.x.d
ips.r.pos.steps = sort(na.omit(ips.ci.r.x.d$r.step.size[which(ips.ci.r.x.d$r.step.size > 0)] ))
ips.r.n = length(ips.r.pos.steps)
ips.r.rank = 1:ips.r.n

mds.r.pos.steps = sort(na.omit(mds.ci.r.x.d$r.step.size[which(mds.ci.r.x.d$r.step.size > 0)] ))
mds.r.n = length(mds.r.pos.steps)
mds.r.rank = 1:mds.r.n

ptt.r.pos.steps = sort(na.omit(ptt.ci.r.x.d$r.step.size[which(ptt.ci.r.x.d$r.step.size > 0)] ))
ptt.r.n = length(ptt.r.pos.steps)
ptt.r.rank = 1:ptt.r.n

mv.r.pos.steps = sort(na.omit(mv.ci.r.x.d$mv.r.step.size[which(mv.ci.r.x.d$mv.r.step.size > 0)] ))
mv.r.n = length(mv.r.pos.steps)
mv.r.rank = 1:mv.r.n

# * 5.2 Finding the RSE ----
for(i in 1:length(a)+1){
  #ips
  ips.r.plotting.position = (ips.r.rank - a[i]) / (ips.r.n - 2*a[i] + 1)
  y.t1 <- sort(na.omit(ips.ci.r.x.d$r.step.size[which(ips.ci.r.x.d$r.step.size > 0)])) ; x.t1 <- -log(1-ips.r.plotting.position) 
  model1 <- lm(y.t1 ~ x.t1)
  ips.rse.t <- summary(model1)
  
  #mds
  mds.r.plotting.position = (mds.r.rank - a[i]) / (mds.r.n - 2*a[i] + 1)
  y.t2 <- sort(na.omit(mds.ci.r.x.d$r.step.size[which(mds.ci.r.x.d$r.step.size > 0)])) ; x.t2 <- -log(1-mds.r.plotting.position) 
  model2 <- lm(y.t2 ~ x.t2)
  mds.rse.t <- summary(model2)
  
  #ptt
  ptt.r.plotting.position = (ptt.r.rank - a[i]) / (ptt.r.n - 2*a[i] + 1)
  y.t3 <- sort(na.omit(ptt.ci.r.x.d$r.step.size[which(ptt.ci.r.x.d$r.step.size > 0)])) ; x.t3 <- -log(1-ptt.r.plotting.position) 
  model3 <- lm(y.t3 ~ x.t3)
  ptt.rse.t <- summary(model3)
  
  #mv
  mv.r.plotting.position = (mv.r.rank - a[i]) / (mv.r.n - 2*a[i] + 1)
  y.t4 <- sort(na.omit(mv.ci.r.x.d$mv.r.step.size[which(mv.ci.r.x.d$mv.r.step.size > 0)])) ; x.t4 <- -log(1-mv.r.plotting.position) 
  model4 <- lm(y.t4 ~ x.t4)
  mv.rse.t <- summary(model4)
  
   #store
  a.rse[i,] <- c(a[i],ips.rse.t$r.squared,mds.rse.t$r.squared,ptt.rse.t$r.squared,mv.rse.t$r.squared)
}

# 5.3 Plotting the R2 Values ----
par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,0,0))
#Plot the individual traits using different points for each line
plot(x = a.rse$alpha, y = a.rse[,2], type = "b", lty = 3, pch=15, col='#66c2a5', xlab='',ylab='', las=1, bty='n', ylim = c(0, 1))
lines(a.rse$alpha, a.rse[,3], type = "b", lty = 3,pch=17, col='#8da0cb')
lines(a.rse$alpha, a.rse[,4], type = "b", lty = 3,pch=20, col='#fc8d62')
lines(a.rse$alpha, a.rse[,5], type = "b", lty = 3,pch=20, col='#fc8d62') #need another color and shape for multivariate
legend(0.2,0.2, legend=c("Pelvic Score", "Dorsal Spine Number", 'Touching Pterygiophore Num.', 'Multivariate'),
       col=c("#66c2a5", "#8da0cb", '#fc8d62','#fc8d62'), lty = 3, pch=c(15,17,20), cex=0.8, bty='n')
#add in the labels
mtext('R-Squared',side=2, line = 3)
mtext('alpha',side=1, line=3)

#conclusion: effect of alpha on R2 estimates is pretty flat from 0 < a < 0.75. Even when alpha has an effect, R2 values remain high. Thus, there is always support for a fit to exponential.

# 6. Potential Lilliefors-Corrected Kolmogorov-Smirnov test for Exponential ----
#https://www.imsbio.co.jp/RGM/R_rdfile?f=KScorrect/man/LcKS.Rd&d=R_CC
# * 6.1 Raw Steps KS Test ----
par(mfrow = c(2,4), mar = c(2,4,1,1), oma = c(2,3,0,0))
r.steps <- list(r.ips, r.mds, r.ptt, r.mv) #store steps
r.names <- c('ips.r','mds.r','ptt.r','mv.r')
r.KS.values <- c()
for (i in 1:length(r.steps)){ #perform LcKS test on steps from each trait
  ks.temp <- LcKS(r.steps[[i]], 'pexp')
  hist(ks.temp$D.sim, main=paste('Histogram of', r.names[i], '| p-value', ks.temp$p.value)) #create a histogram for each trait based on the D.sim
  abline(v = ks.temp$D.obs, lty=2) #mark the D test statistic for this set
  r.KS.values <- c(r.KS.values, ks.temp$p.value) #store p-values for each trait
}

# * 6.2 Standardized Steps KS Test ----
x.steps <- list(x.ips, x.mds, x.ptt, x.mv) #store steps
x.names <- c('ips.x','mds.x','ptt.x','mv.x')
x.KS.values <- c()
for (i in 1:length(x.steps)){ #perform LcKS test on steps from each trait
  ks.temp <- LcKS(x.steps[[i]], 'pexp')
  hist(ks.temp$D.sim, main=paste('Histogram of', x.names[i], '| p-value', ks.temp$p.value)) #create a histogram for each trait based on the D.sim
  abline(v = ks.temp$D.obs, lty=2) #mark the D test statistic for this set
  x.KS.values <- c(x.KS.values, ks.temp$p.value) #store p-values for each trait
}

# end of analysis
#==============================================================================#
