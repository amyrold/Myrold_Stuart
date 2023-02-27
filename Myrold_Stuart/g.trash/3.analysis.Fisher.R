##### Table of Contents #####
#0) Import data and housekeeping
  #0a) Import
  #0b) Optimum definition. From Hunt et al. 2008
#1) Univariate analyses
#2) Multivariate analyses

##### 0) Import Data and housekeeping #####
##### 0a) Import Data #####
#individual level l.series with dorsal spines and touching pterygiophores
#These data are what were analyzed in Bell et al. 2006 and Hunt et al. 2008. 
data.l.sample.sbf <- read.csv(paste(path.data.cleaned.Fisher, "l.series_from_FossilIntervalData_recreatesAppendix1Bell2006.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

# find only lineage II fish. Lineage II fish were the lineage to which evolutionary models were fit to infer evolutionary mode subsequent to the transition.
tapply(X = data.l.sample.sbf$ips, INDEX = data.l.sample.sbf$interval, FUN = mean, na.rm = T) # pelvic score jumps at 5000 years, suggesting a mix of lineage I and lineage II fish.
data.l.sample.sbf$ips[data.l.sample.sbf$interval == 5000] #mix of 1.0s and 3.0s. Keep the 3.0s as lineage II fish.

#subsample to only the lineage II fish.
data.l.sample.sbf.final <- rbind(data.l.sample.sbf[which(data.l.sample.sbf$ips == 3 & data.l.sample.sbf$interval == 5000), ], data.l.sample.sbf[data.l.sample.sbf$interval > 5000, ])
data.l.sample.sbf.final <- data.l.sample.sbf.final[-which(data.l.sample.sbf.final$interval == 23250) , ]

#means by sample. Useful for univariate and multivariate analysis.
sample.means <- apply(data.l.sample.sbf.final[ , c(2, 3, 5)], MAR = 2, FUN = tapply, data.l.sample.sbf.final$interval, mean, na.rm = TRUE)

#data were logged in Bell et al. 2006
data.l.sample.sbf.final.logged <- data.frame(data.l.sample.sbf.final$interval, data.l.sample.sbf.final$specID, log(data.l.sample.sbf.final[, 2:5]+1)); colnames(data.l.sample.sbf.final.logged) <- c("interval", "specID", "log.ips", "log.mds", "log.mpt", "log.ptt")
logged.sample.means <- apply(data.l.sample.sbf.final.logged[ , c(3,4,6)], MAR = 2, FUN = tapply, data.l.sample.sbf.final$interval, mean, na.rm = TRUE)


##### 0b) Optimum definition #####
#based on the theta estimates from Hunt et al. 2008
#optimum number of dorsal spines
#ln(x + 1) 
theta.mds.transformed <- 0.79
theta.mds <- exp(0.79) - 1

#optimum number of touching pterygiophores
#ln(x + 1)
theta.ptt.transformed <- 0.82
theta.ptt <- exp(0.82) - 1 + 1 #the plus one here gets optimum in line with actual data. 
##### NB: something is wrong here. optimum is below all the data but one point. #####
# Note well: sample means do not re-create Appendix 1 from Bell et al. 2006, which was used directly for analyses in Hunt et al. 2008

#optimum number pelvic score
#ln(x + 1)
theta.ips.transformed <- 0.73
theta.ips <- exp(0.73) - 1

optimum <- c(theta.ips, theta.mds, theta.ptt)

# end 0) 

##### 1) Univariate analysis #####
#ips
ips.t0 <- sample.means[1:(nrow(sample.means)-1) , 1]
ips.t1 <- sample.means[2:nrow(sample.means), 1]
ips.movements <- ips.t1-ips.t0
ips.distance.from.optimum <- sample.means[ , 1] - optimum[1]

par(mfrow = c(2, 1), mar = c(4,4,1,0), oma = c(0, 0, 0, 0))
#trait.
plot(y = sample.means[ , 1], x = seq(5000-5000, 23000-5000, by = 250), pch = 16, bty = "n", xaxt = "n", xlab = "Years since lineage II appeared", ylab = "Pelvic score", ylim = c(0, 3), type = "b", lty = 3)
axis(side = 1, at = seq(0, 18000, by = 1000))
segments(x0 = 0, x1 = 18000, y0 = optimum[1], y1 = optimum[1], lty = 2)

#change. Expect mostly negative changes. larger earlier
color.ips <- vector()
for(i in 1:length(ips.movements)) if(ips.movements[i] < 0) color.ips <- c(color.ips, "red") else color.ips <- c(color.ips, "black")
plot(y = ips.movements, x = seq(5000-5000, 22750-5000, by = 250), pch = 16, bty = "n", xaxt = "n", xlab = "Start of time interval (in years since Lineage II appeared)", ylab = "Pelvic score change", col = color.ips, las = 1, lty = 3)
abline(a = 0, b = 0)
segments(x0 = seq(5000-5000, 22750-5000, by = 250), x1 = seq(5000-5000, 22750-5000, by = 250), y0 = rep(0, length(seq(5000-5000, 22750-5000, by = 250))), y1 = ips.movements)
axis(side = 1, at = seq(0, 18000, by = 1000))

#plot distribution steps toward optimum irrespective of order.
#to do: limit only to Hunts' t1/2 * 2
ips.movements.neg <- ips.movements[which(ips.movements<0)]
hist(abs(ips.movements.neg[which(as.numeric(names(ips.movements.neg)) %in% seq(6750, 14250, by = 250))]))

#plot distribution of steps toward optimum ordered.
ips.movements.neg.temporal <- ips.movements.neg[which(as.numeric(names(ips.movements.neg)) %in% seq(6750, 14250, by = 250))]
plot(abs(ips.movements.neg.temporal) ~ as.numeric(names(ips.movements.neg.temporal)), pch = 16 )

#rapid adaptation from standing variation (small early steps) while population waits for large effect mutations

#mds
mds.t0 <- sample.means[1:(nrow(sample.means)-1) , 2]
mds.t1 <- sample.means[2:nrow(sample.means), 2]
mds.movements <- mds.t1-mds.t0

par(mfrow = c(2, 1), mar = c(4,4,1,0), oma = c(0, 0, 0, 0))
#trait.
plot(y = sample.means[ , 2], x = seq(5000-5000, 23000-5000, by = 250), pch = 16, bty = "n", xaxt = "n", xlab = "Years since lineage II appeared", ylab = "Dorsal Spine Number", ylim = c(0, 3), type = "b", lty = 3)
axis(side = 1, at = seq(0, 18000, by = 1000))
segments(x0 = 0, x1 = 18000, y0 = optimum[2], y1 = optimum[2], lty = 2)

#change. Expect mostly negative changes. larger earlier
color.mds <- vector()
for(i in 1:length(mds.movements)) if(mds.movements[i] < 0) color.mds <- c(color.mds, "red") else color.mds <- c(color.mds, "black")
plot(y = mds.movements, x = seq(5000-5000, 22750-5000, by = 250), pch = 16, bty = "n", xaxt = "n", xlab = "Start of time interval (in years since Lineage II appeared)", ylab = "Dorsal Spine number change", col = color.mds, las = 1, type = "b", lty = 3)
abline(a = 0, b = 0)
axis(side = 1, at = seq(0, 18000, by = 1000))

#ptt
ptt.t0 <- sample.means[1:(nrow(sample.means)-1) , 3]
ptt.t1 <- sample.means[2:nrow(sample.means), 3]
ptt.movements <- ptt.t1-ptt.t0

par(mfrow = c(2, 1), mar = c(4,4,1,0), oma = c(0, 0, 0, 0))
#trait.
plot(y = sample.means[ , 3], x = seq(5000-5000, 23000-5000, by = 250), pch = 16, bty = "n", xaxt = "n", xlab = "Years since lineage II appeared", ylab = "Touching pterygiophore number", ylim = c(0,6), type = "b", lty = 3)
axis(side = 1, at = seq(0, 18000, by = 1000))
segments(x0 = 0, x1 = 18000, y0 = optimum[3], y1 = optimum[3], lty = 2) #something wrong with this optimum.

#change. Expect mostly negative changes. larger earlier
color.ptt <- vector()
for(i in 1:length(ptt.movements)) if(ptt.movements[i] < 0) color.ptt <- c(color.ptt, "red") else color.ptt <- c(color.ptt, "black")
plot(y = ptt.movements, x = seq(5000-5000, 22750-5000, by = 250), pch = 16, bty = "n", xaxt = "n", xlab = "Start of time interval (in years since Lineage II appeared)", ylab = "ptt number change", col = color.ptt, las = 1, type = "b", lty = 3)
abline(a = 0, b = 0)
axis(side = 1, at = seq(0, 18000, by = 1000))


##### 2) Multivariate analysis  #####
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



