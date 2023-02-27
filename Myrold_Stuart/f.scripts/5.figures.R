# Table of Contents ----
#   1. Data visualization. Phyletic evolution, adaptive steps, step distributions
#     1.1 ips. confidence interval method
#     1.2 mds confidence interval method
#     1.3 ptt confidence interval method
#     1.4 Multivariate Data Visualization
#   2. log-linear plots of standardized step sizes
#     2.1 setting up the parameters
#     2.2 creating the plots
#   3. Fitting distributions against data
#     3.1 Fitting step size (r) against exponential dist
#     3.2 Fitting standardized step size (x) against exponential dist
#   4. using qqplot methods described in the comments below
#     4.1 raw mv step size (r)
#     4.2 standardized mv step size (x)
# end Table of Contents
#==============================================================================#

# 1. FIGURE 1: Data visualization. Phyletic evolution, adaptive steps, step distributions ----
#xlabel = "Relative time of deposition (in years since Lineage II appeared)"
xlabel = "Relative time of deposition in years since lineage II appeared"

# * 1.1 ips. confidence interval method ----
ips.ci.r.x.d
par(mfcol = c(3, 4), mar = c(4,4,1,0), oma = c(2, 2 , 2, 1))

#ips trend
plot(y = ips.ci.r.x.d$sample.mean.step.start, x = ips.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", ylim = c(0, 3), type = "b", lty = 3, las = 1)
axis(side = 1, at = ips.ci.r.x.d$linII.time)
segments(x0 = ips.ci.r.x.d$linII.time[1], x1 = ips.ci.r.x.d$linII.time[length(ips.ci.r.x.d$linII.time)], y0 = optimum[1], y1 = optimum[1], lty = 2)
segments(x0 = ips.ci.r.x.d$linII.time, x1 = ips.ci.r.x.d$linII.time, y0 = ips.ci.r.x.d$sample.mean.step.start, y1 = ips.ci.r.x.d$sample.mean.step.start - data.sbf.f$ips.CI[data.sbf.f$linII.time %in% ips.ci.r.x.d$linII.time])

#ips step sizes, r
color.ips <- vector()
for(i in 1:(length(ips.ci.r.x.d$r.step.size)-1)) if(ips.ci.r.x.d$r.step.size[i] < 0) color.ips <- c(color.ips, "red") else color.ips <- c(color.ips, "black")
plot(y = ips.ci.r.x.d$r.step.size, x = ips.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.ips, las = 1, lty = 3, ylim = c(-0.6, 0.6))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = ips.ci.r.x.d$linII.time, x1 = ips.ci.r.x.d$linII.time, y0 = rep(0, length(ips.ci.r.x.d$r.step.size)), y1 = ips.ci.r.x.d$r.step.size, col = "grey50")
axis(side = 1, at = ips.ci.r.x.d$linII.time)
#ips r stepsizes, bootstrap mean and bootstrap standard deviation. 
segments(x0 = ips.boot.t$linII.time + 50, x1 = ips.boot.t$linII.time + 100, y0 = ips.boot.t$r.step.size.bs, y1 = ips.boot.t$r.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = ips.boot.t$linII.time + 75, x1 = ips.boot.t$linII.time + 75, y0 = ips.boot.t$r.step.size.bs-ips.boot.t$r.step.size.ci.bs, y1 = ips.boot.t$r.step.size.bs+ips.boot.t$r.step.size.ci.bs, col = "grey60",lwd=1)

#ips step sizes, x
color.ips <- vector()
for(i in 1:(length(ips.ci.r.x.d$x.step.size)-1)) if(ips.ci.r.x.d$x.step.size[i] < 0) color.ips <- c(color.ips, "red") else color.ips <- c(color.ips, "black")
plot(y = ips.ci.r.x.d$x.step.size, x = ips.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.ips, las = 1, lty = 3, ylim = c(-0.7, 0.7))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = ips.ci.r.x.d$linII.time, x1 = ips.ci.r.x.d$linII.time, y0 = rep(0, length(ips.ci.r.x.d$x.step.size)), y1 = ips.ci.r.x.d$x.step.size, col = "grey50")
axis(side = 1, at = ips.ci.r.x.d$linII.time)
text(x = 3500, y = 0.5, label = "toward optimum", col = "black", adj = 0) ; text(x = 3500, y = -0.5, label = "away from optimum", col = "red", adj = 0)
#ips x stepsizes, bootstrap mean and bootstrap standard deviation. 
segments(x0 = ips.boot.t$linII.time + 50, x1 = ips.boot.t$linII.time + 100, y0 = ips.boot.t$x.step.size.bs, y1 = ips.boot.t$x.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = ips.boot.t$linII.time + 75, x1 = ips.boot.t$linII.time + 75, y0 = ips.boot.t$x.step.size.bs-ips.boot.t$x.step.size.ci.bs, y1 = ips.boot.t$x.step.size.bs+ips.boot.t$r.step.size.ci.bs, col = "grey60",lwd=1)

# * 1.2 mds confidence interval method ----
mds.ci.r.x.d

#mds trend
plot(y = mds.ci.r.x.d$sample.mean.step.start, x = mds.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", ylim = c(0, 3), type = "b", lty = 3, las = 1)
axis(side = 1, at = mds.ci.r.x.d$linII.time)
segments(x0 = mds.ci.r.x.d$linII.time[1], x1 = mds.ci.r.x.d$linII.time[length(mds.ci.r.x.d$linII.time)], y0 = optimum[2], y1 = optimum[2], lty = 2)
segments(x0 = mds.ci.r.x.d$linII.time, x1 = mds.ci.r.x.d$linII.time, y0 = mds.ci.r.x.d$sample.mean.step.start, y1 = mds.ci.r.x.d$sample.mean.step.start - data.sbf.f$mds.CI[data.sbf.f$linII.time %in% mds.ci.r.x.d$linII.time])

#mds step sizes, r
color.mds <- vector()
for(i in 1:(length(mds.ci.r.x.d$r.step.size)-1)) if(mds.ci.r.x.d$r.step.size[i] < 0) color.mds <- c(color.mds, "red") else color.mds <- c(color.mds, "black")
plot(y = mds.ci.r.x.d$r.step.size, x = mds.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.mds, las = 1, lty = 3, ylim = c(-0.6, 0.6))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = mds.ci.r.x.d$linII.time, x1 = mds.ci.r.x.d$linII.time, y0 = rep(0, length(mds.ci.r.x.d$r.step.size)), y1 = mds.ci.r.x.d$r.step.size, col = "grey50")
axis(side = 1, at = mds.ci.r.x.d$linII.time)
#mds r stepsizes, bootstrap mean and bootstrap standard deviation. 
segments(x0 = mds.boot.t$linII.time + 25, x1 = mds.boot.t$linII.time + 50, y0 = mds.boot.t$r.step.size.bs, y1 = mds.boot.t$r.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = mds.boot.t$linII.time + 37.5, x1 = mds.boot.t$linII.time + 37.5, y0 = mds.boot.t$r.step.size.bs-mds.boot.t$r.step.size.ci.bs, y1 = mds.boot.t$r.step.size.bs+mds.boot.t$r.step.size.ci.bs, col = "grey60",lwd=1)

#mds standardized step sizes, x
color.mds <- vector()
for(i in 1:(length(mds.ci.r.x.d$x.step.size)-1)) if(mds.ci.r.x.d$x.step.size[i] < 0) color.mds <- c(color.mds, "red") else color.mds <- c(color.mds, "black")
plot(y = mds.ci.r.x.d$x.step.size, x = mds.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.mds, las = 1, lty = 3, ylim = c(-0.7, 0.7))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = mds.ci.r.x.d$linII.time, x1 = mds.ci.r.x.d$linII.time, y0 = rep(0, length(mds.ci.r.x.d$x.step.size)), y1 = mds.ci.r.x.d$x.step.size, col = "grey50")
axis(side = 1, at = mds.ci.r.x.d$linII.time)
#text(x = 0, y = 0.5, label = "toward optimum", col = "black", adj = 0) ; text(x = 0, y = -0.5, label = "away from optimum", col = "red", adj = 0)
#mds x stepsizes, bootstrap mean and bootstrap standard deviation. 
segments(x0 = mds.boot.t$linII.time + 25, x1 = mds.boot.t$linII.time + 50, y0 = mds.boot.t$x.step.size.bs, y1 = mds.boot.t$x.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = mds.boot.t$linII.time + 37.5, x1 = mds.boot.t$linII.time + 37.5, y0 = mds.boot.t$x.step.size.bs-mds.boot.t$x.step.size.ci.bs, y1 = mds.boot.t$x.step.size.bs+mds.boot.t$x.step.size.ci.bs, col = "grey60",lwd=1)

# * 1.3 ptt confidence interval method ----
ptt.ci.r.x.d

#ptt trend
plot(y = ptt.ci.r.x.d$sample.mean.step.start, x = ptt.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", ylim = c(0, max(ptt.ci.r.x.d$sample.mean.step.start)), type = "b", lty = 3, las = 1)
axis(side = 1, at = ptt.ci.r.x.d$linII.time)
segments(x0 = ptt.ci.r.x.d$linII.time[1], x1 = ptt.ci.r.x.d$linII.time[length(ptt.ci.r.x.d$linII.time)], y0 = optimum[3], y1 = optimum[3], lty = 2)
segments(x0 = ptt.ci.r.x.d$linII.time, x1 = ptt.ci.r.x.d$linII.time, y0 = ptt.ci.r.x.d$sample.mean.step.start, y1 = ptt.ci.r.x.d$sample.mean.step.start - data.sbf.f$ptt.CI[data.sbf.f$linII.time %in% ptt.ci.r.x.d$linII.time])

#ptt step sizes, r
color.ptt <- vector()
for(i in 1:(length(ptt.ci.r.x.d$r.step.size)-1)) if(ptt.ci.r.x.d$r.step.size[i] < 0) color.ptt <- c(color.ptt, "red") else color.ptt <- c(color.ptt, "black")
plot(y = ptt.ci.r.x.d$r.step.size, x = ptt.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.ptt, las = 1, lty = 3, ylim = c(-.4, .8))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = ptt.ci.r.x.d$linII.time, x1 = ptt.ci.r.x.d$linII.time, y0 = rep(0, length(ptt.ci.r.x.d$r.step.size)), y1 = ptt.ci.r.x.d$r.step.size, col = "grey50")
axis(side = 1, at = ptt.ci.r.x.d$linII.time)
#ptt r stepsizes, bootstrap mean and bootstrap standard deviation. 
segments(x0 = ptt.boot.t$linII.time + 25, x1 = ptt.boot.t$linII.time + 50, y0 = ptt.boot.t$r.step.size.bs, y1 = ptt.boot.t$r.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = ptt.boot.t$linII.time + 37.5, x1 = ptt.boot.t$linII.time + 37.5, y0 = ptt.boot.t$r.step.size.bs-ptt.boot.t$r.step.size.ci.bs, y1 = ptt.boot.t$r.step.size.bs+ptt.boot.t$r.step.size.ci.bs, col = "grey60",lwd=1)

#ptt standardized step sizes, x
color.ptt <- vector()
for(i in 1:(length(ptt.ci.r.x.d$x.step.size)-1)) if(ptt.ci.r.x.d$x.step.size[i] < 0) color.ptt <- c(color.ptt, "red") else color.ptt <- c(color.ptt, "black")
plot(y = ptt.ci.r.x.d$x.step.size, x = ptt.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.ptt, las = 1, lty = 3, ylim = c(-0.7, 0.7))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = ptt.ci.r.x.d$linII.time, x1 = ptt.ci.r.x.d$linII.time, y0 = rep(0, length(ptt.ci.r.x.d$x.step.size)), y1 = ptt.ci.r.x.d$x.step.size, col = "grey50")
axis(side = 1, at = ptt.ci.r.x.d$linII.time)
#text(x = 0, y = 0.5, label = "toward optimum", col = "black", adj = 0) ; text(x = 0, y = -0.5, label = "away from optimum", col = "red", adj = 0)
#ptt x stepsizes, bootstrap mean and bootstrap standard deviation. 
segments(x0 = ptt.boot.t$linII.time + 25, x1 = ptt.boot.t$linII.time + 50, y0 = ptt.boot.t$x.step.size.bs, y1 = ptt.boot.t$x.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = ptt.boot.t$linII.time + 37.5, x1 = ptt.boot.t$linII.time + 37.5, y0 = ptt.boot.t$x.step.size.bs-ptt.boot.t$x.step.size.ci.bs, y1 = ptt.boot.t$x.step.size.bs+ptt.boot.t$x.step.size.ci.bs, col = "grey60",lwd=1)

# * 1.4 Multivariate Data Visualization ----
# Approach to optimum and step sizes
#dist.to.optimum trend
plot(y = mv.ci.r.x.d$euclid.dist.to.opt.step.start, x = mv.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = '', ylab = "", ylim = c(0, 4), type = "b", lty = 3, las=1)
mtext(text='Euclidean distance to optimum',line=2, side=2,cex=0.6)
abline(h = 0, col = "grey80")
axis(side = 1, at = mv.ci.r.x.d$linII.time)

segments(x0 = mv.ci.r.x.d$linII.time, x1 = mv.ci.r.x.d$linII.time, y0 = mv.ci.r.x.d$euclid.dist.to.opt.step.start, y1 = mv.ci.r.x.d$euclid.dist.to.opt.step.start - mv.ci.r.x.d.bs.f$mv.CI[mv.ci.r.x.d.bs.f$linII.time %in% mv.ci.r.x.d$linII.time])

#  multivariate step sizes, r
color.mv <- vector()
for(i in 1:(length(mv.ci.r.x.d$mv.r.step.size)-1)) if(mv.ci.r.x.d$mv.r.step.size[i] < 0) color.mv <- c(color.mv, "red") else color.mv <- c(color.mv, "black")
plot(y = mv.ci.r.x.d$mv.r.step.size, x = mv.ci.r.x.d$linII.time, pch = 16, bty = "n", xaxt = "n", xlab = "", ylab = "", col = color.mv, las = 1, lty = 3, ylim = c(-0.6, 0.6))
abline(a = 0, b = 0, col = "grey80")
segments(x0 = mv.ci.r.x.d$linII.time, x1 = mv.ci.r.x.d$linII.time, y0 = rep(0, length(mv.ci.r.x.d$linII.time)), y1 = mv.ci.r.x.d$mv.r.step.size, col = "grey50")
#mv r stepsizes, bootstrap mean and bootstrap confidence interval. 
segments(x0 = mv.ci.r.x.d.bs.f$linII.time + 50, x1 = mv.ci.r.x.d.bs.f$linII.time + 125, y0 = mv.ci.r.x.d.bs.f$mv.r.step.size.bs, y1 = mv.ci.r.x.d.bs.f$mv.r.step.size.bs, col = "grey60")
segments(x0 = mv.ci.r.x.d.bs.f$linII.time + 87.5, x1 = mv.ci.r.x.d.bs.f$linII.time + 87.5, y0 = mv.ci.r.x.d.bs.f$mv.r.step.size.bs-mv.ci.r.x.d.bs.f$r.step.size.ci.bs, y1 = mv.ci.r.x.d.bs.f$mv.r.step.size.bs+mv.ci.r.x.d.bs.f$r.step.size.ci.bs, col = "grey60")
#
axis(side = 1, at = mv.ci.r.x.d$linII.time)

# multivariate standardized step sizes, x
color.mv <- vector()
for(i in 1:(length(mv.ci.r.x.d$mv.x.step.size)-1)) if(mv.ci.r.x.d$mv.x.step.size[i] < 0) color.mv <- c(color.mv, "red") else color.mv <- c(color.mv, "black")
plot(y = mv.ci.r.x.d$mv.x.step.size, x = mv.ci.r.x.d$linII.time, pch = 16, bty = "n", yaxt = 'n', xaxt = "n", xlab = '', ylab = "", col = color.mv, las = 1, lty = 3, ylim = c(-0.9, 0.6))
abline(a = 0, b = 0, col = "grey80",lwd=1.5)
segments(x0 = mv.ci.r.x.d$linII.time, x1 = mv.ci.r.x.d$linII.time, y0 = rep(0, length(mv.ci.r.x.d$linII.time)), y1 = mv.ci.r.x.d$mv.x.step.size, col = "grey50")
#mv x stepsizes, bootstrap mean and bootstrap confidence interval. 
segments(x0 = mv.ci.r.x.d.bs.f$linII.time + 50, x1 = mv.ci.r.x.d.bs.f$linII.time + 125, y0 = mv.ci.r.x.d.bs.f$mv.x.step.size.bs, y1 = mv.ci.r.x.d.bs.f$mv.x.step.size.bs, col = "grey60",lwd=1.5)
segments(x0 = mv.ci.r.x.d.bs.f$linII.time + 87.5, x1 = mv.ci.r.x.d.bs.f$linII.time + 87.5, y0 = mv.ci.r.x.d.bs.f$mv.x.step.size.bs-mv.ci.r.x.d.bs.f$x.step.size.ci.bs, y1 = mv.ci.r.x.d.bs.f$mv.x.step.size.bs+mv.ci.r.x.d.bs.f$x.step.size.ci.bs, col = "grey60",lwd=1)
#
axis(side = 1, at = mv.ci.r.x.d$linII.time)
axis(side = 2, at = seq(-1,1, by=0.2),las=1)

mtext(text = xlabel, side =1, outer = TRUE)
mtext(text = c('Standardized Step Size (x)','Step Size (r)','Trait Mean (95% C.I.)'), at=c(.19,.53,.85), line=0,side=2,outer=TRUE)
mtext(text = c('Pelvic Score (ips)','Dorsal Spine Num. (mds)','Touching Pterygiophore Num. (ptt)','Multivariate (mv)'), at=c(.14,.39,.64,.89 ), line=0,side =3, outer=TRUE)
#==============================================================================#

#  2. log-linear plots of standardized step sizes ----
#https://stats.stackexchange.com/questions/76994/how-do-i-check-if-my-data-fits-an-exponential-distribution
#An exponential distribution will plot as a straight line against −ln(1−plotting position) where plotting position is (rank −𝑎)/(𝑛−2𝑎+1), rank is 1 for lowest value, 𝑛 is sample size, and popular choices for 𝑎 include 1/2. That gives an informal test which can be as or more useful than any formal test. – Nick Cox Nov 19 '13 at 13:42 

# * 2.1 setting up the parameters ----
# ** 2.1.1 raw step size ----
a = 1/2
#ips
ips.ci.r.x.d
ips.r.pos.steps = sort(na.omit(ips.ci.r.x.d$r.step.size[which(ips.ci.r.x.d$r.step.size > 0)] ))
ips.r.n = length(ips.r.pos.steps)
ips.r.rank = 1:ips.r.n
ips.r.plotting.position = (ips.r.rank - a) / (ips.r.n - 2*a + 1)

#mds
mds.r.pos.steps = sort(na.omit(mds.ci.r.x.d$r.step.size[which(mds.ci.r.x.d$r.step.size > 0)] ))
mds.r.n = length(mds.r.pos.steps)
mds.r.rank = 1:mds.r.n
mds.r.plotting.position = (mds.r.rank - a) / (mds.r.n - 2*a + 1)

#ptt
ptt.r.pos.steps = sort(na.omit(ptt.ci.r.x.d$r.step.size[which(ptt.ci.r.x.d$r.step.size > 0)] ))
ptt.r.n = length(ptt.r.pos.steps)
ptt.r.rank = 1:ptt.r.n
ptt.r.plotting.position = (ptt.r.rank - a) / (ptt.r.n - 2*a + 1)

#mv
mv.r.pos.steps = sort(na.omit(mv.ci.r.x.d$mv.r.step.size[which(mv.ci.r.x.d$mv.r.step.size > 0)] ))
mv.r.n = length(mv.r.pos.steps)
mv.r.rank = 1:mv.r.n
mv.r.plotting.position = (mv.r.rank - a) / (mv.r.n - 2*a + 1)

# ** 2.1.2 standardized step size ----
a = 1/2
#ips
ips.ci.r.x.d
ips.pos.steps = sort(na.omit(ips.ci.r.x.d$x.step.size[which(ips.ci.r.x.d$x.step.size > 0)] ))
ips.n = length(ips.pos.steps)
ips.rank = 1:ips.n
ips.plotting.position = (ips.rank - a) / (ips.n - 2*a + 1)

#mds
mds.pos.steps = sort(na.omit(mds.ci.r.x.d$x.step.size[which(mds.ci.r.x.d$x.step.size > 0)] ))
mds.n = length(mds.pos.steps)
mds.rank = 1:mds.n
mds.plotting.position = (mds.rank - a) / (mds.n - 2*a + 1)

#ptt
ptt.pos.steps = sort(na.omit(ptt.ci.r.x.d$x.step.size[which(ptt.ci.r.x.d$x.step.size > 0)] ))
ptt.n = length(ptt.pos.steps)
ptt.rank = 1:ptt.n
ptt.plotting.position = (ptt.rank - a) / (ptt.n - 2*a + 1)

#mv
mv.pos.steps = sort(na.omit(mv.ci.r.x.d$mv.x.step.size[which(mv.ci.r.x.d$mv.x.step.size > 0)] ))
mv.n = length(mv.pos.steps)
mv.rank = 1:mv.n
mv.plotting.position = (mv.rank - a) / (mv.n - 2*a + 1)

# * 2.2 creating the plots ----
# ** 2.2.1 FIGURE S2: raw step size ----
par(mfrow = c(4,1), mar = c(2,4,1,1), oma = c(2,3,0,0))

#ips
plot(x = -log(1-ips.r.plotting.position), y = sort(na.omit(ips.ci.r.x.d$r.step.size[which(ips.ci.r.x.d$r.step.size > 0)])), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Pelvic Score", type = "b", lty = "dashed", ylim = c(0, 0.8), xlim = c(0, 4))
y.t1 <- sort(na.omit(ips.ci.r.x.d$r.step.size[which(ips.ci.r.x.d$r.step.size > 0)])) ; x.t1 <- -log(1-ips.r.plotting.position) 
model1 <- lm(y.t1 ~ x.t1)
summary(model1)
max.x.1 <- max(-log(1-ips.r.plotting.position)) ; y.1 <- max.x.1 * 0.1796 + 0.0133
segments(x0 = 0, y0 = 0.0133, x1 = max.x.1, y1 = y.1)
text(x = 0, y = 0.4, label = expression(italic("R"^2) * " = 0.94; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model, suggesting exponentially distributed data.

#mds
plot(x = -log(1-mds.r.plotting.position), y = sort(na.omit(mds.ci.r.x.d$r.step.size[which(mds.ci.r.x.d$r.step.size > 0)])), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Dorsal Spine Number", type = "b", lty = "dashed", ylim = c(0, 0.5), xlim = c(0, 4))
y.t2 <- sort(na.omit(mds.ci.r.x.d$r.step.size[which(mds.ci.r.x.d$r.step.size > 0)])) ; x.t2 <- -log(1-mds.r.plotting.position) 
model2 <- lm(y.t2 ~ x.t2)
summary(model2)
max.x.2 <- max(-log(1-mds.r.plotting.position)) ; y.2 <- max.x.2 * 0.1784 + 0.0558
segments(x0 = 0, y0 = 0.0558, x1 = max.x.2, y1 = y.2)
text(x = 0, y = 0.4, label = expression(italic("R"^2) * " = 0.90; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model, 

#ptt
plot(x = -log(1-ptt.r.plotting.position), y = sort(na.omit(ptt.ci.r.x.d$r.step.size[which(ptt.ci.r.x.d$r.step.size > 0)])), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Touching Pterygiophore Num.", type = "b", lty = "dashed", ylim = c(0, 0.8), xlim = c(0, 4))
y.t3 <- sort(na.omit(ptt.ci.r.x.d$r.step.size[which(ptt.ci.r.x.d$r.step.size > 0)])) ; x.t3 <- -log(1-ptt.r.plotting.position) 
model3 <- lm(y.t3 ~ x.t3)
summary(model3)
max.x.3 <- max(-log(1-ptt.r.plotting.position)) ; y.3 <- max.x.3 * 0.3072 + 0.09198
segments(x0 = 0, y0 = 0.09198, x1 = max.x.3, y1 = y.3)
text(x = 0, y = 0.4, label = expression(italic("R"^2) * " = 0.88; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model,

#mv
plot(x = -log(1-mv.r.plotting.position), y = sort(na.omit(mv.ci.r.x.d$mv.r.step.size[which(mv.ci.r.x.d$mv.r.step.size > 0)] )), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Multivariate", type = "b", lty = "dashed", ylim = c(0, 0.6), xlim = c(0, 4))
y.t4 <- sort(na.omit(mv.ci.r.x.d$mv.r.step.size[which(mv.ci.r.x.d$mv.r.step.size > 0)])) ; x.t4 <- -log(1-mv.r.plotting.position) 
model4 <- lm(y.t4 ~ x.t4)
summary(model4)
max.x.4 <- max(-log(1-mv.r.plotting.position)) ; y.4 <- max.x.4 * 0.1375 + 0.06441
segments(x0 = 0, y0 = 0.06441, x1 = max.x.4, y1 = y.4)
text(x = 0, y = 0.5, label = expression(italic("R"^2) * " = 0.95; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model,

#outer margin text
mtext(side = 2, line = 1, text = "RAW STEP SIZE (x)", outer = TRUE)
mtext(side = 1, line = 1, text = "PLOTTING POSITION", outer = TRUE)

# ** 2.2.2 FIGURE 3: standardized step size ----
par(mfrow = c(4,1), mar = c(2,4,1,1), oma = c(2,3,0,0))

#ips
plot(x = -log(1-ips.plotting.position), y = sort(na.omit(ips.ci.r.x.d$x.step.size[which(ips.ci.r.x.d$x.step.size > 0)])), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Pelvic Score", type = "b", lty = "dashed", ylim = c(0, 0.5), xlim = c(0, 4))
y.t1 <- sort(na.omit(ips.ci.r.x.d$x.step.size[which(ips.ci.r.x.d$x.step.size > 0)])) ; x.t1 <- -log(1-ips.plotting.position) 
model1 <- lm(y.t1 ~ x.t1)
summary(model1)
max.x.1 <- max(-log(1-ips.plotting.position)) ; y.1 <- max.x.1 * 0.1007 + 0.0101
segments(x0 = 0, y0 = 0.0101, x1 = max.x.1, y1 = y.1)
text(x = 0, y = 0.4, label = expression(italic("R"^2) * " = 0.97; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model, suggesting exponentially distributed data.

#mds
plot(x = -log(1-mds.plotting.position), y = sort(na.omit(mds.ci.r.x.d$x.step.size[which(mds.ci.r.x.d$x.step.size > 0)])), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Dorsal Spine Number", type = "b", lty = "dashed", ylim = c(0, 0.5), xlim = c(0, 4))
y.t2 <- sort(na.omit(mds.ci.r.x.d$x.step.size[which(mds.ci.r.x.d$x.step.size > 0)])) ; x.t2 <- -log(1-mds.plotting.position) 
model2 <- lm(y.t2 ~ x.t2)
summary(model2)
max.x.2 <- max(-log(1-mds.plotting.position)) ; y.2 <- max.x.2 * 0.1163 + 0.0033
segments(x0 = 0, y0 = 0.0033, x1 = max.x.2, y1 = y.2)
text(x = 0, y = 0.4, label = expression(italic("R"^2) * " = 0.87; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model, 

#ptt
plot(x = -log(1-ptt.plotting.position), y = sort(na.omit(ptt.ci.r.x.d$x.step.size[which(ptt.ci.r.x.d$x.step.size > 0)])), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Touching Pterygiophore Number", type = "b", lty = "dashed", ylim = c(0, 0.5), xlim = c(0, 4))
y.t3 <- sort(na.omit(ptt.ci.r.x.d$x.step.size[which(ptt.ci.r.x.d$x.step.size > 0)])) ; x.t3 <- -log(1-ptt.plotting.position) 
model3 <- lm(y.t3 ~ x.t3)
summary(model3)
max.x.3 <- max(-log(1-ptt.plotting.position)) ; y.3 <- max.x.3 * 0.18631 - 0.01162
segments(x0 = 0, y0 = 0.01162, x1 = max.x.3, y1 = y.3)
text(x = 0, y = 0.4, label = expression(italic("R"^2) * " = 0.99; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model,

#mv
plot(x = -log(1-mv.plotting.position), y = sort(na.omit(mv.ci.r.x.d$mv.x.step.size[which(mv.ci.r.x.d$mv.x.step.size > 0)] )), bty = "n", pch = 16, las = 1, xlab = "", ylab = "Multivariate", type = "b", lty = "dashed", ylim = c(0, 0.6), xlim = c(0, 4))
y.t4 <- sort(na.omit(mv.ci.r.x.d$mv.x.step.size[which(mv.ci.r.x.d$mv.x.step.size > 0)])) ; x.t4 <- -log(1-mv.plotting.position) 
model4 <- lm(y.t4 ~ x.t4)
summary(model4)
max.x.4 <- max(-log(1-mv.plotting.position)) ; y.4 <- max.x.4 * 0.1507 - 0.0068
segments(x0 = 0, y0 = 0.0068, x1 = max.x.4, y1 = y.4)
text(x = 0, y = 0.5, label = expression(italic("R"^2) * " = 0.96; " * italic(p) * "< 0.0001"), adj = 0) #explained pretty well by a linear model,

#outer margin text
mtext(side = 2, line = 1, text = "STANDARDIZED STEP SIZE (x)", outer = TRUE)
mtext(side = 1, line = 1, text = "PLOTTING POSITION", outer = TRUE)

#==============================================================================#

# 3. Fitting distributions against data ----
# * 3.1 FIGURE S3: Fitting step size (r) against exponential dist ----
par(mfrow = c(4,1), mar = c(4, 1 ,1,0), oma = c(1, 3.25, 0, 1))
# ips step size histogram
ips.r.rate <- unname(eexp(ips.ci.r.x.d$r.step.size[ips.ci.r.x.d$r.step.size > 0])$parameters)
#hist(ips.ci.r.x.d$r.step.size[ips.ci.r.x.d$r.step.size > 0], xlab = "ips step size (r)", main = "", prob = TRUE)
hist(ips.ci.r.x.d$r.step.size, ylab = "", xlab = "", main = "", prob = TRUE, xlim = c(-0.6, 0.6), ylim = c(0, 4), breaks = 12, col = c(rep("grey30", 5), rep("grey70", 6)), las = 1)
curve(dexp(x, rate = ips.r.rate), from = 0, to = 0.55, col = 2, lty = 2, lwd = 2, add = TRUE)
mtext(text = "ips step size (r)", side = 1, line = 2)
#mtext(text = "count", side = 2, line = 2.5)

# mds step size histogram
mds.r.rate <- unname(eexp(mds.ci.r.x.d$r.step.size[mds.ci.r.x.d$r.step.size > 0])$parameters)
#hist(mds.ci.r.x.d$r.step.size[mds.ci.r.x.d$r.step.size > 0], xlab = "mds step size (r)", main = "", prob = TRUE)
hist(mds.ci.r.x.d$r.step.size, ylab = "", xlab = "", main = "", prob = TRUE, breaks = 10, xlim = c(-0.4, 0.5), col = c(rep("grey30", 4), rep("grey70", 5)), las = 1)
curve(dexp(x, rate = mds.r.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.45)
mtext(text = "mds step size (r)", side = 1, line = 2)
#mtext(text = "count", side = 2, line = 2.5)

# ptt step size histogram
ptt.r.rate <- unname(eexp(ptt.ci.r.x.d$r.step.size[ptt.ci.r.x.d$r.step.size > 0])$parameters)
#hist(ptt.ci.r.x.d$r.step.size[ptt.ci.r.x.d$r.step.size > 0], xlab = "ptt step size (r)", main = "", prob = TRUE)
hist(ptt.ci.r.x.d$r.step.size, ylab = "", xlab = "", main = "", prob = TRUE, breaks = 12, col = c(rep("grey30", 4), rep("grey70", 8)), las = 1)
curve(dexp(x, rate = ptt.r.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.75)
mtext(text = "ptt step size (r)", side = 1, line = 2.5)
#mtext(text = "count", side = 2, line = 2.5)

# mv step size histogram
mv.r.rate <- unname(eexp(mv.ci.r.x.d$mv.r.step.size[mv.ci.r.x.d$mv.r.step.size > 0])$parameters)
hist(mv.ci.r.x.d$mv.r.step.size, xlab = "", ylab = "", main = "", las = 1, breaks = 10, col = c(rep("grey30", 5), rep("grey70", 12)), ylim = c(0, 9)) 
curve(dexp(x, rate = mv.r.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.55)
mtext(text = "mv step size (r)", side = 1, line = 2.5)
mtext(text = "count", side = 2, line = 1.75, outer=TRUE)

# * 3.2 FIGURE 2: Fitting standardized step size (x) against exponential dist ----
par(mfrow = c(4,1),mar = c(4,1,1,0), oma = c(1, 3, 0, 1))
# ips standardized step size histogram
ips.x.rate <- unname(eexp(ips.ci.r.x.d$x.step.size[ips.ci.r.x.d$x.step.size > 0])$parameters)
hist(ips.ci.r.x.d$x.step.size, ylab = "", xlab = "", main = "", breaks = 12, col = c(rep("grey30", 8), rep("grey70", 4)), las = 1)
curve(dexp(x, rate = ips.x.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.35)
mtext(text = "ips standardized step size (x)", side = 1, line = 2.5)
#mtext(text = "count", side = 2, line = 2.5)

# mds standardized step size histogram
mds.x.rate <- unname(eexp(mds.ci.r.x.d$x.step.size[mds.ci.r.x.d$x.step.size > 0])$parameters)
hist(mds.ci.r.x.d$x.step.size, xlab = "", ylab = "", main = "", las = 1, breaks = 8, col = c(rep("grey30", 4), rep("grey70", 4)))
curve(dexp(x, rate = mds.x.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.35)
mtext(text = "mds standardized step size (x)", side = 1, line = 2.5)
#mtext(text = "count", side = 2, line = 2.5)

#ptt standardized step size histogram
ptt.x.rate <- unname(eexp(ptt.ci.r.x.d$x.step.size[ptt.ci.r.x.d$x.step.size > 0])$parameters)
hist(ptt.ci.r.x.d$x.step.size, xlab = "", ylab = "", main = "", las = 1, breaks = 10, col = c(rep("grey30", 4), rep("grey70", 6)), ylim = c(0, 4))
curve(dexp(x, rate = ptt.x.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.45)
mtext(text = "ptt standardized step size (x)", side = 1, line = 2.5)
#mtext(text = "count", side = 2, line = 2.5)

# mv standardized step size histogram
mv.x.rate <- unname(eexp(mv.ci.r.x.d$mv.x.step.size[mv.ci.r.x.d$mv.x.step.size > 0])$parameters)
hist(mv.ci.r.x.d$mv.x.step.size, xlab = "", ylab = "", main = "", las = 1, breaks = 10, col = c(rep("grey30", 9), rep("grey70", 12)), ylim = c(0, 12)) 
curve(dexp(x, rate = mv.x.rate), col = 2, lty = 2, lwd = 2, add = TRUE, from = 0.05, to = 0.55)
mtext(text = "mv standardized step size (x)", side = 1, line = 2.5)
mtext(text = "count", side = 2, outer=TRUE, line = 1.5,)

# end Section 3
#==============================================================================#

