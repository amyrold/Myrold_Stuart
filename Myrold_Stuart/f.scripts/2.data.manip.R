# Table of Contents ----
# 1. Import Data and Housekeeping
# * 1.1 Import
# * 1.2 Housekeeping
# * 1.3 Optimum Definition. From Hunt et al. 2008
# * 1.4 Length of adaptive walk in years Definition. From Hunt et al. 2008.
# end Table of Contents

#==============================================================================#

# 1. Import Data and Housekeeping ----
# * 1.1 Import Data ----
#individual level l.series with dorsal spines (mds), pelvic score (ips), and touching pterygiophores (ptt)
#These data are what were analyzed in Bell et al. 2006 and Hunt et al. 2008. 
data.l.sample.sbf <- read.csv(paste(p.data.clean, "l.series_from_FossilIntervalData_recreatesAppendix1Bell2006.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)

#==============================================================================#

# * 1.2 Housekeeping ----
# To Do: outlier check and removal ----
#==============================================================================#

# * 1.3 Multivariate adaptive optimum - global variable ----
#based on the theta estimates from Hunt et al. 2008
#optimum number of dorsal spines
#ln(x + 1) 
theta.mds.transformed <- 0.79
theta.mds <- exp(0.79) - 1

#optimum number of touching pterygiophores
#ln(x + 1)
theta.ptt.transformed <- 0.82
# NB: inferred optimum from Hunt et al. 2008 is below all the data but one point in the data set we are using.
# NB: sample means do not re-create Appendix 1 from Bell et al. 2006, which was used directly for analyses in Hunt et al. 2008
theta.ptt <- exp(0.82) - 1 + 1 
#the plus one here gets optimum in line with actual data. We think that Bell et al. 2006 might have run log(x + 1) for mds and ips but merely log(x) for ptt. No one from Hunt et al. remembers what was done.

#optimum pelvic score
#ln(x + 1)
theta.ips.transformed <- 0.73
theta.ips <- exp(0.73) - 1

optimum <- c(theta.ips, theta.mds, theta.ptt)

# end Data Manip 

#==============================================================================#

