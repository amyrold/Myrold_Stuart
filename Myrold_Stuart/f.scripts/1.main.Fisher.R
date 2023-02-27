# General ----
# Working Title: Stickles back Fisher OR Stickleback fish back Fisher
# Authors: Yoel E. Stuart, Aaron Myrold, and Michael A. Bell

# Questions and Motivation ----
#This paper is motivated by Fisher's Geometric Model of Adaptation (1930), discussed in Orr (2005, Nature Reviews Genetics) and citations therein. The original model and subsequent modifications make two main predictions: 
#1) The distribution of effect sizes of mutations that are substituted during an adaptive walk should be approximately exponential. That is, adaptation in Fisher's model should involve a few mutations of relatively large phenotypic effect and many of relatively small effect. 
#2) The mean effect sizes of mutations substituted at the first versus the second substitution, and so on for subsequent substitutions, should fall off by a constant proportion, approximating a geometric sequence. That is, larger-effect substitutions should be substituted early on in the adaptive walk, with smaller effects substituted later. 

#An interesting problem we'll need to investigate, understand, and explain in any paper is whether it's fair to equate time-averaged phenotypic steps toward the adaptive peak (i.e., what we can observe from the fossils) with steps toward the adaptive peak caused by mutation substitutions. Fisher's theorem is about substitution effect sizes. We will be reporting phenotypic change, at 250yr intervals. Is it likely that multiple substitutions are fixed during those intervals? If yes, what can we say about Fisher's theorem? Or, is it fair to assume that beneficial mutations are rare and that we'd only expect one adaptive mutation to fix in 250 years? Can this be parameterized by existing stickleback data?

#Another interesting problem we'll need to investigate, understand, and explain, is whether and how Fisher's model treats standing genetic variation in the direction of the peak. We know that populations can partially adapt over short time scales because they have genetic variation pre-adapted to the environmental change. What is the theory underlying the effect sizes of substitution of standing variation, and what is the prediction for phenotypic steps during adaptation from standing genetic variation? For example, Barrett and Schluter (2008) suggest that adaptation from standing genetic variations should fix more substitutions of small effect than adaptation from de novo mutation. Then, How does this distribution overlap with the Fisher distribution, and what does their joint distribution look like, if adaptation proceeds through both? (Does the distribution for pelvic score, which proceeded from a de novo mutation look different for the other traits that had standing variation?)

#We have the equivalent of the mutational effect size (given some assumptions) represented by the phenotypic steps (standardized values x against the distance to origin) taken by the population as it evolves toward the new optimum. Again, given some assumptions, we would expect these step sizes to be exponentially distributed, with most small  and a few large. This is derived in Orr 1998, Conclusions 3 and 7.

#Can we detect a joint distribution of standing genetic variation and de novo mutation by partitioning the step size distributions into sections that match the exponential and sections that don't. We might also be able to provide an empirical joint distribution with proposed splits between the de novo and standing, that theorists can use to parameterize their models.

# R version
R.version.string # "R version 4.0.2 (2020-06-22)"
rm(list = ls())

#==============================================================================#

# Script Index ----
# 1.main.R        
# 2.data.manip.R
# 3.analysis.R
# 4.functions.R

#==============================================================================#

# To do ----
# 1) Better document where the raw data are coming from. I think that I got these from Matt Travis.
#2) We want to know if theorists have worked out the distribution of effect sizes of adaptation from standing variation (literature review and writing)
#3) We want to know if theorists have worked out the joint distribution of effect sizes of adaptation from standing variation AND de novo (literature review and writing)
#4) We want to know what empirical studies have revealed. How rare is our question? Microbial and viral evolution studies (Lenski Long Term Evolution Experiment). Peter and Rosemary grant adaptation by Darwin's finches. (literature review and writing)
#5) Figure out what to do with movements away from the optimum. (i.e., what happens when x values are negative).

#Proximate to dos -----
#remove outliers
#add error about the mv distance to optimum in Figure 1, top right panel

#==============================================================================#

# Global Variables ----
# store current user's working directory
wk.dir <- getwd() #/Users/ystuart/Dropbox/Loyola/Projects/Myrold_Fisher"
# plotting parameters
options(max.print = 10000)
clrs <- c("#000000", "#FF0000") # c("black", "red")
# for data manipulation and outlier checking
name.new.fishID <- "fishID.univ" # name used for the universal fishID, used across all projects
outlier.thresh <- 3.5 # how many times the sd a value needs to be to be tagged as an outlier


#==============================================================================#

# Folder Management ----
# names of folders for output data (figures + data output)
folder.names <- c("a.data.raw","b.data.clean", "c.results","d.figures", 'e.manuscripts', 'f.scripts', 'g.trash')
# if folder with name "i" does not exist, create it.
for(i in 1:length(folder.names)){
  if(file.exists(folder.names[i]) == FALSE){
    dir.create(folder.names[i]) 
  }
}

# paths to the folders. The 'p.' indicates the variable is a path.
# make sure the variable names describe the folder.names
p.data.raw <- paste(wk.dir, "/", folder.names[1], "/", sep = "")
p.data.clean <- paste(wk.dir, "/", folder.names[2], "/", sep = "")
p.results <- paste(wk.dir, "/", folder.names[3], "/", sep = "")
p.fig <- paste(wk.dir, "/", folder.names[4], "/", sep = "")

#==============================================================================#

# Libraries ----
# load libraries needed for our analyses
#library(gt)
library(gamlss) #distribution fitting
library(gamlss.dist)
library(gamlss.add)
library(EnvStats)
library(abind)
library('KScorrect')

#==============================================================================#

# Run Scripts ----
source("4.functions.R")
#source("2.data.manip.R")
#source("3.analysis.Fisher.R")
#source("5.figures.R")

# END OF MAIN ----
#==============================================================================#

