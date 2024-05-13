# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4626430/
#https://ebrary.net/192839/health/calculation_icers
####################################################################
####################################################################
# Name: getFrontier.R
# Goal: Find the CEA frontier, up to a given WTP level, by
#       identifying strategies with the highest NMB
# Written by: Sze Suen
# Date: Feb 25, 2015
# Notes:
#    ~Frontier strategies are displaced on the R output screen and
#      plotted in red on the scatter plot.
#
#    ~User needs to provide a csv file of costs and QALYs
#	    (CostQalyInputFile_online_supp.csv) inside the folder specified
#	    below (inputFolder). The CSV should have three columns (labeled
#     in first row) in this order:
#      Strategy number, costs, and QALYs.
#
#    ~User can specify the maximum willingness-to-pay level to
#      consider (maxWTP).  Can be Inf for infinity.
#
#    ~QALY-reducing strategies will be on the frontier if they save
#      enough money (script assumes maximum willingness to save money
#      per QALY lost is equivalent to maximum willingness to pay per QALY
#      gained). If the user does not wish to consider such policies as
#      being on the frontier, do not include strategies with negative
#      QALYs in the input csv file.
#
#    ~Script does not use the one-line code cited in the text
#      as the max function is slow. This implementation is
#      faster and methodologically does the same thing.
#
#    ~May take a few minutes if thousands of strategies and
#       processing resources are low.  Please be patient.
#
#    Please cite article if this code is used.
#
# USER INPUTS:
#inputFolder <- "~/Downloads/"
maxWTP <- 200000        # any positive value or Inf

####################################################################
####################################################################

# load the costs and QALY values for all strategies
#scriptFolder <- getwd();  setwd(inputFolder)
#CEmat <- as.matrix(read.csv("~/Downloads/NIHMS677243-supplement-sup1.csv", header = TRUE))

library(tidyverse)
CEmat <- data.frame(Strategy = c("A","B","C","D","E"),
                Cost = c(16457, 21458, 32877, 24504, 43331),
                QALYs = c(17.335, 17.41, 17.579, 17.491, 17.491)

)  %>%
  column_to_rownames(var = "Strategy") %>%
  as.matrix()


# check for duplicated strategies
dups <- CEmat[c(duplicated(CEmat[,2:3]) | duplicated(CEmat[,2:3], fromLast = TRUE)),1]

# initialize some variables
costsCol <- 1; qalyCol <- 2
numStrat <- dim(CEmat)[1]

# find WTP levels to test so that all strategies on frontier will be captured
# this means testing on either side of all NMB intersections, which are just all the pairwise ICERs
ICERmat <- matrix(1, numStrat, numStrat)
for (i in 1:numStrat ) {
  indexStrat <- matrix(1, numStrat, 2)
  indexStrat[,costsCol] <- indexStrat[,costsCol]*CEmat[i,costsCol]
  indexStrat[,qalyCol] <- indexStrat[,qalyCol]*CEmat[i,qalyCol]
  delCostQalys <- CEmat - indexStrat
  ICERmat[,i] <- delCostQalys[,costsCol] / delCostQalys[,qalyCol]
}
intersections <- sort(unique(c(ICERmat)))
intersections <- intersections[is.finite(intersections)]
WTPtestPoints <- c(0, intersections [intersections >= 0], maxWTP)

# Find the strategy with the max NMB at each of the WTP test points
indiciesOfMax <- vector()
NMBmat <- matrix(0, numStrat, length(WTPtestPoints))
for (i in 1:length(WTPtestPoints) ) {
  NMBmat[,i] <- (WTPtestPoints[i]*CEmat[,qalyCol]) - CEmat[,costsCol]
}
if (is.infinite(maxWTP)) {
  #WTP of infinity means costs are not considered
  NMBmat[,length(WTPtestPoints)] = CEmat[,qalyCol] - (0*CEmat[,costsCol]);
}
maxVals <- apply(NMBmat, 2, max)  #find strategy that maximizes NMB at each WTP
for (i in 1:length(WTPtestPoints) ) {  #find all strategies that match max at each WTP
  indiciesOfMax <- c(indiciesOfMax,which( NMBmat[,i] == maxVals[i]))
}
frontier <- unique(indiciesOfMax)  #find strategy that maximizes NMB at each WTP

# display out: make plot and print to output screen


plot(CEmat[,qalyCol], CEmat[,costsCol])
points(CEmat[frontier,qalyCol], CEmat[frontier,costsCol], col = 'red', pch = 16)

if (length(dups)>0){
  warning("Strategies have the same costs and benefits (displayed above)")
  print(dups)
}
sprintf("Frontier is formed by strategies: %s", paste( sort(CEmat[frontier,1]), collapse=" "))


