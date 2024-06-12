#################### Packages to be used ####################
# install.packages("SuperLearner")      
# install.packages("foreach")     
# install.packages("doParallel") 
############################################################# 

library(foreach)
library(doParallel)
library(SuperLearner)
library(parallel)


###### Sets up the parallel programming ###### 
cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #Use all but 1 core
registerDoParallel(cl)
##############################################
clusterEvalQ(cl, source("SD.R"))

t1 <- Sys.time()
sampnew1 <- parLapply(cl, wardList, function(x) likelihood_new(N = 300, data = x, theta = seq(0.1, 10, by = 0.1)))
t1 <- Sys.time() - t1

t1

t2 <- Sys.time()
sampSD2 <- parLapply(cl, wardList, function(x) likelihoodSD(N = 300, data = x, theta = seq(0.1, 10, by = 0.1)))
t2 <- Sys.time() - t2

SDcur <- prepImps(combSD1, theta = seq(0.1, 10, by = 0.1))
SDp <- plotImps(SDcur)

newcur <- prepImps(combnew1, theta = seq(0.1, 10, by = 0.1))
newp <- plotImps(newcur)

combnew1 <- bindLists(sampnew1, sampnew, theta = seq(0.1,10,by=0.1))
combSD1 <- bindLists(sampSD2, sampSD, theta = seq(0.1,10,by=0.1))

save(combnew1, file = "combnew.Rdata")
save(combSD1, file = "combSD.Rdata")

load("combnew.Rdata")
load("combSD.Rdata")
hej <- prepImps(combSD1, theta = seq(0.1,10,by=0.1))
plotImps(hej)

save(sampSD, file = "sampSD.Rdata")
save(sampnew, file = "sampnew.Rdata")


bla <- load("sampnewtest.Rdata")
load("sampSDtest.Rdata")






########## Ex data ##########
t1 <- Sys.time()
newTest <- likelihood_new_list(exList, theta = seq(0.1, 6, by = 0.2))
t1 <- Sys.time() - t1

save(newTest, file = "newTest.Rdata")
newTest <- NULL

t2 <- Sys.time()
SDTest <- likelihoodSD_list(exList, theta = seq(0.1, 6, by = 0.2))
t2 <- Sys.time() - t2

save(SDTest, file = "SDTest.Rdata")
SDTest <- NULL

# t1 <- Sys.time()
# sampnew2 <- parLapply(cl, exList, function(x) likelihood_new(N = 10, data = x, theta = seq(0.1, 10, by = 0.1)))
# t1 <- Sys.time() - t1
# 
# 
# #t2 <- Sys.time()
# sampSD2 <- parLapply(cl, exList, function(x) likelihoodSD(N = 30000, data = x, theta = seq(0.1, 10, by = 0.1)))
# #t2 <- Sys.time() - t2
# t1 <- Sys.time() - t1
# 
# SDcur2 <- prepImps(sampSD2, theta = seq(0.1, 7.5, by = 0.1))
# SDp2 <- plotImps(SDcur2)
# 
# newcur2 <- prepImps(sampnew2, theta = seq(0.1, 7.5, by = 0.1))
# newp2 <- plotImps(newcur2)
