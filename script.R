library(ggplot2)
library(gridExtra)
source("SD.R")
#source("new.R")


davSD <- hejSD(100, exList, theta = seq(0.1, 5, 0.1))

t1 <- Sys.time()
davnew <- hejnew(10, wardList, theta = seq(0.1, 10, 0.1))
t1 <- Sys.time() - t1

estdavSD <- lapply(davSD, function(x) x$df$est)
totestdavSD <- Reduce("+", estdavSD)
plot(seq(0.1, 5, by = 0.1), totestdavSD)

estdavnew <- lapply(davnew, function(x) x$df$est)
totestdavnew <- Reduce("+", estdavnew)
plot(seq(0.1, 10, by = 0.1), totestdavnew*10^24)

library("profvis")
profvis(likelihoodSD(10, ward, theta = seq(0.1, 10, by = 0.1)))










t1 <- Sys.time()
exLikeSD <- likelihoodSD(30000, exList[[8]], theta = seq(0.1, 10, by = 0.1))
t1 <- Sys.time() - t1
t2 <- Sys.time()
exLikenew <- likelihood_new(30000, exList[[8]], theta = seq(0.1, 10, by = 0.1))
t2 <- Sys.time() - t2

# save(exLikeSD, file = "exLikeSD30.Rdata")
# save(exLikenew, file = "exLikenew30.Rdata")

# load("exLikeSD30.Rdata")
# load("exLikenew30.Rdata")
# 
# load("exLikenew10new.Rdata")
# load("exLikeSD10SD.Rdata")

p1 <- ggplot(exLikeSD$df, aes(x = theta)) + 
  geom_line(aes(y = est*10^14)) + 
  geom_line(aes(y = q975*10^14), linetype = "dashed") + 
  geom_line(aes(y = q025*10^14), linetype = "dashed") + 
  scale_x_continuous(limits = c(0, 5)) +
  scale_y_continuous(name = "Likelihood(*10^11)", limits = c(0, 3)) +
  xlab(latex2exp::TeX("$\\theta$")) + 
  theme_bw()

p2 <- ggplot(exLikenew$df, aes(x = theta)) + 
  geom_line(aes(y = est*10^14)) + 
  geom_line(aes(y = q975*10^14), linetype = "dashed") + 
  geom_line(aes(y = q025*10^14), linetype = "dashed") + 
  scale_x_continuous(limits = c(0, 5)) +
  scale_y_continuous(name = "Likelihood(*10^11)", limits = c(0, 3)) +
  xlab(latex2exp::TeX("$\\theta$")) +
  theme_bw()
  
p3 <- ggplot(exLikeSD$df, aes(x = theta)) +
  geom_line(aes(y = log(est))) +
  xlim(c(0.5, 3.5)) +
  ylim(c(-29, -25))

p4 <- ggplot(exLikeSD$df, aes(x = theta)) +
  geom_line(aes(y = log(est))) +
  xlim(c(0.5, 3.5)) +
  ylim(c(-29, -25))


pOut <- grid.arrange(p1, p2, ncol = 2)
# 
# pdf("SDnewComp.pdf", width = 6, height = 3)
# plot(pOut)
# dev.off()

data.frame(the = exLikenew$df$theta, estNew = exLikenew$df$est, seNew = exLikenew$df$est - exLikenew$df$q025, estSD = exLikeSD$df$est, seSD = exLikeSD$df$est - exLikeSD$df$q025)





t1 <- Sys.time()
exLikeSDward <- likelihoodSD(100, wardList[[1]], theta = seq(0.1, 10, by = 0.1))
t1 <- Sys.time() - t1
t2 <- Sys.time()
exLikenewward <- likelihood_new(100, wardList[[1]], theta = seq(0.1, 10, by = 0.1))
t2 <- Sys.time() - t2

p1 <- ggplot(exLikeSDward$df, aes(x = theta)) + 
  geom_line(aes(y = est)) + 
  geom_line(aes(y = q975), linetype = "dashed") + 
  geom_line(aes(y = q025), linetype = "dashed")

p2 <- ggplot(exLikenewward$df, aes(x = theta)) + 
  geom_line(aes(y = est)) + 
  geom_line(aes(y = q975), linetype = "dashed") + 
  geom_line(aes(y = q025), linetype = "dashed")

# p3 <- ggplot(exLikeSD$df, aes(x = theta)) +
#   geom_line(aes(y = log(est))) +
#   xlim(c(0.5, 3.5)) +
#   ylim(c(-29, -25))
# 
# p4 <- ggplot(exLikeSD$df, aes(x = theta)) +
#   geom_line(aes(y = log(est))) +
#   xlim(c(0.5, 3.5)) +
#   ylim(c(-29, -25))

#grid.arrange(p1, p2, p3, p4, ncol = 2)

grid.arrange(p1, p2, ncol = 2)


