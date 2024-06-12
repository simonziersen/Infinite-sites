library(plyr)
library(dplyr)
library(foreach)
library(doParallel)
library(SuperLearner)
library(ggplot2)
library(gridExtra)
library(tibble)

sTrans <- function(a){
  sum(a * 10^((length(a) - 1):0))
}

bindImps <- function(l, theta){
  tmp <- list()
  for(i in 1:length(theta)){
    for(k in 1:length(l)){
      if(k == 1){
        tmp[[i]] <- l[[k]]$imps[[i]]
      }else{
        tmp[[i]] <- tmp[[i]] + l[[k]]$imps[[i]]
      }
    }
  }
  tmp
}

bindImps2 <- function(l, theta){
  tmp <- list()
  for(i in 1:length(theta)){
    for(k in 1:length(l)){
      if(k == 1){
        tmp[[i]] <- l[[k]][[i]]
      }else{
        tmp[[i]] <- tmp[[i]] + l[[k]][[i]]
      }
    }
  }
  tmp
}

bindLists <- function(l1, l2, theta){
  tmp1 <- l1
  for(i in 1:length(theta)){
    for(k in 1:length(l1)){
      tmp1[[k]]$imps[[i]] <- c(l1[[k]]$imps[[i]], l2[[k]]$imps[[i]])
    }
  }
  tmp1
}

prepImps <- function(samp, theta){
  imps <- bindImps(samp, theta)
  m <- sapply(imps, mean)
  se <- sapply(imps, sderr)
  q025 <- m - se
  q975 <- m + se
  data_frame(theta = theta, est = m, q025 = q025, q975 = q975)
}

plotImps <- function(df){
  p <- ggplot(df, aes(x = theta)) +
    geom_line(aes(y = est*10^20)) + 
    geom_line(aes(y = q975*10^20), linetype = "dashed") + 
    geom_line(aes(y = q025*10^20), linetype = "dashed")
  p
}

genSingleRoot <- function(h){
  ord <- rev(order(apply(h[, - dim(h)[2]], 2, sTrans)))
  h[, c(ord, dim(h)[2])]
}

genRoots <- function(h){
  del <- which(rowSums(h[, -dim(h)[2]]) == 0)
  if(length(del) == 0){
    cap <- 1:(dim(h)[1])
  }else{
    cap <- (1:(dim(h)[1]))[-del]
  }
  alls <- h[, - dim(h)[2]]
  n <- h[, dim(h)[2]]
  roots_tmp1 <- list()
  roots_tmp1[[1]] <- h
  roots_tmp2 <- NULL
  for(i in cap){
    tmp <- abs(t(apply(alls, 1, function(x) x - as.numeric(alls[i, ]))))
    roots_tmp1[[i+1]] <- as_tibble(cbind(tmp, n))
  }
  roots_tmp <- roots_tmp1
  if(length(del) > 0)
    roots_tmp <- roots_tmp[-(del+1)]
  roots <- lapply(roots_tmp, genSingleRoot)
  roots
}

ISdata <- function(data){
  l <- list()
  l[[1]] <- data
  l
}

sk <- function(h, k){
  if(length(k) == 0)
    return(0)
  sapply(k, function(x) sum(h[x, ][- dim(h)[2]]))
}

f <- function(seq){
  #tmp <- seq
  # tmp[max(which(tmp == 1))] <- 0
  # tmp
  seq[max(which(seq == 1))] <- 0
  seq
}
  

fk <- function(h, k){
  tmp <- h
  tmp[k, ][- dim(h)[2]] <- f(tmp[k, ][- dim(h)[2]])
  tmp
}

Rk <- function(h, k, j){
  tmp <- h
  tmp$n[j] <- tmp$n[j] + 1
  tmp1 <- tmp[-k, ]
  tmp1
}

sderr <- function(x){
  sd(x)/sqrt(length(x))
}

fxkxj <- function(h, k){
  h_tmp <- h[, - dim(h)[2]]
  # if(sum(h_tmp[k, ]) == 0)
  #   return(FALSE)
  tmp <- duplicated(rbind(h_tmp, f(h_tmp[k, ])), fromLast = T)
  tmp[- length(tmp)]
  #apply(h_tmp, 1, function(x) all(f(h_tmp[k, ]) == x))
}

Crit1 <- function(h){
  k_tmp <- which(h$n == 1)
  x0 <- apply(h[k_tmp, ][- dim(h)[2]], 1, 
              function(x){
                if(sum(x == 1) == 0){
                  0
                  }else{
                   max(which(x == 1))
                  }
                })
  if(0 %in% x0){
    k_tmp <- k_tmp[- which(x0 == 0)]
    x0 <- x0[- which(x0 == 0)]
  }
  tmp <- apply(h[, x0], 2, sum)
  k_tmp2 <- k_tmp[which(tmp == 1)]
  if(length(k_tmp2) == 0){
    k_tmp3 <- list()
  }else{
    xkxj <- sapply(k_tmp2, function(x) fxkxj(h, x))
    k_tmp3 <- colSums(xkxj) == 0
  }
  if(length(k_tmp3) == 0){
    return(list(x0dist = integer(0), k = integer(0)))
  } else{
    return(list(x0dist = k_tmp2, k = k_tmp2[k_tmp3], xkxj = xkxj))
  }
}

Crit2 <- function(h, Crit){
  if(length(Crit) == 1)
    return(0)
  if(length(Crit$x0dist) == 0)
    return(list(k = integer(0), j = list()))
  k <- Crit$x0dist[!(Crit$x0dist %in% Crit$k)]
  j <- which(Crit$xkxj, arr.ind = T)[, 1]
  list(k = k, j = j)
}

##### SD #####

sim_oneSD <- function(h){
  crit1 <- Crit1(h)
  crit2 <- Crit2(h, crit1)
  q_tmp <- probQSD(h, crit1, crit2)
  q <- q_tmp$q
  n0 <- q_tmp$n0
  k <- sample(1:dim(h)[1], 1, prob = q)
  n <- sum(h$n)
  nk <- NA
  nj <- NA
  k_help <- NA
  if(h$n[k] > 1){
    nk <- h$n[k]
    h$n[k] <- h$n[k] - 1
    scene <- 1
  }
  if(k %in% crit1$k){
    h <- fk(h, k)
    scene <- 2
  }
  if(k %in% crit2$k){
    j <- crit2$j[crit2$k == k]
    nj <- h$n[j]
    k_help <- k
    h <- Rk(h, k, j)
    scene <- 3
  }
  list(h = h, scene = scene, n = n, nk = nk, n0 = n0, nj = nj, k_help = k_help)
}


simSD <- function(h){
  H <- list()
  H[[1]] <- h
  H[[2]] <- sim_oneSD(h)
  for(i in 3:10000){
    if(dim(H[[i-1]]$h)[1] == 1){
      if(H[[i-1]]$h$n == 1)
        break
    }
    H[[i]] <- sim_oneSD(H[[i-1]]$h)
  }
  H
}

wSD <- function(H, theta){
  if(H$scene == 1){
    w <- (H$n0 * (H$nk - 1)) / (H$nk * (H$n - 1 + theta))
  }
  if(H$scene == 2){
    w <- (H$n0 * theta) / (H$n * (H$n - 1 + theta))
  }
  if(H$scene == 3){
    w <- (H$n0 * theta * (H$nj + 1)) / (H$n * (H$n - 1 + theta))
  }
  w  
}

weightsSD <- function(H, theta){
  ws <- sapply(H[-1], function(x) wSD(x, theta))
  W <- prod(ws)
  W
}

ImpSD <- function(N, data, theta){
  W <- list()
  sims <- NULL
  for(i in 1:N){
    sims[[i]] <- simSD(data)
  }
  for(i in seq_along(theta)){
    W[[i]] <- sapply(sims, function(x) weightsSD(x, theta[i]))
  }
  W
}

likelihoodSD <- function(N, data, theta){
  imps <- ImpSD(N, data, theta)
  est <- sapply(imps, mean)
  se <- sapply(imps, sderr)
  q025 <- est - se
  q975 <- est + se
  list(df =data.frame(theta = theta, est = est, q025 = q025, q975 = q975), imps = imps)
}

likelihoodSD_list <- function(data, theta){
  tmp <- list()
  tmp <- parLapply(cl, data, function(x) ImpSD(N = 30000, x, theta = seq(0.1, 5, by = 0.2)))
  wSum <- bindImps2(tmp, theta)
  est <- sapply(wSum, mean)
  se <- sapply(wSum, sderr)
  q025se <- est - se
  q975se <- est + se
  #q025 <- sapply(wSum, function(x) quantile(x, 0.025))
  #q975 <- sapply(wSum, function(x) quantile(x, 0.975))
  tmp <- NULL
  list(df = data.frame(theta = theta, est = est, q025 = q025se, q975 = q975se))
}

probQSD <- function(h, crit1, crit2){
 n0 <- sum(h$n[h$n >= 2]) + sum(length(crit1$k)) + sum(length(crit2$k))
 q <- h$n / n0
 set_zero <- which(h$n == 1)[!(which(h$n == 1) %in% crit1$x0dist)]
 q[set_zero] <- 0
 list(q = q, n0 = n0)
}


##### new #####
sim_one_new <- function(h){
  crit1 <- Crit1(h)
  crit2 <- Crit2(h, crit1)
  q_tmp <- probQ_new(h, crit1, crit2)
  q <- q_tmp$q
  n0 <- q_tmp$n0
  k <- sample(1:dim(h)[1], 1, prob = q)
  n <- sum(h$n)
  s <- sk(h, k)
  nk <- NA
  nj <- NA
  if(h$n[k] > 1){
    nk <- h$n[k]
    h$n[k] <- h$n[k] - 1
    scene <- 1
  }
  if(k %in% crit1$k){
    h <- fk(h, k)
    scene <- 2
  }
  if(k %in% crit2$k){
    j <- crit2$j[crit2$k == k]
    nj <- h$n[j]
    h <- Rk(h, k, j)
    scene <- 3
  }
  list(h = h, scene = scene, n = n, nk = nk, n0 = n0, s = s, nj = nj, qsum = sum(q))
}

sim_new <- function(h){
  H <- list()
  H[[1]] <- h
  H[[2]] <- sim_one_new(h)
  for(i in 3:10000){
    if(dim(H[[i-1]]$h)[1] == 1){
      if(H[[i-1]]$h$n == 1)
        break
    }
    H[[i]] <- sim_one_new(H[[i-1]]$h)
  }
  H
}

w_new <- function(H, theta){
  if(H$scene == 1){
    w <- (H$n0 * (H$nk - 1)) / (H$nk * (H$n - 1 + theta))
  }
  if(H$scene == 2){
    w <- (H$n0 * theta) / (H$s * H$n * (H$n - 1 + theta))
  }
  if(H$scene == 3){
    w <- (H$n0 * theta * (H$nj + 1)) / (H$s * H$n * (H$n - 1 + theta))
  }
  w  
}

weights_new <- function(H, theta){
  ws <- sapply(H[-1], function(x) w_new(x, theta))
  W <- prod(ws)
  W
}

Imp_new <- function(N, data, theta){
  W <- list()
  sims <- NULL
  for(i in 1:N){
    sims[[i]] <- sim_new(data)
  }
  for(i in seq_along(theta)){
    W[[i]] <- sapply(sims, function(x) weights_new(x, theta[i]))
  }
  W
}

likelihood_new_list <- function(data, theta){
  tmp <- list()
  tmp <- parLapply(cl, data, function(x) Imp_new(N = 30000, x, theta = seq(0.1, 6, by = 0.2)))
  wSum <- bindImps2(tmp, theta)
  est <- sapply(wSum, mean)
  se <- sapply(wSum, sderr)
  q025se <- est - se
  q975se <- est + se
  #q025 <- sapply(wSum, function(x) quantile(x, 0.025))
  #q975 <- sapply(wSum, function(x) quantile(x, 0.975))
  tmp <- NULL
  list(df = data.frame(theta = theta, est = est, q025 = q025se, q975 = q975se))
}


likelihood_new <- function(N, data, theta){
  imps <- Imp_new(N, data, theta)
  est <- sapply(imps, mean)
  se <- sapply(imps, sderr)
  q025 <- est - se
  q975 <- est + se
  list(df = data.frame(theta = theta, est = est, q025 = q025, q975 = q975), imps = imps)
}

probQ_new <- function(h, crit1, crit2){
  skcrit1 <- sk(h, crit1$k)
  skcrit2 <- sk(h, crit2$k)
  n0 <- sum(h$n[h$n >= 2]) + sum(skcrit1) + sum(skcrit2)
  q <- h$n / n0
  q[crit1$k] <- skcrit1 / n0
  q[crit2$k] <- skcrit2 / n0
  set_zero <- which(h$n == 1)[!(which(h$n == 1) %in% crit1$x0dist)]
  q[set_zero] <- 0
  list(q = q, n0 = n0)
}

