library(plyr)
library(dplyr)
library(foreach)
library(doParallel)
library(SuperLearner)
library(tibble)

sTrans <- function(a){
  sum(a * 10^((length(a) - 1):0))
}

genSingleRoot <- function(h){
  ord <- rev(order(apply(h[, -dim(h)[2]], 2, sTrans)))
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
  # seqs <- which(rowSums(h[, -dim(h)[2]]) > 1)
  # for(k in seqs){
  #   j <- which(h[k, -dim(h)[2]] == 1)
  #   for(l in j[-1]){
  #     tmp <- alls
  #     tmp[, l] <- abs(tmp[, l] - 1)
  #     roots_tmp2 <- c(roots_tmp2, list(cbind(tmp, n)))
  #   }
  # }
  # roots_tmp2 <- roots_tmp2[-1]
  # roots_tmp <- c(roots_tmp1, roots_tmp2)
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
    k_tmp3 <- sum(xkxj) == 0
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
    h <- Rk(h, k, j)
    scene <- 3
  }
  list(h = h, scene = scene, n = n, nk = nk, n0 = n0, s = s)
}


sk <- function(h, k){
  if(length(k) == 0)
    return(0)
  sapply(k, function(x) sum(h[x, ][- dim(h)[2]]))
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
    w <- (H$n0 * theta) / (H$s * H$n * (H$n - 1 + theta))
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

sderr <- function(x){
  sd(x)/sqrt(length(x))
}

likelihood_new <- function(N, data, theta){
  imps <- Imp_new(N, data, theta)
  est <- sapply(imps, mean)
  se <- sapply(imps, sderr)
  q025 <- est - se
  q975 <- est + se
  list(df = data.frame(theta = theta, est = est, q025 = q025, q975 = q975), imps = imps)
}

### new proposal ###
probQ_new <- function(h, crit1, crit2){
  skcrit1 <- sk(h, crit1$k)
  skcrit2 <- sk(h, crit2$k)
  n0 <- sum(h$n[h$n >= 2]) + sum(skcrit1) + sum(skcrit2)
  q <- h$n / n0
  q[crit1$k] <- skcrit1
  q[crit2$k] <- skcrit2
  set_zero <- which(h$n == 1)[!(which(h$n == 1) %in% crit1$x0dist)]
  q[set_zero] <- 0
  list(q = q, n0 = n0)
}