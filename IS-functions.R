library(dplyr)



#### generate all rooted trees from one unrooted tree ####
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
  if(sum(duplicated(t(h[, -dim(h)[2]]))) > 0){
    s <- which(duplicated(t(h[, -dim(h)[2]])))
    for(i in s){
      tmp <- h
      tmp[, i] <- abs(tmp[, i] - 1)
      roots_tmp2[[which(s == i)]] <- genSingleRoot(tmp)
    }
  }
  roots <- c(roots, roots_tmp2)
  roots
}

genSingleRoot <- function(h){
  ord <- rev(order(apply(h[, - dim(h)[2]], 2, sTrans)))
  h[, c(ord, dim(h)[2])]
}

sTrans <- function(a){
  sum(a * 10^((length(a) - 1):0))
}

#### Generator for rooted trees ####
RT <- function(tree){
  rt <- genSingleRoot(tree)
  c1 <- crit1(rt)
  c2 <- crit2(rt, c1)
  rt$c <- rep(NA, dim(rt)[1])
  rt$c[rt$n > 1] <- "k1"
  rt$c[c1$k] <- "k2"
  rt$c[c2$k] <- "k3"
  rt$s <- rowSums(rt[, -((dim(rt)[2] - 1):dim(rt)[2])])
  n0 <- sum(rt$n[rt$n > 1]) + length(c1$k) + length(c2$k)
  structure(list(rt = rt, n0 = n0), class = "rootedTree")
}

crit1 <- function(h){
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

crit2 <- function(h, Crit){
  if(length(Crit) == 1)
    return(0)
  if(length(Crit$x0dist) == 0)
    return(list(k = integer(0), j = list()))
  k <- Crit$x0dist[!(Crit$x0dist %in% Crit$k)]
  j <- which(Crit$xkxj, arr.ind = T)[, 1]
  list(k = k, j = j)
}

fxkxj <- function(h, k){
  h_tmp <- h[, - dim(h)[2]]
  tmp <- duplicated(rbind(h_tmp, f(h_tmp[k, ])), fromLast = T)
  tmp[- length(tmp)]
}

f <- function(seq){
  seq[max(which(seq == 1))] <- 0
  seq
}

#### Simulating function ####
sim <- function(object, method = "SD"){
  class(object)[2] <- method
  H <- list()
  H[[1]] <- object
  H[[2]] <- sim_one(object)
  for(i in 3:10000){
    if(dim(H[[i-1]]$rt)[1] == 1){
      if(H[[i-1]]$rt$n == 1)
        break
    }
    H[[i]] <- sim_one(H[[i-1]])
  }
  H
}

sim_one <- function(object)
  UseMethod("sim_one")

sim_one.SD <- function(object){
  q <- probQ(object)
  rt <- object$rt
  k <- sample(dim(rt)[1] , 1, prob = q)
  if(rt$c[k] == "k1"){
    object$rt$n[k] <- object$rt$n[k] - 1
    object$scene <- 1
    object$n0 <- object$n0 - 1
    if(rt$n[k] == 2){
      if(rt$s[k] == 0){
        object$rt$c[k] <- NA
        object$n0 <- object$n0 - 1
      } else{
        xk <- which(rt[k, 1:(dim(rt)[2] - 3)] == 1)
        xk0 <- max(xk)
        if(sum(rt[, xk0]) > 1){
          object$rt$c[k] <- NA
          object$n0 <- object$n0 - 1
        } else{
          if(rt$s[k ]== 1){
            if(any(rt$s == 0)){
              object$rt$c[k] <- "k3"
            } else{
              object$rt$c[k] <- "k2"
            }
          } else{
            xk1 <- max(setdiff(xk, xk0))
            dups <- which(rt[, xk1] == 1)
            if(length(dups) > 1){
              for(i in setdiff(dups, k)){
                k_tmp <- FALSE
                if(all(f(rt[k, -((dim(rt)[2] - 2):dim(rt)[2])]) == rt[i, -((dim(rt)[2] - 2):dim(rt)[2])])){
                  k_tmp = TRUE
                  break
                }
              }
              if(k_tmp){
                object$rt$c[k] <- object$rt$c[k] <- "k3"
              } else{
                object$rt$c[k] <- object$rt$c[k] <- "k2" 
              }
            } else{
              object$rt$c[k] <- "k2"                         
            }
          }
        }
      }
    } 
  }
  if(rt$c[k] == "k2"){
    object$rt$s[k] <- object$rt$s[k] - 1
    object$rt <- fk(object$rt, k)
    object$scene <- 2
    scene <- 2
    if(rt$s[k] > 1){
      xk <- which(rt[k, 1:(dim(rt)[2] - 3)] == 1)
      xk1 <- rev(xk)[2]
      if(sum(rt[, xk1]) > 1){
        object$rt$c[k] <- NA
        object$n0 <- object$n0 - 1
      } else{
        tmp <- object$rt[, 1:(dim(object$rt)[2] - 3)]
        tmp2 <- anyDuplicated(rbind(tmp[-k, ], tmp[k, ]), fromLast = T)
        if(tmp2 > 0){
          object$rt$c[k] <- "k3"
        } else{
          object$rt$c[k] <- "k2"
        }
      }
    } else{
      object$rt$c[k] <- NA
      object$n0 <- object$n0 - 1
    }
  }
  if(rt$c[k] == "k3"){
    object$scene <- 3
    tmp <- rt[, 1:(dim(rt)[2] - 3)]
    k_tmp <- f(tmp[k, ])
    j <- anyDuplicated(rbind(tmp, k_tmp), fromLast = T)
    object$rt$c[j] <- "k1"
    object$rt <- Rk(object$rt, k, j)
    if(is.na(rt$c[j])){
      object$n0 <- object$n0 + 1
    } 
  }
  object
}

probQ <- function(object)
  UseMethod("probQ")

probQ.SD <- function(object){
  h <- object$rt
  q <- h$n / object$n0
  q[is.na(h$c)] <- 0
  q
}

fk <- function(h, k){
  tmp <- h
  tmp[k, ][- ((dim(h)[2] - 2):dim(h)[2])] <- f(tmp[k, ][- ((dim(h)[2] - 2):dim(h)[2])])
  tmp
}

Rk <- function(h, k, j){
  tmp <- h
  tmp$n[j] <- tmp$n[j] + 1
  tmp1 <- tmp[-k, ]
  tmp1
}




#### For testing - not yet modified ####
w <- function(object, theta){
  if(object$scene == 1){
    w <- (object$n0 * (object$nk - 1)) / (object$nk * (object$n - 1 + theta))
  }
  if(object$scene == 2){
    w <- (object$n0 * theta) / (object$n * (object$n - 1 + theta))
  }
  if(H$scene == 3){
    w <- (H$n0 * theta * (H$nj + 1)) / (H$n * (H$n - 1 + theta))
  }
  w  
}

weights <- function(H, theta){
  ws <- sapply(H[-1], function(x) wSD(x, theta))
  W <- prod(ws)
  W
}

Imp <- function(N, object, theta){
  W <- list()
  sims <- NULL
  for(i in 1:N){
    sims[[i]] <- sim(object)
  }
  for(i in seq_along(theta)){
    W[[i]] <- sapply(sims, function(x) weightsSD(x, theta[i]))
  }
  W
}

likelihood <- function(N, data, theta){
  imps <- ImpSD(N, data, theta)
  est <- sapply(imps, mean)
  se <- sapply(imps, sderr)
  q025 <- est - se
  q975 <- est + se
  list(df =data.frame(theta = theta, est = est, q025 = q025, q975 = q975), imps = imps)
}









