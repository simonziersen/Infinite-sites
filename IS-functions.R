



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


