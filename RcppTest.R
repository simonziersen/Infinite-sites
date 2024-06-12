library(Rcpp)
library(microbenchmark)

cppFunction('int add(int x, int y, int z){
            int sum = x + y + z;
            return sum;
            }')
add(1, 2, 3)

sumR <- function(x){
  total <- 0
  for(i in seq_along(x)){
    total <- total + x[i]
  }
  total
}

cppFunction('double sumCpp(NumericVector x){
            int n = x.size();
            double total = 0;
            for(int i =0; i < n; i++){
              total += x[i];
            }
            return total;
            }')

cppFunction('double sumC(NumericVector x) {
  int n = x.size();
            double total = 0;
            for(int i = 0; i < n; ++i) {
            total += x[i];
            }
            return total;
            }')

x <- runif(1e3)
microbenchmark(
  sumR(x),
  sumCpp(x),
  sumC(x),
  sum(x)
)

pdistR <- function(x, ys){
  sqrt((x - ys)^2)
}

cppFunction('NumericVector pdistC(double x, NumericVector ys) {
  int n = ys.size();
  NumericVector out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = sqrt(pow(ys[i] - x, 2.0));
  }
  return out;
}')

cppFunction('NumericVector pdistCpp(double x, NumericVector y){
            int n = y.size();
            NumericVector out(n);
            for(int i = 0; i < n; i++){
              out[i] = sqrt(pow(x - y[i], 2.0));
            }
            return out;
            }')
x <- 9.7
y <- runif(1e5)

microbenchmark(
  pdistR(x, y),
  pdistCpp(x, y),
  pdistC(x, y)
)

cppFunction('NumericVector rowSumsCpp(NumericMatrix x){
            int nrow = x.nrow(), ncol = x.ncol();
            NumericVector out(nrow);
            
            for(int i = 0; i < nrow; i++){
              out [i] = 0;
              for(int j = 0; j < ncol; j++){
                out[i] += x(i, j);
              }
            }
            
            return out;
            }')

cppFunction('NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}')

set.seed(1014)
x <- matrix(sample(100), 10)

microbenchmark(
  rowSumsCpp(x),
  rowSumsC(x),
  rowSums(x)
)

x <- runif(1e5)
sourceCpp("cppFunctions.cpp")

meanR <- function(x){
  y <- 0
  for(i in seq_along(x)){
   y <- y + x[i] 
  }
  y / length(x)
}

microbenchmark(
  meanR(x),
  mean(x),
  meanC(x),
  meanCto(x)
)





