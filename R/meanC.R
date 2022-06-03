
# C Function to calculate mean - note: faster, but compromises on numerical accuracy
Rcpp::cppFunction('double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}')
