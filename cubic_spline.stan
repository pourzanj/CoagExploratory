functions {
  
  // solves tridiagonal matrix Ax=f with 2 across
  // across the diagonal, a on the lower diag,
  // and c on the upper diag
  real[] thomas_solve(real[] a, real[] c, real[] f) {
    
    int N = size(f);
    
    real alpha[N];
    real beta[N-1];
    real y[N];
    real x[N];
    
    alpha[1] = 2;
    y[1] = f[1];
    
    for(i in 2:N) {
      beta[i-1] = a[i-1]/alpha[i-1];
      alpha[i] = 2 - beta[i-1]*c[i-1];
  
      y[i] = f[i] - beta[i-1]*y[i-1];
    }
  
    x[N] = y[N]/alpha[N];
    for(i in (N-1):1) {
      x[i] = (y[i]-c[i]*x[i+1])/alpha[i];
      print(x[i]);
    }
    
    return(x);
  } 
  
}

data {
  
  int N;
  real a[N-1];
  real c[N-1];
  real f[N];
  
}

parameters {
  real x;
}

model {
  
  real solved_x[N];
  
  solved_x = thomas_solve(a,c,f);
  
  print("a: ", a);
  print("c: ", c);
  print("f: ", f);
  print("x: ", solved_x);

  x ~ normal(0,1);
}