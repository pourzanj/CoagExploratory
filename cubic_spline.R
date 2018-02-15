fit_cubic_spline <- function(x, f) {
  
  N <- length(x)-1
  h <- diff(x)
  
  mu <- h[1:(N-1)] / (h[1:(N-1)] + h[2:N])
  lambda <- h[2:N] / (h[1:(N-1)] + h[2:N])
  d <- 6/(h[1:(N-1)]+h[2:N]) * ((f[3:(N+1)]-f[2:N])/h[2:N]-(f[2:N]-f[1:(N-1)])/h[1:(N-1)])
  
  # extend vectors because there are still unknowns
  M <- thomas_solve(c(mu,1), c(1,lambda), c(d[1],d,d[N-1]))
  
  C <- (f[2:(N+1)]-f[1:N])/h - (h/6)*(M[2:(N+1)]-M[1:N])
  Ct <- f[1:N] - M[1:N]*(h^2/6)
  
  function(t) {
    i <- 1
    if(t <= x[1]) i <- 1
    else if(t >= x[N+1]) i <- N
    else i <- min(which(t < x)) - 1
    M[i]*(x[i+1]-t)^3/(6*h[i]) + M[i+1]*(t-x[i])^3/(6*h[i]) + C[i]*(t-x[i]) + Ct[i]
  }
  
}

fit_cubic_spline <- function(x, f) {
  
  N <- length(x)-1
  h <- diff(x)
  
  mu <- h[1:(N-1)] / (h[1:(N-1)] + h[2:N])
  lambda <- h[2:N] / (h[1:(N-1)] + h[2:N])
  d <- 6/(h[1:(N-1)]+h[2:N]) * ((f[3:(N+1)]-f[2:N])/h[2:N]-(f[2:N]-f[1:(N-1)])/h[1:(N-1)])
  
  # extend vectors because there are still unknowns
  M <- thomas_solve(c(mu,1), c(1,lambda), c(d[1],d,d[N-1]))
  
  C <- (f[2:(N+1)]-f[1:N])/h - (h/6)*(M[2:(N+1)]-M[1:N])
  Ct <- f[1:N] - M[1:N]*(h^2/6)
  
  function(t) {
    i <- 1
    if(t <= x[1]) i <- 1
    else if(t >= x[N+1]) i <- N
    else i <- min(which(t < x)) - 1
    
    M[i]*(x[i+1]-t)^3/(6*h[i]) + M[i+1]*(t-x[i])^3/(6*h[i]) + C[i]*(t-x[i]) + Ct[i]
    a <- -M[i]/(6*h[i]) + M[i+1]/(6*h[i])
    b <- (M[i]*x[i+1])/(2*h[i]) - (M[i+1]*x[i])/(2*h[i])
    c <- -(M[i]*x[i+1]^2)/(2*h[i]) + (M[i+1]*x[i]^2)/(2*h[i]) + C[i]
    d <- (M[i]*x[i+1]^3)/(6*h[i]) - (M[i+1]*x[i]^3)/(6*h[i]) - C[i]*x[i] + Ct[i]
    
    return(list(a = a, b=b, c=c, d=d))
  }
  
}


thomas_solve <- function(b, c, f) {
  
  N <- length(f)
  alpha <- rep(NA, N)
  beta <- rep(NA, N-1)
  y <- rep(NA,N)
  x <- rep(NA,N)
  
  alpha[1] <- 2
  y[1] <- f[1]
  for(i in 2:N) {
    beta[i-1] <- b[i-1]/alpha[i-1]
    alpha[i] <- 2 - beta[i-1]*c[i-1]
    
    y[i] <- f[i] - beta[i-1]*y[i-1]
  }
  
  x[N] <- y[N]/alpha[N]
  for(i in (N-1):1) x[i] <- (y[i]-c[i]*x[i+1])/alpha[i]
  
  return(x)
}

# get_cubic_root <- function(a,b,c,d) {
#   Delta0 <- b^2 - 3*a*c
#   Delta1 <- 2*b^3 - 9*a*b*c + 27*a^2*d
#   
#   C1 <- ((Delta1 + sqrt(as.complex(Delta1^2-4*Delta0^3)))/2)^(1/3)
#   C2 <- (-1/2 + 1/2*sqrt(3)*complex(imaginary=1))*C1
#   C3 <- (-1/2 - 1/2*sqrt(3)*complex(imaginary=1))*C1
#   
#   return(list(-1/(3*a)*(b+C1+Delta0/C1), -1/(3*a)*(b+C2+Delta0/C2), -1/(3*a)*(b+C3+Delta0/C3)))
# }

# to test, compare against the following
library(limSolve)
nn   <- 20                          # nr rows and columns of A
aa   <- runif(nn-1)
bb   <- rep(2,20)
cc   <- runif(nn-1)
B    <- rnorm(nn)
Solve.tridiag(aa, bb, cc, B)
thomas_solve(aa,cc,B)
stan(file = "cubic_spline.stan", data = list(N = 20, a = aa, c = cc, f = B), iter = 1, chains = 1)

Fnt <- fit_cubic_spline(out$time, out$Fn)
#for(t in out$time) print(Fnt(t))
out %>% select(time, Fn) %>% ggplot(aes(time,Fn)) + geom_point() + stat_function(fun = Vectorize(Fnt), color = "red") +
  stat_function(fun = function(x) -0.001109547*x^3 + 0.02269766*x^2 + 0.1780197*x + -0.09963753) +
  xlim(0,2)