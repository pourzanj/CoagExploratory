library(deSolve)
library(tidyverse)

v <- function(z, c, eta) (z*c)^eta / (1 + (z*c)^eta)

adi <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    r1 <- k_init*v(TF, c_init, eta_init)*FII/(K_init+FII)
    r2 <- k_amp*v(FX, c_amp, eta_amp)*FII/(K_amp+FII)
    r3 <- k_inh*FIIa^gamma*AT
    r4 <- k_apc*v(FIIa, c_apc, eta_apc)*Tm*PC/(K_apc+PC)
    r5 <- k_clot*v(FIIa, c_clot, eta_clot)*Fg/(K_clot+Fg)
    r6 <- k_lys*v(tPA, c_lys, eta_lys)*FnII2/(K_lys+FnII2)
    r7 <- k_pai*tPA*PAI
    
    #print(paste(t,",", v(FX, c_amp, eta_amp), ",", FII/(K_amp+FII), ",", r2))
    dFII <- - r2
    dFIIa <- r2 - r3
    dPC <- -r4
    dAPC <- r4
    dAT <- -r3
    dFg <- -r5
    dFnII2 <- r5 - r6
    dFDP <- r6
    dtPA <- -r7
    dPAI <- -r7
        
    list(c(dFII, dFIIa, dPC, dAPC, dAT, dFg, dFnII2, dFDP, dtPA, dPAI))
  })
}

parameters <- c(TF = 15e-12, FX = 1.6e-7, Tm = 1e-9,
                k_init = 1e-8, K_init = 1e-6, c_init = 1, eta_init = 1,
                k_amp = 1e-6, K_amp = 1e-6, c_amp = 1e5, eta_amp = 1,
                k_inh = 1e5, gamma = 1.26,
                k_apc = 1e1, K_apc = 1e-8, c_apc = 1e3, eta_apc = 1,
                k_clot = 1e-3, K_clot = 1e-9, c_clot = 1e1, eta_clot = 1,
                k_lys = 1e-1, K_lys = 1e-8, c_lys = 1, eta_lys = 1,
                k_pai = 1e6)

parameters <- c(TF = 15e-12, FX = 1.6e-7, Tm = 1e-9,
                k_init = 1e-8, K_init = 1e-6, c_init = 1, eta_init = 1,
                k_amp = 9.881137e-06, K_amp = 3.503041e-05, c_amp = 1e3, eta_amp = 1.049349,
                k_inh = 1e3, gamma = 1.26,
                k_apc = 10.02335, K_apc = 9.998233e-09, c_apc = 1029.202, eta_apc = 0.6666018,
                k_clot = 1e0, K_clot = 1e-5, c_clot = 1e3, eta_clot = 0.998769,
                k_lys = 0.09999705, K_lys = 1.000029e-08, c_lys = 0.9999705, eta_lys = 1.000564,
                k_pai = 449351.6)

state      <- c(FII = 1.4e-6, FIIa = 0, PC = 7e-8, APC = 0, AT = 3.4e-6,
                Fg = 9e-6, FnII2 = 0, FDP = 0,
                tPA = 5e-9, PAI = 4e-10)
times      <- seq(0, 1800, by = 1)

out_adi <- ode(y = state, times = times, func = adi, parms = parameters, atol = 1e-6, rtol = 1e-18) %>% as.data.frame %>% as_tibble

out_adi %>% gather(state,value,-time) %>%
  ggplot(aes(time, value)) + geom_line() + facet_grid(state ~ ., scales = "free")

##############
out <- out %>% select(time, FII, FIIa, PC, APC, AT, Fg, FnII2, FDP, tPA, PAI)
y0 <- as.matrix(out)[1,2:11]
out_stan <- (out[2:nrow(out),] %>% filter(time %% 10 == 0))

dat <- list(N = nrow(out_stan), ts = out_stan$time, x_r = c(15e-12, 1.6e-7, 1e-9),
            sigma = c(1e-8, 1e-11, 1e-9, 1e-10, 1e-7, 1e-7, 1e-7, 1e-8, 1e-10, 1e-11),
            y0, y = select(out_stan,-time) %>% as.matrix)
fit <- stan(file = "adi_coagulation.stan", data = dat, iter = 2000, chains = 1, refresh = 1,
            init = list(list(k_init = 1e-8, K_init = 1e-6, c_init = 1, eta_init = 1,
                             k_amp = 9.881137e-08, K_amp = 3.503041e-07, c_amp = 98841.25, eta_amp = 1.049349,
                             k_inh = 101298.5, gamma = 0.9347609,
                             k_apc = 10.02335, K_apc = 9.998233e-09, c_apc = 1029.202, eta_apc = 0.6666018,
                             k_clot = 0.001000066, K_clot = 1e-9, c_clot = 10.00066, eta_clot = 0.998769,
                             k_lys = 0.09999705, K_lys = 1.000029e-08, c_lys = 0.9999705, eta_lys = 1.000564,
                             k_pai = 449351.6)))

log_prob(fit, upars = unconstrain_pars(fit, pars = list(k_init = 1e-8, K_init = 1e-6, c_init = 1, eta_init = 1,
                                                       k_amp = 9.881137e-08, K_amp = 3.503041e-07, c_amp = 98841.25, eta_amp = 1.049349,
                                                       k_inh = 101298.5, gamma = 0.9347609,
                                                       k_apc = 10.02335, K_apc = 9.998233e-09, c_apc = 1029.202, eta_apc = 0.6666018,
                                                       k_clot = 0.001000066, K_clot = 1e-9, c_clot = 10.00066, eta_clot = 0.998769,
                                                       k_lys = 0.09999705, K_lys = 1.000029e-08, c_lys = 0.9999705, eta_lys = 1.000564,
                                                       k_pai = 449351.6)))



ll <- function(theta) log_prob(fit, upars = unconstrain_pars(fit, pars = list(k_init = theta[1], K_init = theta[2], c_init = theta[3], eta_init = theta[4],
                                                                         k_amp = theta[5], K_amp = theta[6], c_amp = theta[7], eta_amp = theta[8],
                                                                         k_inh = theta[9], gamma = theta[10],
                                                                         k_apc = theta[11], K_apc = theta[12], c_apc = theta[13], eta_apc = theta[14],
                                                                         k_clot = theta[15], K_clot = theta[16], c_clot = theta[17], eta_clot = theta[18],
                                                                         k_lys = theta[19], K_lys = theta[20], c_lys = theta[21], eta_lys = theta[22],
                                                                         k_pai = theta[23])))

gradll <- function(theta) grad_log_prob(fit, upars = unconstrain_pars(fit, pars = list(k_init = theta[1], K_init = theta[2], c_init = theta[3], eta_init = theta[4],
                                                                              k_amp = theta[5], K_amp = theta[6], c_amp = theta[7], eta_amp = theta[8],
                                                                              k_inh = theta[9], gamma = theta[10],
                                                                              k_apc = theta[11], K_apc = theta[12], c_apc = theta[13], eta_apc = theta[14],
                                                                              k_clot = theta[15], K_clot = theta[16], c_clot = theta[17], eta_clot = theta[18],
                                                                              k_lys = theta[19], K_lys = theta[20], c_lys = theta[21], eta_lys = theta[22],
                                                                              k_pai = theta[23])))

ll <- function(utheta) log_prob(fit, upars = utheta)
gradll <- function(utheta) grad_log_prob(fit, upars = utheta)

pos <- psoptim(rep(1e-1,23), fn = ll, gr = gradll, lower = 0)

ll <- function(theta) {
  p <- list(TF = 15e-12, FX = 1.6e-7, Tm = 1e-9,
            k_init = theta[1], K_init = theta[2], c_init = theta[3], eta_init = theta[4],
       k_amp = theta[5], K_amp = theta[6], c_amp = theta[7], eta_amp = theta[8],
       k_inh = theta[9], gamma = theta[10],
       k_apc = theta[11], K_apc = theta[12], c_apc = theta[13], eta_apc = theta[14],
       k_clot = theta[15], K_clot = theta[16], c_clot = theta[17], eta_clot = theta[18],
       k_lys = theta[19], K_lys = theta[20], c_lys = theta[21], eta_lys = theta[22],
       k_pai = theta[23])
  
  yi <- ode(y = state, times = times, func = adi, parms = p, atol = 1e-6, rtol = 1e-18) %>% as.data.frame %>% as.matrix
  
  res <- t(((as.matrix(out)-yi)[,2:11])^2)/c(1e-8, 1e-11, 1e-9, 1e-10, 1e-7, 1e-7, 1e-7, 1e-8, 1e-10, 1e-11)
  ret <- sum(res)/(1800*10)
  if(is.nan(ret)) ret <- Inf
  
  return(ret)
  
}

m <- stan_model(file = "adi_coagulation.stan")
v <- vb(m, data = dat,
        init = list(k_init = 1e-8, K_init = 1e-6, c_init = 1, eta_init = 1,
                         k_amp = 9.881137e-08, K_amp = 3.503041e-07, c_amp = 98841.25, eta_amp = 1.049349,
                         k_inh = 101298.5, gamma = 0.9347609,
                         k_apc = 10.02335, K_apc = 9.998233e-09, c_apc = 1029.202, eta_apc = 0.6666018,
                         k_clot = 0.001000066, K_clot = 1e-9, c_clot = 10.00066, eta_clot = 0.998769,
                         k_lys = 0.09999705, K_lys = 1.000029e-08, c_lys = 0.9999705, eta_lys = 1.000564,
                         k_pai = 449351.6))
