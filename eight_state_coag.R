library(deSolve)
library(tidyverse)
library(rstan)

#########
# FINAL MODEL
#########
eight <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    cascade_delay <- 1/(1+exp(c-b*t))
    tfpi <- 1/(1+exp(-ct+bt*t))
    prothrombinase <- (prothrombinase_max)/(1+exp(ctt-btt*t))
    r_FIIa <- k_FIIa*cascade_delay*tfpi*((FII/g_FII)^c_FIIa)/(K_FIIa + (FII/g_FII)^c_FIIa)
    r_prothrombinase <- k_prothrombinase*prothrombinase*(FII/g_FII)
    r_AT <- k_AT*FIIa*(AT/g_AT)
    r_clot <- k_clot*FIIa*((Fg/9e-6)^c_clot)/(K_clot + (Fg/9e-6)^c_clot)
    
    r_lys <- k_lys*tPA*(Fn^c_lys)/(K_lys + Fn^c_lys)
    r_PAI <- k_PAI*tPA*PAI
    
    dFII <- g_FII*(-r_FIIa-r_prothrombinase)                 # pct activity
    dFIIa <- r_FIIa -r_AT           # 10s of nmol/L
    dAT <- -g_AT*r_AT            # pct activity
    dFg <- -r_clot                  # 1 mg/dl = 29.41 nmol/L
    dFn <- 1e7*(r_clot - r_lys)           # mm of clot
    dtPA <- -r_PAI                  # 1 ng/mL = 14.29 pmol/L     
    dPAI <- -r_PAI                  # 1 ng/mL = 23.26 pmol/L
    
    list(c(dFII, dFIIa, dAT, dFg, dFn, dtPA, dPAI))
  })
}

parameters <- c(g_FII = 100/1.4e-6, g_AT = 100/3.4e-6,
                c = 5.0, b = 0.4e-1,
                ct = 10, bt = 0.15e-1,
                k_FIIa = 3.5e-9, c_FIIa = 1, K_FIIa = 1.4e-6,
                k_prothrombinase = 1e8, prothrombinase_max = 1.5e-10, ctt = 7, btt = 6e-3,
                k_AT = 1.6e4,
                k_clot = 3.0, c_clot = 1, K_clot = 0.75,
                k_lys = 1, c_lys = 1, K_lys = 0.5,
                k_PAI = 4.5e5)
state      <- c(FII = 1e2, FIIa = 0, AT = 100, Fg = 9e-6, Fn = 0, tPA = 7e-11, PAI = 4e-10)
times      <- seq(0, 1800, by = 30)

out <- ode(y = state, times = times, func = eight, parms = parameters, method = "bdf",  atol = 1e-6, rtol = 1e-5) %>% as.data.frame %>% as_tibble

out %>%
  gather(state, value, -time) %>%
  mutate(state = factor(state, levels = c("FII", "FIIa", "AT", "Fg", "Fn", "tPA", "PAI"))) %>%
  ggplot(aes(time,value)) + geom_line() + facet_grid(state ~ ., scales = "free") +
  geom_line(aes(time,value), color = "red", data = mit2)

#paste("R:", round(out$time[min(which(out$Fn > 0.001))],1), "K:", round(out$time[min(which(out$Fn > 0.2))],1), "MA:", round(max(out$Fn),3)*100, "Ly30:", round(100*(1-tail(out$Fn,1)/max(out$Fn)),1))

#ggplot(data.frame(x = c(0, 1800)), aes(x)) + stat_function(fun = function(t) 1/(1+exp(5-0.4e-1*t)) * 1/(1+exp(-10.0+0.15e-1*t)), geom = "line")

# fit log error
#out %>% gather(state, value, -time) %>% left_join(mit2, by = c("state","time")) %>%
#  mutate(logSqError = log10(value.x)-log10(value.y)) %>% ggplot(aes(logSqError)) + geom_histogram() + facet_grid(. ~ state,scales = "free")

# simulate trajectories from the prior
# N <- 1000
# possible_pars <- tibble(c = rnorm(N,5,5e-1), b = rnorm(N,4e-2,1e-2),
#                         ct = rnorm(N,10,1), bt = rnorm(N,1.5e-2,5e-3),
#                         k_FIIa = rnorm(N,3.5e-9,5e-10), c_FIIa = rnorm(N,1,5e-2), K_FIIa = rnorm(N,1.4e-6,1e-7),
#                         k_clot = rnorm(N,3.0,1e-1),c_clot = rnorm(N,1,1e-2),K_clot = rnorm(N,0.75,1e-2),
#                         k_lys = rnorm(N,1,1e-1),c_lys = rnorm(N,1,1e-2), K_lys = rnorm(N,0.5,1e-2))
# 
# sim_system <- function(c,b,ct,bt,k_FIIa,c_FIIa,K_FIIa,k_clot,c_clot,K_clot,k_lys,c_lys,K_lys) {
#   
#   parameters <- c(g_FII = 100/1.4e-6, g_AT = 100/3.4e-6,
#                   c = c, b = b,
#                   ct = ct, bt = bt,
#                   k_FIIa =k_FIIa, c_FIIa = c_FIIa, K_FIIa = K_FIIa,
#                   k_AT = 1.6e4,
#                   k_clot = k_clot, c_clot = c_clot, K_clot = K_clot,
#                   k_lys = k_lys, c_lys = c_lys, K_lys = K_lys,
#                   k_PAI = 4.5e5)
#   state      <- c(FII = 1e2, FIIa = 0, AT = 100, Fg = 9e-6, Fn = 0, tPA = 7e-11, PAI = 4e-10)
#   times      <- seq(0, 1800, by = 30)
#   
#   ode(y = state, times = times, func = eight, parms = parameters, method = "bdf",  atol = 1e-5, rtol = 1e-6) %>% as.data.frame %>% as_tibble
# }
# 
# possible_pars %>%
#   pmap(sim_system) %>%
#   map2(1:N, function(x,y) x %>% mutate(draw = y)) %>%
#   bind_rows() %>%
#   gather(state,value,-time,-draw) %>%
#   ggplot(aes(time, value, group = draw)) +
#   geom_line(alpha =0.1) +
#   facet_grid(state ~ ., scales ="free") +
#   geom_line(aes(time,value, group=NULL), color = "red", data = mit2)
  
# fit in Stan
mit_stan <- mit %>% filter((time %% 30) == 0)
Nt <- nrow(mit_stan) -1
ts <- mit_stan %>% .$time %>% tail(Nt)
y0 <- mit_stan %>% select(-time) %>% head(1) %>% as.matrix %>% as.vector
y <- mit_stan %>% select(-time) %>% tail(Nt) %>% as.matrix
sigma <- c(1e-2,8e-15, 1e-3, 1e-9, 5e-6, 5e-11, 1e-11)
dat <- list(Nt = Nt, ts = ts, y0 = y0, y = y, sigma = sigma)

fit <- stan(file = "fit_coag_ode.stan", data = dat, chains = 1, iter = 100, refresh = 1, control = list(max_treedepth = 5),
            init = list(list(c_cascade = 5.0, b_cascade = 0.4e-1,
                             c_tfpi = 10, b_tfpi = 0.15e-1,
                             k_FIIa = 3.5e-9, c_FIIa = 1, K_FIIa = 1.4e-6,
                             k_clot = 3.0, c_clot = 1, K_clot = 0.75,
                             k_lys = 1, c_lys = 1, K_lys = 0.5,
                             c_prothrombinase = 7, b_prothrombinase = 6e-3)))

z <- extract(fit, pars = "z")$z %>% melt %>% as_tibble %>% arrange(iterations, Var3,Var2) %>%
  mutate(time = rep(seq(30,1800,by=30),7*50)) %>%
  mutate(state = factor(Var3, labels = c("FII","FIIa","AT","Fg","Fn","tPA","PAI")))

z %>% ggplot(aes(time,value,group=iterations)) + geom_line() + facet_grid(state ~ ., scales = "free") + geom_line(aes(group = NULL), color = "red", data = mit2)

# return estimates
estimates <- extract(fit, pars = c("c_cascade","b_cascade", "c_tfpi", "b_tfpi", "k_FIIa", "K_FIIa","K_FIIa","k_clot","c_clot","K_clot","k_lys","c_lys","K_lys","c_prothrombinase","b_prothrombinase")) %>% map(function(x) x[1])

# m <- stan_model(file = "fit_coag_ode.stan")
# optimizing(m, data = dat, iter = 1, save_iterations = TRUE,
#            init = list(c = 5.0, b = 0.4e-1,
#                        ct = 10, bt = 0.15e-1,
#                        k_FIIa = 3.5e-9, c_FIIa = 1, K_FIIa = 1.4e-6,
#                        k_clot = 3.0, c_clot = 1, K_clot = 0.75,
#                        k_lys = 1, c_lys = 1, K_lys = 0.5))
