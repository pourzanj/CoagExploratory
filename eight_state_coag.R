library(deSolve)
library(tidyverse)

#########
# FINAL MODEL
#########
eight <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    cascade_delay <- 1/(1+exp(c-b*t))
    r_FIIa <- k_FIIa*cascade_delay*(FII^c_FIIa)/(K_FIIa + FII^c_FIIa)
    r_AT <- k_AT*FIIa*AT
    r_clot <- k_clot*FIIa*(Fg^c_clot)/(K_clot + Fg^c_clot)
    r_lys <- k_lys*tPA*(Fn^c_lys)/(K_lys + Fn^c_lys)
    r_PAI <- k_PAI*tPA*PAI
    
    dFII <- -g_FII*r_FIIa                 # pct activity
    dFIIa <- r_FIIa -r_AT           # 10s of nmol/L
    dAT <- -g_AT*r_AT            # pct activity
    dFg <- -r_clot                  # 1 mg/dl = 29.41 nmol/L
    dFn <- r_clot - r_lys           # mm of clot
    dtPA <- -r_PAI                  # 1 ng/mL = 14.29 pmol/L     
    dPAI <- -r_PAI                  # 1 ng/mL = 23.26 pmol/L
    
    list(c(dFII, dFIIa, dAT, dFg, dFn, dtPA, dPAI))
  })
}

parameters <- c(c = 0, b = 100,
                k_FIIa = 1e-8, c_FIIa = 1, K_FIIa = 1, g_FII = 1e10,
                k_AT = 1e-1, pctM_AT = 1e-2, g_AT = 1e2,
                k_clot = 1e-2, c_clot = 1, K_clot = 1,
                k_lys = 1, c_lys = 10, K_lys = 0.5,
                k_PAI = 4e7)
state      <- c(FII = 100, FIIa = 0, AT = 100, Fg = 1, Fn = 0, tPA = 1.429e-10, PAI = 2.326e-10)
times      <- seq(0, 30, by = 0.5)

out <- ode(y = state, times = times, func = eight, parms = parameters, atol = 1e-2, rtol = 1e-2) %>% as.data.frame %>% as_tibble

out %>%
  gather(state, value, -time) %>%
  mutate(state = factor(state, levels = c("FII", "FIIa", "AT", "Fg", "Fn", "tPA", "PAI"))) %>%
  ggplot(aes(time,value)) + geom_line() + facet_grid(state ~ ., scales = "free")

paste("R:", round(out$time[min(which(out$Fn > 0.001))],1), "K:", round(out$time[min(which(out$Fn > 0.2))],1), "MA:", round(max(out$Fn),3)*100, "Ly30:", round(100*(1-tail(out$Fn,1)/max(out$Fn)),1))

