library(deSolve)
library(tidyverse)

adi <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # set rates
    delay_term = 1 - 1/(1+exp(0.1*(t-60)))
    
    Ftotal <- Fibrin + Fiber
    k_cat_plasminogen_tPA <- k_cat_plasminogen_Fibrin_tPA*(Ftotal)/(K_plasminogen_Fibrin_tPA+Ftotal)
    Km_plasminogen_tPA <- Km_plasminogen_Fibrin_tPA*(Ka_plasminogen_Fibrin_tPA+Ftotal)/(K_plasminogen_Fibrin_tPA+Ftotal)
    
    r <- rep(1, 18)
    
    r[1] <- (k_trigger*TRIGGER*(FII/(K_FII_trigger + FII)))*delay_term
    r[2] <- k_amplification*FIIa*(FII/(K_FII_amplification + FII))
    r[3] <- k_APC_formation*TM*(PC/(K_PC_formation + PC))
    r[4] <- k_inhibition_ATIII*(ATIII)*(FIIa^1.26)
    r[5] <- (FIIa*k_cat_Fibrinogen*Fibrinogen)/(Km_Fibrinogen+Fibrinogen)                           # Cleavage of fibrinopeptides A and/or B to form fibrin monomer
    r[6] <- k_fibrin_monomer_association*(Fibrin_monomer^2)                                         # Protofibril formation through association of fibrin monomers
    r[7] <- k_protofibril_association*(Protofibril^2)                                               # Association of protofibril-protofibril to form fibers
    r[8] <- k_protofibril_monomer_association*(Protofibril)*(Fibrin_monomer)                        # Fibril association with monomer forms protofibril again
    r[9] <- k_fibrin_formation*(Fiber);                                                             # Fibrin growth
    r[10] <- (tPA*k_cat_plasminogen_tPA*Plasminogen)/(Km_plasminogen_tPA+Plasminogen)
    r[11] <- (uPA*k_cat_plasminogen_Fibrin_uPA*Plasminogen)/(Km_plasminogen_Fibrin_uPA+Plasminogen)
    r[12] <- (Plasmin*k_cat_Fibrin*Fibrin)/(Km_Fibrin+Fibrin)
    r[13] <- k_inhibit_plasmin*Plasmin*(antiplasmin)
    r[14] <- k_inhibit_PAI_1_APC*(APC)*(PAI_1)
    r[15] <- k_inhibit_tPA_PAI_1*(tPA)*PAI_1
    r[16] <- k_inhibit_uPA_PAI_1*(uPA)*PAI_1
    r[17] <- (Plasmin*k_cat_fiber*Fiber)/(Km_fiber+Fiber)                                           # Plasmin degrading fiber
    r[18] <- (Plasmin*k_cat_fibrinogen_deg*Fibrinogen)/(Km_fibrinogen_deg+Fibrinogen)
    
    # set controled rates (rv)
    c <- rep(1, 18)
    
    # Initiation control -
    initiation_trigger_term <- (alpha_trigger_activation*TRIGGER^order_trigger_activation)/(1 + (alpha_trigger_activation*TRIGGER^order_trigger_activation))
    initiation_TFPI_term <- 1 - (alpha_trigger_inhibition_TFPI*TFPI^order_trigger_inhibition_TFPI)/(1 + (alpha_trigger_inhibition_TFPI*TFPI^order_trigger_inhibition_TFPI))
    c[1] <- min(initiation_trigger_term,initiation_TFPI_term)
    
    # Amplification FIIa feedback  -
    inhibition_term <- 1 - (alpha_amplification_APC*APC^order_amplification_APC)/(1 + (alpha_amplification_APC*APC^order_amplification_APC))
    inhibition_term_TFPI <- 1 - (alpha_amplification_TFPI*TFPI^order_amplification_TFPI)/(1 +(alpha_amplification_TFPI*TFPI^order_amplification_TFPI))
    factor_product <- FV*FX*FVIII*FIX
    factor_amplification_term <- (0.1*factor_product^2)/(1+(0.1*factor_product^2))
    c[2] <- min(inhibition_term,inhibition_term_TFPI,factor_amplification_term)
    
    # Shutdown term -
    shutdown_term <- (alpha_shutdown_APC*FIIa^order_shutdown_APC)/(1 + (alpha_shutdown_APC*FIIa^order_shutdown_APC))
    c[3] <- shutdown_term
    
    # control_term_fibrinolysis_TAFI -
    control_term_fibrinolysis_TAFI <- 1 - ((alpha_fib_inh_TAFI*TAFI)^order_fib_inh_TAFI)/(1+((alpha_fib_inh_TAFI*TAFI)^order_fib_inh_TAFI))
    if(control_term_fibrinolysis_TAFI > 0.001) c[10] = control_term_fibrinolysis_TAFI
    
    control_term_fibrinolysis_TAFI <- 1 - ((alpha_fib_inh_TAFI*TAFI)^order_fib_inh_TAFI)/(1+((alpha_fib_inh_TAFI*TAFI)^order_fib_inh_TAFI))
    if(control_term_fibrinolysis_TAFI>0.001) c[11] = control_term_fibrinolysis_TAFI
    
    # control_term_fibrinolysis_fXIII -
    control_term_fibrinolysis_fXIII <- 1 - ((alpha_fib_inh_fXIII*fXIII)^order_fib_inh_fXIII)/(1+((alpha_fib_inh_fXIII*fXIII)^order_fib_inh_fXIII))
    if(control_term_fibrinolysis_fXIII>0.001) c[12] = control_term_fibrinolysis_fXIII
    
    #####################
    # get rv and set ODE equations
    rv <- r*c
    
    dFII <- -rv[2] - rv[1]
    dFIIa <- rv[2] + rv[1] - rv[4]
    dPC <- -rv[3]
    dAPC <- rv[3] - rv[14]
    dAT <- -rv[5]
    dTm <- 0.0
    dTrig <- 0.0
    dFn <- rv[9] - rv[12]
    dPn <- rv[10] + rv[11] - rv[13]
    dFg <- -rv[5] - rv[18]
    dPg <- -rv[10] - rv[11]
    dtPA <- -rv[15]
    duPA <- -rv[16]
    dFnI <- rv[5] - rv[6]
    dFnI2 <- rv[6] + rv[8] - rv[7]
    dAP <- -rv[13]
    dPAI <- -rv[14] - rv[15] -rv[16]
    dFnII <- rv[7] - rv[9] - rv[17]
    
    list(c(dFII, dFIIa, dPC, dAPC, dAT, dTm, dTrig, dFn, dPn, dFg, dPg, dtPA, duPA, dFnI, dFnI2, dAP, dPAI, dFnII))
  })
}

parameters <- c(k_trigger = 64.8512, K_FII_trigger = 235.6090, k_amplification = 5.0436e-7,                 # Kinetic parameters for thrombin generation and regulation
                K_FII_amplification = 24.1885, k_APC_formation = 0.0094, K_PC_formation = 6.0898, k_inhibition_ATIII = 2.9463e-6,
                
                k_cat_Fibrinogen = 0.5219 Km_Fibrinogen = 12.5155, k_fibrin_monomer_association = 0.0855,             # Kinetic parameters for fibrin generation
                k_protofibril_association = 0.0011, k_protofibril_monomer_association = 0.0118,
                
                k_cat_plasminogen_Fibrin_tPA = 0.0430, Km_plasminogen_Fibrin_tPA = 0.7994,
                k_cat_plasminogen_Fibrin_uPA = 0.0003, Km_plasminogen_Fibrin_uPA = 0.8068,
                k_cat_Fibrin = 1.5303, Km_Fibrin  = 0.4169, k_inhibit_plasmin = 7.6994e-6,
                k_inhibit_PAI_1_APC = 6.8239e-5, k_inhibit_tPA_PAI_1 = 0.1094, k_inhibit_uPA_PAI_1 = 0.0116,             # Kinetic parameters for fibrinolysis
                
                k_fibrin_formation = 0.0841,                                                              # Kinetic parameter for fibrin formation
                
                k_cat_fiber = 1.4263, Km_fiber = 0.1824,                                                        # Fiber degradation
                
                Ka_plasminogen_Fibrin_tPA = 0.0138, K_plasminogen_Fibrin_tPA = 0.0598,                          # Plasmin activation
                k_cat_fibrinogen_deg = 0.0256, Km_fibrinogen_deg = 1.9721,
                
                alpha_trigger_activation = 0.0413, order_trigger_activation = 0.0162,                                       # initiation
                alpha_trigger_inhibition_TFPI = 0.0203, order_trigger_inhibition_TFPI = 0.1627,
                
                alpha_amplification_APC = 2.1841e-6, order_amplification_APC = 0.0699,                                         # Amplification
                alpha_amplification_TFPI = 0.4244, order_amplification_TFPI = 0.0023,
                
                alpha_shutdown_APC = 0.0247, order_shutdown_APC = 0.0030,                                          # APC generation
                alpha_fib_inh_fXIII = 0.0702, order_fib_inh_fXIII = 0.1229,
                alpha_fib_inh_TAFI = 0.0151, order_fib_inh_TAFI = 0.6620,
                
                alpha_tPA_inh_PAI_1 = 0.2508, order_tPA_inh_PAI_1 = 0.4010,                                      # Control constants for plasmin activation term
                alpha_uPA_inh_PAI_1 = 0.1148, order_uPA_inh_PAI_1 = 0.0346)


state      <- c(dFII = 64.8512, dFIIa = 235.6090, dPC = 5.0436e-7, dAPC = 24.1885, dAT,
                dTm, dTrig, dFn, dPn, dFg,
                dPg, dtPA, duPA, dFnI, dFnI2,
                dAP, dPAI, dFnII)
times      <- seq(0, 1800, by = 1)

out_adi <- ode(y = state, times = times, func = adi, parms = parameters, atol = 1e-6, rtol = 1e-18) %>% as.data.frame %>% as_tibble

out_adi %>% gather(state,value,-time) %>%
  ggplot(aes(time, value)) + geom_line() + facet_grid(state ~ ., scales = "free")