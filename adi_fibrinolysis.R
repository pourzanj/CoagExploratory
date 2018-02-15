library(deSolve)
library(tidyverse)

adi <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #############
    # set rates
    delay_term = 1 - 1/(1+exp(0.1*(t-60)))
    
    Ftotal <- Fn + FnII
    k_cat_plasminogen_tPA <- k_cat_plasminogen_Fibrin_tPA*(Ftotal)/(K_plasminogen_Fibrin_tPA+Ftotal)
    Km_plasminogen_tPA <- Km_plasminogen_Fibrin_tPA*(Ka_plasminogen_Fibrin_tPA+Ftotal)/(K_plasminogen_Fibrin_tPA+Ftotal)
    
    r <- rep(1, 18)
    
    r[1] <- (k_trigger*Trig*(FII/(K_FII_trigger + FII)))*delay_term
    r[2] <- k_amplification*FIIa*(FII/(K_FII_amplification + FII))
    r[3] <- k_APC_formation*Tm*(PC/(K_PC_formation + PC))
    r[4] <- k_inhibition_ATIII*(AT)*(FIIa^1.26)
    r[5] <- (FIIa*k_cat_Fibrinogen*Fg)/(Km_Fibrinogen+Fg)                           # Cleavage of fibrinopeptides A and/or B to form fibrin monomer
    r[6] <- k_fibrin_monomer_association*(FnI^2)                                         # Protofibril formation through association of fibrin monomers
    r[7] <- k_protofibril_association*(FnI2^2)                                               # Association of protofibril-protofibril to form fibers
    r[8] <- k_protofibril_monomer_association*(FnI2)*(FnI)                        # Fibril association with monomer forms protofibril again
    r[9] <- k_fibrin_formation*(FnII);                                                             # Fibrin growth
    r[10] <- (tPA*k_cat_plasminogen_tPA*Pg)/(Km_plasminogen_tPA+Pg)
    r[11] <- (uPA*k_cat_plasminogen_Fibrin_uPA*Pg)/(Km_plasminogen_Fibrin_uPA+Pg)
    r[12] <- (Pn*k_cat_Fibrin*Fn)/(Km_Fibrin+Fn)
    r[13] <- k_inhibit_plasmin*Pn*(AP)
    r[14] <- k_inhibit_PAI_1_APC*(APC)*(PAI)
    r[15] <- k_inhibit_tPA_PAI_1*(tPA)*PAI
    r[16] <- k_inhibit_uPA_PAI_1*(uPA)*PAI
    r[17] <- (Pn*k_cat_fiber*FnII)/(Km_fiber+FnII)                                           # Plasmin degrading fiber
    r[18] <- (Pn*k_cat_fibrinogen_deg*Fg)/(Km_fibrinogen_deg+Fg)
    
    #print(paste(FIIa, (FII/(K_FII_amplification + FII)),  r[2]))
    
    #######################
    # set controled rates (v)
    c <- rep(1, 18)
    
    # Initiation control -
    initiation_trigger_term <- (alpha_trigger_activation*Trig^order_trigger_activation)/(1 + (alpha_trigger_activation*Trig^order_trigger_activation))
    initiation_TFPI_term <- 1 - (alpha_trigger_inhibition_TFPI*TFPI^order_trigger_inhibition_TFPI)/(1 + (alpha_trigger_inhibition_TFPI*TFPI^order_trigger_inhibition_TFPI))
    c[1] <- min(initiation_trigger_term,initiation_TFPI_term)
    
    # Amplification FIIa feedback  -
    inhibition_term <- 1 - (alpha_amplification_APC*APC^order_amplification_APC)/(1 + (alpha_amplification_APC*APC^order_amplification_APC))
    inhibition_term_TFPI <- 1 - (alpha_amplification_TFPI*TFPI^order_amplification_TFPI)/(1 +(alpha_amplification_TFPI*TFPI^order_amplification_TFPI))
    factor_product <- FV*FX*FVIII*FIX
    factor_amplification_term <- (0.1*factor_product^2)/(1+(0.1*factor_product^2))
    c[2] <- min(inhibition_term,inhibition_term_TFPI,factor_amplification_term)
    
    #print(paste("~~~~~",inhibition_term,(alpha_amplification_APC*APC^order_amplification_APC)))
    
    # Shutdown term -
    shutdown_term <- (alpha_shutdown_APC*FIIa^order_shutdown_APC)/(1 + (alpha_shutdown_APC*FIIa^order_shutdown_APC))
    c[3] <- shutdown_term
    
    #print(paste("--------",r[3], c[3], (PC/(K_PC_formation + PC))))
    
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

tpa_level <- 200.0
dilution_factor <- 0.84
plasminogen_mean_level <- 1.23
mixing_delay_time <- 372.0

parameters <- c(k_trigger = 64.8512, K_FII_trigger = 235.6090, k_amplification = 5.0436e-7,                 # Kinetic parameters for thrombin generation and regulation
                K_FII_amplification = 24.1885, k_APC_formation = 0.0094, K_PC_formation = 6.0898, k_inhibition_ATIII = 2.9463e-6,
                
                k_cat_Fibrinogen = 0.5219, Km_Fibrinogen = 12.5155, k_fibrin_monomer_association = 0.0855,             # Kinetic parameters for fibrin generation
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
                alpha_uPA_inh_PAI_1 = 0.1148, order_uPA_inh_PAI_1 = 0.0346,
                
                TFPI = 0.0025*1000*dilution_factor, FV = 0.02*1000*dilution_factor,
                FVIII = 0.0007*1000*dilution_factor, FIX = 0.09*1000*dilution_factor,
                FX = 0.17*1000*dilution_factor, fXIII = 0.093*1000*dilution_factor,
                TAFI = 51*dilution_factor)

parameters <- c(k_trigger = 130, K_FII_trigger = 1400, k_amplification = 1e-4,                 # Kinetic parameters for thrombin generation and regulation
                K_FII_amplification = 100, k_APC_formation = 0.001, K_PC_formation = 60, k_inhibition_ATIII = 1e-7,
                
                k_cat_Fibrinogen = 0.01, Km_Fibrinogen = 1e3, k_fibrin_monomer_association = 0.01,             # Kinetic parameters for fibrin generation
                k_protofibril_association = 0.01, k_protofibril_monomer_association = 0.01,
                
                k_cat_plasminogen_Fibrin_tPA = 2, Km_plasminogen_Fibrin_tPA = 10,
                k_cat_plasminogen_Fibrin_uPA = 0.1, Km_plasminogen_Fibrin_uPA = 10,
                k_cat_Fibrin = 0, Km_Fibrin  = 10, k_inhibit_plasmin = 0.01,
                k_inhibit_PAI_1_APC = 0.01, k_inhibit_tPA_PAI_1 = 0.01, k_inhibit_uPA_PAI_1 = 0.01,             # Kinetic parameters for fibrinolysis
                
                k_fibrin_formation = 100,                                                              # Kinetic parameter for fibrin formation
                
                k_cat_fiber = 10, Km_fiber = 0.01,                                                        # Fiber degradation
                
                Ka_plasminogen_Fibrin_tPA = 100, K_plasminogen_Fibrin_tPA = 10,                          # Plasmin activation
                k_cat_fibrinogen_deg = 0.01, Km_fibrinogen_deg = 0.01,
                
                alpha_trigger_activation = 1, order_trigger_activation = 10,                                       # initiation
                alpha_trigger_inhibition_TFPI = 0.1, order_trigger_inhibition_TFPI = 2,
                
                alpha_amplification_APC = 0.1, order_amplification_APC = 2,                                         # Amplification
                alpha_amplification_TFPI = 0.1, order_amplification_TFPI = 2,
                
                alpha_shutdown_APC = 0.1, order_shutdown_APC = 2,                                          # APC generation
                alpha_fib_inh_fXIII = 0.1, order_fib_inh_fXIII = 2,
                alpha_fib_inh_TAFI = 0.1, order_fib_inh_TAFI = 2,
                
                alpha_tPA_inh_PAI_1 = 0.1, order_tPA_inh_PAI_1 = 2,                                      # Control constants for plasmin activation term
                alpha_uPA_inh_PAI_1 = 0.1, order_uPA_inh_PAI_1 = 2,
                
                TFPI = 0.0025*1000*dilution_factor, FV = 0.02*1000*dilution_factor,
                FVIII = 0.0007*1000*dilution_factor, FIX = 0.09*1000*dilution_factor,
                FX = 0.17*1000*dilution_factor, fXIII = 0.093*1000*dilution_factor,
                TAFI = 51*dilution_factor)

# parameters <- c(k_trigger = 31.522365931674972, K_FII_trigger = 66.81795370664223, k_amplification = 19248779084764403e-23,                 # Kinetic parameters for thrombin generation and regulation
#                 K_FII_amplification = 17.308506361938452, k_APC_formation = .006199294881536407, K_PC_formation = 5.0610534447782305, k_inhibition_ATIII = 2310156904318893e-21,
#                 
#                 k_cat_Fibrinogen = .16188976648279985, Km_Fibrinogen = 24.45543734015879, k_fibrin_monomer_association = .25655649184268553,             # Kinetic parameters for fibrin generation
#                 k_protofibril_association = .004451943771731852, k_protofibril_monomer_association = .027706018198975522,
#                 
#                 k_cat_plasminogen_Fibrin_tPA = .3880304033521856, Km_plasminogen_Fibrin_tPA = .4878425451398645,
#                 k_cat_plasminogen_Fibrin_uPA = .0031162176415766166, Km_plasminogen_Fibrin_uPA = 2.298742758964309,
#                 k_cat_Fibrin = .2131564706323077, Km_Fibrin  = .17911424963906494, k_inhibit_plasmin = 6336828083835846e-20,
#                 k_inhibit_PAI_1_APC = .0015166007520712082, k_inhibit_tPA_PAI_1 = .2714678874308177, k_inhibit_uPA_PAI_1 = .03486133882977271,             # Kinetic parameters for fibrinolysis
#                 
#                 k_fibrin_formation = .10159009689784514,                                                              # Kinetic parameter for fibrin formation
#                 
#                 k_cat_fiber = 4.660009283955812, Km_fiber = .45241905666115917,                                                        # Fiber degradation
#                 
#                 Ka_plasminogen_Fibrin_tPA = .23394261195757574, K_plasminogen_Fibrin_tPA = .22879945768429272,                          # Plasmin activation
#                 k_cat_fibrinogen_deg = .3743746360155661, Km_fibrinogen_deg = .6828207498270956,
#                 
#                 alpha_trigger_activation = .1493648990054381, order_trigger_activation = .23129096317811593,                                       # initiation
#                 alpha_trigger_inhibition_TFPI = .004465018936605084, order_trigger_inhibition_TFPI = .5108810133563065,
#                 
#                 alpha_amplification_APC = 4525360645172537e-21, order_amplification_APC = .6669339709099888,                                         # Amplification
#                 alpha_amplification_TFPI = .5356018392355524, order_amplification_TFPI = .10369765372469193,
#                 
#                 alpha_shutdown_APC = .0466641204015549, order_shutdown_APC = .001215905574245854,                                          # APC generation
#                 alpha_fib_inh_fXIII = .3930243765689481, order_fib_inh_fXIII =.1649112891014183,
#                 alpha_fib_inh_TAFI = .11427276943665514, order_fib_inh_TAFI = .30327044473478815,
#                 
#                 alpha_tPA_inh_PAI_1 = .3465978388748791, order_tPA_inh_PAI_1 = .5606724802749996,                                      # Control constants for plasmin activation term
#                 alpha_uPA_inh_PAI_1 = .10827243009061041, order_uPA_inh_PAI_1 = .4407980060233813,
#                 
#                 TFPI = 0.0025*1000*dilution_factor, FV = 0.02*1000*dilution_factor,
#                 FVIII = 0.0007*1000*dilution_factor, FIX = 0.09*1000*dilution_factor,
#                 FX = 0.17*1000*dilution_factor, fXIII = 0.093*1000*dilution_factor,
#                 TAFI = 51*dilution_factor)

state      <- c(FII = 1.4*1000*dilution_factor, FIIa = 0, PC = 60*dilution_factor, APC = 0.0*dilution_factor , AT = 3.4*1000*dilution_factor,
                Tm = 0.001*1000*dilution_factor, Trig = 0.005*1000*dilution_factor, Fn = 0*dilution_factor, Pn = 0*dilution_factor, Fg = 10.00*1000*dilution_factor,
                Pg = 1.4*1000*dilution_factor, tPA = tpa_level*dilution_factor, uPA = 0, FnI = 0, FnI2 = 0,
                AP = 1.18*1000*dilution_factor, PAI = 0.00056*1000*dilution_factor, FnII = 0)
times      <- seq(0, 1800, by = 1)

out_adi <- ode(y = state, times = times, func = adi, parms = parameters, atol = 1e-5, rtol = 1e-5, maxsteps = 1e5) %>% as.data.frame %>% as_tibble %>%
  mutate(TEG = Fn + FnI + FnI2 + FnII)

out_adi %>% gather(state,value,-time) %>%
  ggplot(aes(time, value)) + geom_line() + facet_grid(state ~ ., scales = "free")
