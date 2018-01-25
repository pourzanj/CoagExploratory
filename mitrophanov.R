library(deSolve)
library(tidyverse)


mitraphanov <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # TF related rates
    r1   <- k2*TF*FVII - k1*TF_FVII
    r2   <- k4*TF*FVIIa - k3*TF_FVIIa
    r3   <- k5*TF_FVIIa*FVII
    r4   <- k6*FXa*FVII
    r5   <- k7*FIIa*FVII
    r6A  <- k9*TF_FVIIa*FX - k8*TF_FVIIa_FX
    r6B  <- k10*TF_FVIIa_FX
    r7   <- k12*TF_FVIIa*FXa - k11*TF_FVIIa_FXa
    r8A  <- k14*TF_FVIIa*FIX - k13*TF_FVIIa_FIX
    r8B  <- k15*TF_FVIIa*FIXa
    r9   <- k16*FXa*FII
    r10  <- k17*FIIa*FVIII
    r11  <- k19*FVIIIa*FIXa - k18*FIXa_FVIIIa
    r12A <- k21*FIXa_FVIIIa*FX - k20*FIXa_FVIIIa_FX
    r12B <- k22*FIXa_FVIIIa_FX
    r13  <- k24*FVIIIa - k23*FVIIIa1L*FVIIIa2
    r14  <- k25*FIXa_FVIIIa_FX
    r15  <- k25*FIXa_FVIIIa
    r16  <- k26*FIIa*FV
    r17  <- k28*FXa*FVa - k27*FXa_FVa
    r18A <- k30*FXa_FVa*FII - k29*FXa_FVa_FII
    r18B <- k31*FXa_FVa_FII
    r19  <- k32*mIIa*FXa_FVa
    r20  <- k34*FXa*TFPI - k33*FXa_TFPI
    r21  <- k36*TF_FVIIa_FXa*TFPI - k35*TF_FVIIa_FXa_TFPI
    r22  <- k37*TF_FVIIa*FXa*TFPI
    r23  <- k38*FXa*AT
    r24  <- k39*mIIa*AT
    r25  <- k40*FIXa*AT
    r26  <- k41*FIIa*AT
    r27  <- k42*TF_FVIIa*AT
    r28  <- k43*FIXa*FX
    r29  <- k44*mIIa*FV
    r30A <- k46*Fg*FIIa - k45*Fg_FIIa
    r30B <- k47*Fg_FIIa
    r31A <- k46*FnI*FIIa - k48*FnI_FIIa
    r31B <- k49*FnI_FIIa
    r32  <- k51*FnI*FnI - k50*FnI2
    r33A <- k46*FnI2*FIIa - k52*FnI2_FIIa
    r33B <- k53*FnI2_FIIa
    r34  <- k46*FnII*FIIa - k54*FnII_FIIa
    r35  <- k55*FnI2_FIIa*AT
    r36  <- k55*FnI_FIIa*AT
    r37  <- k56*FnII_FIIa*AT
    r38  <- k57*Pn*AP
    r39  <- k58*tPA*PAI
    
    CFnII <- FnII + 2*FnII2
    r40  <- (k59*tPA*Pg*CFnII/(K1+CFnII))/(Pg+K2*(K3+CFnII)/(K1+CFnII))
    r41  <- (k60*Pn*FnII)/(K4 + FnII)
    r42  <- (k60*Pn*FnII2)/(K4+FnII2)
    
    r43  <- k62*Tm*FIIa - k61*Tm_FIIa
    r44A <- k64*Tm_FIIa*PC - k63*Tm_FIIa_PC
    r44B <- k65*Tm_FIIa_PC
    r45  <- k66*Tm_FIIa*AT
    r46  <- k68*APC*FVa - k67*APC_FVa
    r47  <- k69*APC_FVa
    
    r48  <- k70*APC_FVa
    r49A <- k68*APC*FVa5 - k67*APC_FVa5
    r49B <- k70*APC_FVa5
    r50A <- k68*APC*FVa3 - k67*APC_FVa3
    r50B <- k69*APC*FVa53
    r51  <- k71*FVa3
    r52  <- k71*FVa53
    r53  <- k68*APC*LCA1 - k67*APC_LCA1
    r54  <- k64*APC*Tm_FIIa - k63*Tm_FIIa_APC
    r55  <- k28*FXa*FVa5 - k72*FXa_FVa5
    r56  <- k28*FXa*FVa3 - k72*FXa_FVa3
    r57A <- k30*FXa_FVa5*FII - k29*FXa_FVa5_FII
    r57B <- k73*FXa_FVa5_FII
    r58A <- k30*FXa_FVa3*FII - k29*FXa_FVa3_FII
    r58B <- k74*FXa_FVa3_FII
    r59  <- k75*FXa_FVa5*mIIa
    r60  <- k76*FXa_FVa3*mIIa
    r61  <- k77*FXa_FVa3
    r62  <- k77*FXa_FVa3_FII
    r63  <- k62*Tm*mIIa - k61*Tm_mIIa
    r64A <- k64*Tm_mIIa*PC - k63*Tm_mIIa_PC
    r64B <- k65*Tm_mIIa_PC
    r65  <- k66*Tm_mIIa*AT
    r66  <- k28*FXa*FVa53
    r67A <- k30*FXa_FVa53*FII - k29*FXa_FVa53_FII
    r67B <- k74*FXa_FVa53_FII
    r68  <- k76*FXa_FVa53*mIIa
    r69  <- k77*FXa_FVa53
    r70  <- k77*FXa_FVa53_FII
    r71  <- k79*FII*FVa - k78*FII_FVa
    r72  <- k80*FXa_FVa5*APC
    
    # cascade
    dTF <- -r1 - r2
    dFVII <- -r1 - r4 - r5
    dTF_FVII <- r1
    dTF_FVIIa <- r2 - r6A - r7 - r8A + r8B - r22 - r27
    dFVIIa <- -r2 + r3 + r4 + r5
    dFX <- -r6A - r12A + r14- r28
    dFXa <- -r7 + r12B - r17 - r20 - r23 + r28 - r55 - r56 + r61 + r62 - r66 + r69 + r70
    dFIIa <- r9 + r19 - r26 - r30A + r30B - r31A + r31B - r33A + r33B - r34 - r43 + r59 + r60 + r68
    dTF_FVIIa_FX <- r6A - r6B
    dTF_FVIIa_FXa <- r6B + r7 - r21
    dFIX <- -r8A
    dFIXa <- r8B - r11 + r14 + r15 - r25
    dTF_FVIIa_FIX <- r8A - r8B
    dFII <- -r9 - r18A - r57A - r58A + r62 - r67A + r70 - r71
    dFVIII <- -r10
    dFVIIIa <- r10 - r11 - r13
    dFIXa_FVIIIa <- r11 - r12A + r12B - r15
    dFIXa_FVIIIa_FX <- r12A - r12B - r14
    dFVIIIa1L <- r13 + r14 + r15
    dFVIIIa2 <- r13 + r14 + r15
    dFV <- -r16 - r29
    dFVa <- r16 - r17 + r29 - r46 - r71
    dFXa_FVa <- r17 - r18A + r18B
    dFXa_FVa_FII <- r18A - r18B
    dmIIa <- r18B - r19 - r24 + r57B + r58B - r59 - r60 - r63 - r68
    dTFPI <- -r20 - r21
    dFXa_TFPI <- r20 - r22
    dTF_FVIIa_FXa_TFPI <- r21 + r22
    dAT <- -r23 - r24 - r25 -r26 - r27 - r35 - r36 - r37 - r45 - r65
    dFg <- -r30A
    dFg_FIIa <- r30A - r30B
    dFnI <- r30B - r31A -2*r32
    dFnI_FIIa <- r31A - r31B - r36
    dFnII <- r31B + -r34 - r41
    dFnII_FIIa <- r34 - r37
    dFnI2 <- r32 - r33A
    dFnI2_FIIa <- r33A - r33B - r35
    dFnII2 <- r33B - r42
    dPn <- -r38 + r40
    dAP <- -r38
    dtPA <- -r39
    dPAI <- -r39
    dPg <- -r40
    dFDP <- r41 + 2*r42
    dTm <- -r43 + r45 - r63 + r65
    dTm_FIIa <- r43 - r44A + r44B - r45 - r54
    dPC <- -r44A - r64A
    dTm_FIIa_PC <- r44A - r44B
    dAPC <- r44B - r46 + r47 + r48 - r49A + r49B - r50A + r50B - r53 - r54 + r64B
    dAPC_FVa <- r46 - r47 - r48
    dFVa5 <- r47 - r49A - r55
    dFVa3 <- r48 - r50A - r51 - r56
    dAPC_FVa5 <- r49A - r49B
    dAPC_FVa3 <- r50A - r50B
    dFVa53 <- r49B + r50B - r52 - r66
    dLCA1 <- r51 + r52 - r53 + r61 + r62 + r69 + r70
    dFXa_FVa5 <- r55 - r57A + r57B - r72
    dFXa_FVa3 <- r56 - r58A + r58B - r61
    dFXa_FVa5_FII <- r57A - r57B
    dFXa_FVa3_FII <- r58A - r58B - r62
    dTm_mIIa <- r63 - r64A + r64B - r65
    dTm_mIIa_PC <- r64A - r64B
    dFXa_FVa53 <- r66 - r67A + r67B - r69 + r72
    dFXa_FVa53_FII <- r67A - r67B - r70
    
    dAPC_LCA1 <- r53
    dTm_FIIa_APC <- r54
    dFII_FVa <- r71 
    
    list(c(dTF, dFVII, dTF_FVII, dTF_FVIIa, dFVIIa, dFX, dFXa, dFIIa, dTF_FVIIa_FX,
           dTF_FVIIa_FXa, dFIX, dFIXa, dTF_FVIIa_FIX, dFII, dFVIII, dFVIIIa, dFIXa_FVIIIa,
           dFIXa_FVIIIa_FX, dFVIIIa1L, dFVIIIa2, dFV, dFVa, dFXa_FVa, dFXa_FVa_FII,
           dmIIa, dTFPI, dFXa_TFPI, dTF_FVIIa_FXa_TFPI, dAT, dFg, dFg_FIIa, dFnI,
           dFnI_FIIa, dFnII, dFnII_FIIa, dFnI2, dFnI2_FIIa, dFnII2, dPn, dAP, dtPA,
           dPAI, dPg, dFDP, dTm, dTm_FIIa, dPC, dTm_FIIa_PC, dAPC, dAPC_FVa,
           dFVa5, dFVa3, dAPC_FVa5, dAPC_FVa3, dFVa53, dLCA1, dFXa_FVa5, dFXa_FVa3,
           dFXa_FVa5_FII, dFXa_FVa3_FII, dTm_mIIa, dTm_mIIa_PC, dFXa_FVa53, dFXa_FVa53_FII,
           dAPC_LCA1, dTm_FIIa_APC, dFII_FVa))
  })
}



parameters <- c(k1 = 3.1e-3, k2 = 3.2e6,
                k3 = 3.1e-3, k4 = 2.3e7,
                k5 = 4.4e5,
                k6 = 1.3e7,
                k7 = 2.3e4,
                k8 = 1.05, k9 = 2.5e7, k10 = 6,
                k11 = 19, k12 = 2.2e7,
                k13 = 2.4, k14 = 1e7, k15 = 1.8,
                k16 = 7.5e3,
                k17 = 2e7,
                k18 = 5e-3, k19 = 1e7,
                k20 = 1e-3, k21 = 1e8, k22 = 8.2,
                k23 = 2.2e4, k24 = 6e-3,
                k25 = 1e-3,
                k26 = 2e7,
                k27 = 0.2, k28 = 4e8,
                k29 = 103, k30 = 1e8, k31 = 63.5,
                k32 = 2.3e8,
                k33 = 3.6e-4, k34 = 9e5,
                k35 = 1.1e-4, k36 = 3.2e8,
                k37 = 5e7,
                k38 = 4.2e3,
                k39 = 7.1e3,
                k40 = 4.9e2,
                k41 = 1.6e4,
                k42 = 2.3e2,
                k43 = 5.7e3,
                k44 = 3e6,
                k45 = 7.2e2, k46 = 1e8, k47 = 84,
                k48 = 7.5e2, k49 = 7.4,
                k50 = 0.064, k51 = 1e6,
                k52 = 7.5e2, k53 = 49,
                k54 = 1e3,
                k55 = 1.6e4,
                k56 = 1e4,
                k57 = 3e6,
                k58 = 4e7,
                k61 = 4.628e-1, k62 = 8.038e7,
                k63 = 2.026e2, k64 = 3.377e7, k65 = 2.545e-1,
                k66 = 1.495e4,
                k67 = 2.012, k68 = 3.426e7,
                k69 = 4.4e-1,
                k70 = 1.381e-1,
                k71 = 1.25e-2,
                k72 = 2.947e-1,
                k73 = 2.369e1,
                k74 = 5.994,
                k75 = 1.317e8,
                k76 = 1.042e8,
                k77 = 2.93e-3,
                k78 = 210, k79 = 3.33e7,
                k80 = 1.707e6,
                k59 = 0.09, k60 = 0.47, K1 = 7.7e-8, K2 = 4.1e-7, K3 = 3e-7, K4 = 2.1e-6)

state      <- c(TF = 15e-12, FVII = 1e-8, TF_FVII = 0, TF_FVIIa = 0, FVIIa = 1e-10, FX = 1.6e-7, FXa = 0, FIIa = 0, TF_FVIIa_FX = 0,
                TF_FVIIa_FXa = 0, FIX = 9e-8, FIXa = 0, TF_FVIIa_FIX = 0, FII = 1.4e-6, FVIII = 7e-10, FVIIIa = 0, FIXa_FVIIIa = 0, 
                FIXa_FVIIIa_FX = 0, FVIIIa1L = 0, FVIIIa2 = 0, FV = 2e-8, FVa = 0, FXa_FVa = 0, FXa_FVa_FII = 0,
                mIIa = 0, TFPI = 2.5e-9, FXa_TFPI = 0, TF_FVIIa_FXa_TFPI = 0, AT = 3.4e-6, Fg = 9e-6, Fg_FIIa = 0, FnI = 0,
                FnI_FIIa = 0, FnII = 0, FnII_FIIa = 0, FnI2 = 0, FnI2_FIIa = 0, FnII2 = 0, Pn = 0, AP = 1e-6, tPA = 7e-11,
                PAI = 4e-10, Pg = 2e-6, FDP = 0, Tm = 1e-9, Tm_FIIa = 0, PC = 7e-8, Tm_FIIa_PC = 0, APC = 0, APC_FVa = 0,
                FVa5 = 0, FVa3 = 0, APC_FVa5 = 0, APC_FVa3 = 0, FVa53 = 0, LCA1 = 0, FXa_FVa5 = 0, FXa_FVa3 = 0,
                FXa_FVa5_FII = 0, FXa_FVa3_FII = 0, Tm_mIIa = 0, Tm_mIIa_PC = 0, FXa_FVa53 = 0, FXa_FVa53_FII = 0,
                APC_LCA1 = 0, Tm_FIIa_APC = 0, FII_FVa = 0)

times      <- seq(0, 200, by = 0.01)

out <- ode(y = state, times = times, func = mitraphanov, parms = parameters) %>% as.data.frame %>% as_tibble

#########
out %>% select(time, TF, FVIIa, FVII, TF_FVIIa) %>% gather(state,value,-time) %>% ggplot(aes(time, value)) + geom_line() + facet_grid(state ~ ., scales = "free")
