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
    r50B <- k69*APC_FVa53
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
    dTF_FVIIa <- r2 - r6 - r7 - r8A + r8B - r22 - r27
    dFVIIa <- -r2 + r3 + r4 + r5
    dFX <- -r6A - r12A + r14- r28
    dFXa <- -r7 + r12 - r17 - r20 - r23 + r28 - r55 - r56 + r61 + r62 - r66 + r69 + r70
    dFIIa <- r9 + r19 - r26 - r30A + r30B - r31A + r31B - r33A + r33B - r34 - r43 + r59 + r60 + r68
    dTF_FVIIa_FX <- r6A - r6B
    dTF_FVIIa_FXa <- r6B + r7 - r21
    dFIX <- -r8A
    dFIXa <- r8B - r11 + r14 + r15 - r25
    dTF_FVIIa_FIX <- r8
    dFII <- -r9 - r18 - r57 -r58 + r62 - r67 + r70 - r71
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
    dmIIa <- r18B - r19 - r24 + r57 + r58 - r59 - r60 - r63 - r68
    dTFPI <- -r20 - r21
    dFXa_TFPI <- r20 - r22
    dTF_FVIIa_FXa_TFPI <- r21 + r22
    dAT <- -r23 - r24 - r25 -r26 - r27 - r35 - r36 - r37 - r45 - r65
    dFg <- -r30A
    dFg_FIIa <- r30A - r30B
    dFnI <- r30 - r31 -2*r32
    dFnI_FIIa <- r31A - r31B - r36
    dFnII <- r31 + -r34 - r41
    dFnII_FIIa <- r34 - r37
    dFnI2 <- r32 - r33A
    dFnI2_FIIa <- r33 - r35
    dFnII2 <- r33B - r42
    dPn <- -r38 + r40
    dAP <- -r38
    dtPA <- -r39
    dPAI <- -r39
    dPg <- -r40
    dFDP <- r41 + 2*r42
    dTm <- -r43 + r45 - r63 + r65
    dTm_FIIa <- r43 - r44A + r44B - r45 - r54
    dPC <- -r44 - r64
    dTm_FIIa_PC <- r44A - r44B
    dAPC <- r44 - r46 + r47 + r48 - r49A + r49B - r50A + r50B - r53 - r54 + r64
    dAPC_FVa <- r46 - r47 - r48
    dFVa5 <- r47 - r49 - r55
    dFVa3 <- r48 - r50 - r51 - r56
    dAPC_FVa5 <- r49A - r49B
    dAPC_FVa3 <- r50A - r50B
    dFVa53 <- r49 + r50 - r52 - r66
    dLCA1 <- r51 + r52 - r53 + r61 + r62 + r69 + r70
    dFXa_FVa5 <- r55 - r57A + r57B - r72
    dFXa_FVa3 <- r56 - r58A + r58B - r61
    dFXa_FVa5_FII <- r57A - r57B
    dFXa_FVa3_FII <- r58A - r58B - r62
    dTm_mIIa <- r63 - r64A + r64B - r65
    dTm_mIIa_PC <- r64A - r64B
    dFXa_FVa53 <- r66 - r67A + r67B - r69 + r72
    dFXa_FVa53_FII <- r67A - r67B - r70
    
    list(c(dTF, dFVII, dTF_FVII, dTF_FVIIa, dFVIIa, dFX, dFXa, dFIIa, dTF_FVIIa_FX,
           dTF_FVIIa_FXa, dFIX, dFIXa, dTF_FVIIa_FIX, dFII, dFVIII, dFVIIIa, dFIXa_FVIIIa,
           dFIXa_FVIIIa_FX, dFVIIIa1L, dFVIIIa2, dFV, dFVa, dFXa_FVa, dFXa_FVa_FII,
           dmIIa, dTFPI, dFXa_TFPI, dTF_FVIIa_FXa_TFPI, dAT, dFg, dFg_FIIa, dFnI,
           dFnI_FIIa, dFnII, dFnII_FIIa, dFnI2, dFnI2_FIIa, dFnII2, dPn, dAP, dtPA,
           dPAI, dPg, dFDP, dTM, dTm_FIIa, dPC, dTm_FIIa_PC, dAPC, dAPC_FVa,
           dFVa5, dFVa3, dAPC_FVa5, dAPC_FVa3, dFVa53, dLCA1, dFXa_FVa5, dFXa_FVa3,
           dFXa_FVa5_FII, dFXa_FVa3_FII, dTm_mIIa, dTm_mIIa_PC, dFXa_FVa53, dFXa_FVa53_FII))
  })
}



parameters <- c(k41 = 1, k46, k45, k47, k46, k48, k49,
                k52, k53, k54, k55, k56,
                k61, k62, k63, k64, k65, k66)
state      <- c(X = 1, Y = 1, Z = 1)
times      <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

#########
