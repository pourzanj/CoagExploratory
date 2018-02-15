functions {
  
  // coag system
  real[] dz_dt(real t, real[] z, real[] theta, real[] x_r, int[] x_i) {
    
    real FII = z[1];
    real FIIa = z[2];
    real AT = z[3];
    real Fg = z[4];
    real Fn = z[5];
    real tPA = z[6];
    real PAI = z[7];

    real g_FII = theta[1];
    real g_AT = theta[2];
    real c_cascade = theta[3];
    real b_cascade = theta[4];
    real c_tfpi = theta[5];
    real b_tfpi = theta[6];
    
    real k_FIIa = theta[7];
    real c_FIIa = theta[7];
    real K_FIIa = theta[8];
    
    real k_AT = theta[10];
    
    real k_clot = theta[11];
    real c_clot = theta[12];
    real K_clot = theta[13];
    
    real k_lys = theta[14];
    real c_lys = theta[15];
    real K_lys = theta[16];
    
    real k_PAI = theta[17];
    
    real k_prothrombinase = theta[18];
    real prothrombinase_max = theta[19];
    real c_prothrombinase = theta[20];
    real b_prothrombinase = theta[21];
    
    real cascade_delay = 1/(1+exp(c_cascade-b_cascade*t));
    real tfpi = 1/(1+exp(-c_tfpi+b_tfpi*t));
    real prothrombinase = prothrombinase_max/(1+exp(c_prothrombinase-b_prothrombinase*t));
    
    real r_FIIa = k_FIIa*cascade_delay*tfpi*((FII/g_FII)^c_FIIa)/(K_FIIa + (FII/g_FII)^c_FIIa);
    real r_prothrombinase = k_prothrombinase*prothrombinase*(FII/g_FII);
    real r_AT = k_AT*FIIa*(AT/g_AT);
    real r_clot = k_clot*FIIa*((Fg/9e-6)^c_clot)/(K_clot + (Fg/9e-6)^c_clot);
    real r_lys = k_lys*tPA*(Fn^c_lys)/(K_lys + Fn^c_lys);
    real r_PAI = k_PAI*tPA*PAI;

    real dFII = g_FII*(-r_FIIa-r_prothrombinase);
    real dFIIa = r_FIIa -r_AT;
    real dAT = -g_AT*r_AT;
    real dFg = -r_clot;
    real dFn = 1e7*(r_clot - r_lys);
    real dtPA = -r_PAI;
    real dPAI = -r_PAI;

    return { dFII, dFIIa, dAT, dFg, dFn, dtPA, dPAI };
  }
  
}

data {
  int<lower = 0> Nt; // total num time points to capture ODE
  real ts[Nt];       // time points to capture ODE
  
  real y0[7];
  real y[Nt,7];       // time series from mitraphanov
  
  real sigma[7];
}

parameters { 

    real<lower=0> c_cascade;
    real<lower=0> b_cascade;
    real<lower=0> c_tfpi;
    real<lower=0> b_tfpi;
    
    real<lower=0> k_FIIa;
    real<lower=0.95> c_FIIa;
    real<lower=0> K_FIIa;
    
    real<lower=0> k_clot;
    real<lower=0.95> c_clot;
    real<lower=0> K_clot;
    
    real<lower=0> k_lys;
    real<lower=0.95> c_lys;
    real<lower=0> K_lys;
    
    real<lower=0> c_prothrombinase;
    real<lower=0> b_prothrombinase; 
    
    //real<lower=0> sigma[7];
}

transformed parameters {
  real g_FII = 100/1.4e-6;
  real g_AT = 100/3.4e-6;
  real k_AT = 1.6e4;
  real k_PAI = 4.5e5;
  real k_prothrombinase = 1e8;
  real prothrombinase_max = 1.5e-10;
  
  real theta[21] = {g_FII, g_AT, c_cascade, b_cascade, c_tfpi, b_tfpi, k_FIIa, c_FIIa, K_FIIa, k_AT, k_clot, c_clot, K_clot, k_lys, c_lys, K_lys, k_PAI,k_prothrombinase,prothrombinase_max, c_prothrombinase, b_prothrombinase}; 
  
  // all 6 soln curves at the T time points. This is just a temp. variable (see below)
  real z0[7] = {1e2, 0, 1e2, 9e-6, 0, 7e-11, 4e-10};
  real z[Nt,7]; 

  // get TEG readings and ending FDP value for each patient
  z = integrate_ode_bdf(dz_dt, z0, 0, ts, theta, rep_array(0.0, 0), rep_array(0, 0), 1e-6, 1e-5, 1e3);
}

model {

  c_cascade ~ normal(5, 5);
  b_cascade ~ normal(4e-2, 5e-1);
  c_tfpi ~ normal(10,10);
  b_tfpi ~ normal(1.5e-2,1e-3);

  k_FIIa ~ normal(3.5e-9,1e-8);
  c_FIIa ~ normal(1,5e-1);
  K_FIIa ~ normal(1.4e-6,1e-5);

  k_clot ~ normal(3,5e-0);
  c_clot ~ normal(1,5e-2);
  K_clot ~ normal(0.75,3e-0);

  k_lys ~ normal(1,1e-1);
  c_lys ~ normal(1,5e-1);
  K_lys ~ normal(0.5,2e-1);
  
  c_prothrombinase ~ normal(7,7);
  b_prothrombinase ~ normal(6e-3,6e-3);

  for (k in 1:7) {
    //y[,k] ~ lognormal(log(to_vector(z[, k])+ 1e-12), sigma[k]);
    y[,k] ~ normal(z[, k], sigma[k]);
  }
}
