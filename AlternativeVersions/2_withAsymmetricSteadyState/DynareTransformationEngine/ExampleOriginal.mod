var Y, C, R, PI, W, L, PI_STAR, NU, MC, AUX1, AUX2, A, G, beta, M, Sg;
varexo epsilon_a, epsilon_m, epsilon_g, epsilon_b;

parameters beta_STEADY, A_STEADY, Sg_STEADY, PI_STEADY, varepsilon, theta, phi_pi, phi_y, rho_a, rho_b, rho_g, sigma_g, sigma_b, sigma_a, sigma_m, vartheta, psi;

beta_STEADY = 0.994;
Sg_STEADY = 0.2;
A_STEADY = 1;
PI_STEADY = 1.005;
psi = 2;
vartheta = 1;
theta = 0.75;
varepsilon = 6;
phi_pi = 1.5;
phi_y = 0.25;
rho_a = 0.9;
rho_g = 0.8;
rho_b = 0.8;
sigma_a = 0.0025;
sigma_m = 0.0025;
sigma_b = 0.26;
sigma_g = 0.0032;

model;
    1 = R * beta(+1) * ( C / C(+1) ) / PI(+1);
    #error = ( R * beta(+1) * ( C / C(+1) ) / PI(+1) - 1 ) / 1;
    W = psi * L^vartheta*C / ( 1 - psi * L^(1+vartheta) / ( 1+vartheta) ) * ( 1 - psi * STEADY_STATE(L)^(1+vartheta) / ( 1+vartheta) );
    MC = W/A;
    varepsilon * AUX1 = (varepsilon - 1) * AUX2;
    AUX1 = MC * (Y/C) + theta * beta(+1) * PI(+1)^(varepsilon) * AUX1(+1);
    AUX2 = PI_STAR * ((Y/C) + theta * beta(+1) * ((PI(+1)^(varepsilon-1))/PI_STAR(+1)) * AUX2(+1));
    R = ( ( PI_STEADY / beta_STEADY ) * ((PI/STEADY_STATE(PI))^phi_pi) * ((Y/STEADY_STATE(Y))^phi_y) * M );
    G = Sg*Y;
    1 = theta * (PI^(varepsilon-1)) + (1 - theta) * PI_STAR^(1 - varepsilon);
    NU = theta * (PI^varepsilon) * NU(-1) + (1 - theta) * PI_STAR^(-varepsilon);
    Y = C + G;
    Y = (A/NU) * L;
    
    A = A_STEADY^(1 - rho_a) * A(-1) ^(rho_a) * exp(sigma_a * epsilon_a);
    M = exp(-sigma_m * epsilon_m);
    log( beta / ( 1 - beta ) ) = (1 - rho_b) * log( beta_STEADY / ( 1 - beta_STEADY ) ) + rho_b * log( beta(-1) / ( 1 - beta(-1) ) ) + sigma_b*epsilon_b;
    log( Sg / ( 1 - Sg ) ) = (1 - rho_g) * log( Sg_STEADY / ( 1 - Sg_STEADY ) ) + rho_g * log( Sg(-1) / ( 1 - Sg(-1) ) ) + sigma_g*epsilon_g;
end;

steady_state_model;
    A = A_STEADY;
    Sg = Sg_STEADY;
    beta = beta_STEADY;
    M = 1;
    PI = PI_STEADY;
    PI_STAR = ( (1 - theta * (1 / PI)^(1 - varepsilon) ) / (1 - theta) )^(1/(1 - varepsilon));
    NU = ( ( 1 - theta) / (1 - theta * PI ^varepsilon) ) * PI_STAR ^(-varepsilon);
    W = A_STEADY * PI_STAR * ((varepsilon - 1) / varepsilon) * ( (1 - theta * beta * PI ^varepsilon)/(1 - theta * beta * PI ^(varepsilon-1)));
    C = (W /(psi * ((1/(1 - Sg)) * NU/ A_STEADY)^vartheta))^(1/(1 + vartheta));
    Y = (1 / (1 - Sg)) * C;
    G = Sg * Y;
    L = Y *NU / A_STEADY;
    MC = W / A_STEADY;
    R = PI / beta;
    AUX1 = W / A_STEADY * (Y /C)/(1 - theta * beta * PI ^varepsilon);
    AUX2 = PI_STAR * (Y /C)/(1 - theta * beta * PI ^(varepsilon-1));
end;

shocks;
    var epsilon_a = 1;
    var epsilon_g = 1;
    var epsilon_b = 1;
    var epsilon_m = 1;
end;

steady;
check;

@#ifdef dynareOBC
    stoch_simul( order = 3, irf = 0, periods = 1100, irf_shocks = ( epsilon_b ), replic = 100 ) error;
@#else
    stoch_simul( order = 3, irf = 0, periods = 1100, irf_shocks = ( epsilon_b ), replic = 100 );
@#endif
