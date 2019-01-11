@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1
@#define EndoVariables = EndoVariables + [ "PI", "0", "theta^(1/(1-varepsilon))", "1" ]
// 1 is the growth rate of the stochastic trend of PI.
@#define EndoVariables = EndoVariables + [ "L", "0", "((1+vartheta)/psi)^(1/(1+vartheta))", "GN" ]
// GN is the growth rate of the stochastic trend of L.
// Whatever goes in that slot must be an MLV, with the usual _LEAD1, _LAG1, etc. defined for it
// This happens automatically if it is defined using the transformation engine, as here.
@#define PureTrendEndoVariables = PureTrendEndoVariables + [ "N", "GN" ]
// Defines N as a variable with stochastic trend N, whose growth rate is GN.
@#define EndoVariables = EndoVariables + [ "NU", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "AUX1", "0", "Inf", "1" ]
@#define ShockProcesses = ShockProcesses + [ "A", "0", "Inf", "A_STEADY", "rho_a", "sigma_a" ]
@#define ShockProcesses = ShockProcesses + [ "M", "0", "Inf", "1", "0", "-sigma_m" ]
@#define ShockProcesses = ShockProcesses + [ "beta", "0", "1", "beta_STEADY", "rho_b", "sigma_b" ]
@#define ShockProcesses = ShockProcesses + [ "Sg", "0", "1", "Sg_STEADY", "rho_g", "sigma_g" ]
@#define ShockProcesses = ShockProcesses + [ "GN", "0", "Inf", "1.00", "0.5", "0.01" ]
@#include "CreateShocks.mod"
@#include "ClassifyDeclare.mod"

parameters beta_STEADY, psi, A_STEADY, Sg_STEADY, PI_STEADY, vartheta, sigma, varepsilon, theta, phi_pi, phi_y, rho_a, rho_b, rho_g, sigma_g, sigma_b, sigma_a, sigma_m;

beta_STEADY = 0.994;
Sg_STEADY = 0.2;
A_STEADY = 1;
PI_STEADY = 1.005;
psi = 2;
vartheta = 1;
theta = 0.75;
varepsilon = 6;
sigma = 1;
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
    @#include "InsertNewModelEquations.mod"
    #Y = (A/NU) * L;
    #Y_LEAD = (A_LEAD/NU_LEAD) * L_LEAD;
    #G = Sg*Y;
    #G_LEAD = Sg_LEAD*Y_LEAD;
    #PI_STAR = (( 1 - theta * (PI^(varepsilon-1)) ) / (1 - theta))^(1/(1-varepsilon));
    #PI_STAR_LEAD = (( 1 - theta * (PI_LEAD^(varepsilon-1)) ) / (1 - theta))^(1/(1-varepsilon));
    #C = Y - G;
    #C_LEAD = Y_LEAD - G_LEAD;
    #W = psi * L^vartheta*C / ( 1 - psi * L^(1+vartheta) / ( 1+vartheta) ) * ( 1 - psi * STEADY_STATE(L)^(1+vartheta) / ( 1+vartheta) );
    #MC = W/A;
    #AUX2 = varepsilon / (varepsilon - 1) * AUX1;
    #AUX2_LEAD = varepsilon / (varepsilon - 1) * AUX1_LEAD;
    #R = ( ( PI_STEADY / beta_STEADY ) * ((PI/STEADY_STATE(PI))^phi_pi) * ((Y/STEADY_STATE(Y))^phi_y) * M );
    1 = R * beta_LEAD * ( C / C_LEAD ) / PI_LEAD;
    AUX1 = MC * (Y/C) + theta * beta_LEAD * PI_LEAD^(varepsilon) * AUX1_LEAD;
    AUX2 = PI_STAR * ((Y/C) + theta * beta_LEAD * ((PI_LEAD^(varepsilon-1))/PI_STAR_LEAD) * AUX2_LEAD);
    log( NU ) = log( theta * (PI^varepsilon) * NU_LAG + (1 - theta) * PI_STAR^(-varepsilon) );
    #error = ( R * beta_LEAD * ( C / C_LEAD ) / PI_LEAD - 1 ) / 1;
end;

steady_state_model;
    @#include "InsertNewStartSteadyStateEquations.mod"
    PI_ = PI_STEADY;
    PI_STAR_ = ( (1 - theta * (1 / PI_)^(1 - varepsilon) ) / (1 - theta) )^(1/(1 - varepsilon));
    NU_ = ( ( 1 - theta) / (1 - theta * PI_ ^varepsilon) ) * PI_STAR_ ^(-varepsilon);
    W_ = A_ * PI_STAR_ * ((varepsilon - 1) / varepsilon) * ( (1 - theta * beta_ * PI_ ^varepsilon)/(1 - theta * beta_ * PI_ ^(varepsilon-1)));
    C_ = (W_ /(psi * ((1/(1 - Sg_)) * NU_/ A_)^vartheta))^(1/(1 + vartheta));
    Y_ = (1 / (1 - Sg_)) * C_;
    G_ = Sg_ * Y_;
    L_ = Y_ *NU_ / A_;
    MC_ = W_ / A_;
    R_ = PI_ / beta_;
    AUX1_ = W_ / A_ * (Y_ /C_)/(1 - theta * beta_ * PI_ ^varepsilon);
    AUX2_ = PI_STAR_ * (Y_ /C_)/(1 - theta * beta_ * PI_ ^(varepsilon-1));
    @#include "InsertNewEndSteadyStateEquations.mod"
end;

shocks;
    @#include "InsertNewShockBlockLines.mod"
end;

steady;
check;

@#ifdef dynareOBC
    stoch_simul( order = 3, irf = 0, periods = 1100, irf_shocks = ( epsilon_beta ), replic = 100 ) error;
@#else
    stoch_simul( order = 3, irf = 0, periods = 1100, irf_shocks = ( epsilon_beta ), replic = 100 );
@#endif
