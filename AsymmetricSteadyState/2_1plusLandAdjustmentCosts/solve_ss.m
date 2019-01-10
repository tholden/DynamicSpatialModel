function E_by_F_1_ = solve_ss(SpatialPointsPerDimension)

global par

SpatialDimensions = 2;
SpatialNorm = 2;

param_thetaN = 5; %10
param_Omega = 2.5; %3
param_lambda = 0.1; %0.1

param_Phi2 = 4; %4
param_PhiL = 2; %2
param_zeta = 8; %8

alpha = 0.3;
gamma = 0.5;
kappa = 0.5;
nu = 2;
varsigma = 1.5;
zeta = param_zeta;
lambda = param_lambda;
deltaJ = 0.01;
deltaK = 0.03;
Phi2 = param_Phi2;
PhiL = param_PhiL;

thetaC = 4;
thetaF = 1;
thetaL = 0.25 / 0.75 * thetaF * gamma;
thetaH = 4;
thetaN = param_thetaN;
psi1 = 0.5;
psi2 = 0.5;
psi3 = psi1 * 0.02 / ( 1 - 0.02 );

UtilityParamSum = thetaC + thetaF + thetaL + thetaH + thetaN + psi1 + psi2 + psi3;

thetaC = thetaC / UtilityParamSum;
thetaF = thetaF / UtilityParamSum;
thetaL = thetaL / UtilityParamSum;
thetaH = thetaH / UtilityParamSum;
thetaN = thetaN / UtilityParamSum;
psi1 = psi1 / UtilityParamSum;
psi2 = psi2 / UtilityParamSum;
psi3 = psi3 / UtilityParamSum;

Gamma = 1;
Omega = param_Omega;


dBar = ( SpatialDimensions * ( 0.5 ^ ( SpatialNorm ) ) ) ^ ( 1 / ( SpatialNorm ) );

GA_ = 1.005;
GN_ = 1.0025;
beta_ = 0.99;
phi_ = 1;
tau_ = 1;

GYTrend_    = ( GA_ * GN_ ) ^ ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) / ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) - lambda ) );
GZTrend_    = GYTrend_ ^ ( 1 / ( 1 + lambda ) );
GJTrend_    = GYTrend_ ^ ( 1 / ( 1 + lambda ) );
GFTrend_    = GYTrend_ ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
GPTrend_    = GYTrend_ ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
GQTrend_    = GYTrend_ ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
GSRKTrend_  = GYTrend_ ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
GSPTrend_   = GYTrend_ ^ ( - gamma / ( 1 + lambda ) );
GUTrend_    = GYTrend_ ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN_ ^ ( thetaC + thetaF + thetaL );
GmuNTrend_  = GYTrend_ ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN_ ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );

Xi_LEAD_ = beta_ * GN_ / GFTrend_ * GUTrend_ ^ ( 1 - varsigma );
R_ = 1 / Xi_LEAD_;
N_ = 1;
N_LAG_ = 1/GN_;

par.thetaC = thetaC;
par.thetaF = thetaF;
par.Phi2 = Phi2;
par.GYTrend_ = GYTrend_;
par.Xi_LEAD_ = Xi_LEAD_;
par.GSRKTrend_ = GSRKTrend_;

E_by_F_1_ = GetE_1_by_F_1_homotopy( 1 , SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega, ...
               lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, ...
               GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1, ...
               psi2, psi3, tau_, dBar, beta_, varsigma );
           
end
           
           
           
           
           