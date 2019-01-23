function E_by_F_1_ = solve_ss(SpatialPointsPerDimension,param)

global par E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_

SpatialDimensions = 2;
SpatialNorm = 2;


alpha = 0.3;
gamma = 0.5;
kappa = 0.5;
nu = 2;
varsigma = 1.5;
zeta = param.zeta;
lambda = param.lambda;
deltaJ = 0.01;
deltaK = 0.03;
Phi2 = param.Phi2;
PhiL = param.PhiL;

thetaC = 4;
thetaF = 1;
thetaL = 0.25 / 0.75 * thetaF * gamma;
thetaH = 4;
thetaN = param.thetaN;
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
Omega = param.Omega;


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

if isempty(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_)
    E_by_F_1_ = GetE_1_by_F_1_homotopy( 1 , SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega, ...
                   lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, ...
                   GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1, ...
                   psi2, psi3, tau_, dBar, beta_, varsigma );
else
    E_by_F_1_ = GetE_1_by_F_1( 1 , SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega, ...
                   lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, ...
                   GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1, ...
                   psi2, psi3, tau_, dBar, beta_, varsigma );
end
           
end
           
           
           
           
           