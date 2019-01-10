clear; clear global;

alpha = 0.3;
gamma = 0.5;
kappa = 0.5;
nu = 2;
varsigma = 1.5;
zeta = 8;
lambda = 0.1;
deltaJ = 0.01;
deltaK = 0.03;
Phi2 = 4;

thetaC = 4;
thetaF = 1; % food off-premises, food services + clothing = about 20% of ( PCE minus housing ) https://www.bea.gov/iTable/iTable.cfm?reqid=19&step=2#reqid=19&step=3&isuri=1&1910=x&0=-9&1921=survey&1903=65&1904=2015&1905=2017&1906=a&1911=0
thetaL = 0.25 / 0.75 * thetaF * gamma; % target of 0.75 for steady-state land use in agriculture, following data from https://www.ers.usda.gov/data-products/major-land-uses/
thetaH = 4;
thetaN = 10;
psi1 = 0.5;
psi2 = 0.5;
psi3 = psi1 * 0.02 / ( 1 - 0.02 ); % http://eyeonhousing.org/2013/01/latest-study-shows-average-buyer-expected-to-stay-in-a-home-13-years/

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
Omega = 3; % pop/km^2 for the contiguous US is 41.5. for wyoming it is 2.33 for new jersey it is 470. correspond to abs log ratios of 2.88 and 2.43 respectively.

dBar = ( 2 ) ^ ( 1 / ( 2 ) );

GA_  = 1.005;
GN_  = 1.0025;
tau_  = 1;
phi_  = 1;
beta_  = 0.99;

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

SpatialPointsPerDimension = 3;
shape = 'T';

SpatialPoints = SpatialPointsPerDimension*SpatialPointsPerDimension;

E_by_F_ = zeros( SpatialPoints,1 );
N_x = zeros( SpatialPoints,1 );
F_x = zeros( SpatialPoints,1 );
K_x = zeros( SpatialPoints,1 );
H_x = zeros( SpatialPoints,1 );
Q_x = zeros( SpatialPoints,1 );
SN_xx = zeros( SpatialPoints,SpatialPoints );

% solver fmincon || fmincon_global || fsolve
solver = 'fmincon';

E_by_F_(1) = GetE_1_by_F_1_new( solver, shape, 1 , SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega,...
    lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_,...
    GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_ , psi1,...
            psi2, psi3, tau_, dBar, beta_, varsigma);
% E_by_F_(1) = GetE_1_by_F_1_alt(  1 , SpatialPointsPerDimension  , psi3 , psi1 , GN_);
        
% for x=1:SpatialPoints
% E_by_F_(x) = GetE_x_by_F_x( x );
% N_x(x) = GetN_x( x );
% F_x(x) = GetF_x( x );
% K_x(x) = GetK_x( x );
% H_x(x) = GetH_x( x );
% Q_x(x) = GetQ_x( x );
% SN_xx(: , x) = GetSN_xx( x );
% end

%     A_x_ = 1;
%     N_x_LAG_ = N_x_ ./ GN_; 
%     P_Over_Q_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GSRKTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;
%     
%     util_SN_x_ = zeros( SpatialPoints,1 );
%     SD_x_ = zeros( SpatialPoints,1 );
%     SN_x_out = zeros( SpatialPoints,1 );
%     SN_x_in = zeros( SpatialPoints,1 );
%     integral_N = zeros( SpatialPoints,1 );
%     for x=1:SpatialPoints 
%         integral_util_SN = zeros( SpatialPoints,1 );
%         integral_SD = zeros( SpatialPoints,1 );
%         integral_SN_out = zeros( SpatialPoints,1 );
%         integral_SN_in = zeros( SpatialPoints,1 );
%         for tildex=1:SpatialPoints
%             integral_util_SN( tildex ) = Weight * GN_ * N_x_LAG_( tildex ) * log( SNxx_( x , tildex ) / N_x_LAG_( x ) );
%             integral_SD( tildex ) = Weight * d( x , tildex ) * SNxx_( x , tildex  );
%             integral_SN_out( tildex ) = Weight * ( SNxx_( x , tildex ) );
%             integral_SN_in( tildex ) = Weight * ( SNxx_( tildex , x ) );
%         end
%         util_SN_x_( x ) = sum( integral_util_SN );
%         SD_x_( x ) = sum( integral_SD );
%         SN_x_out( x ) = sum( integral_SN_out );
%         SN_x_in( x ) = sum( integral_SN_in );
%         integral_N( x ) = Weight * ( N_x_( x ) );
%     end 
%     
%     SN_x_ = SN_x_out;
%     E_x_ = E_by_F_ .* F_x_;    
%     L_x_ = thetaF .* gamma ./ ( thetaL .* E_by_F_ + thetaF .* gamma );
%     ZF_x_ = ( F_x_ ./ L_x_ .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) );
%     SP_x_ = ( 1 - gamma ) .* F_x_ ./ ZF_x_;
%     P_x_ = P_Over_Q_x_ .* Q_x_;
%     Z_x_ = ( ( ( K_x_ ./ GYTrend_ ) .^ alpha .* ( A_ .* H_x_ ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* SP_x_ ./ P_x_ ) .^ kappa ) .^ ( 1 / ( 1 - kappa ) );
%     SRK_x_ = ( 1 - kappa ) .* alpha .* SP_x_ .* Z_x_ ./ ( K_x_ ./ GYTrend_ );
%     W_x_ = ( 1 - kappa ) .* ( 1 - alpha ) .* SP_x_ .* Z_x_ ./ H_x_;
%     C_x_ = thetaC .* F_x_ ./ ( thetaF .* P_x_ );
%     I_x_ = K_x_ .* ( 1 - ( 1 - deltaK ) ./ GYTrend_ ) ./ ( 1 - Phi2 ./ 2 .* ( GYTrend_ - 1 ) .^ 2 );
%     M_x_ = kappa .* SP_x_ .* Z_x_ ./ P_x_;
%     Y_x_ = C_x_ + I_x_ + M_x_;
%     
%     YBar_x_ = zeros( SpatialPoints,1 );
%     for x=1:SpatialPoints 
%         integral_YBar = zeros( 1,SpatialPoints );
%         for tildex=1:SpatialPoints
%             integral_YBar( tildex ) = Weight * Y_x_( tildex ) * ( P_x_( tildex ) )^((1+lambda)/lambda) * exp( -tau_ * d( x , tildex ) / lambda );
%         end
%         YBar_x_( x ) = sum( integral_YBar );
%     end
%     
%     J_x_ = ( Z_x_ - ZF_x_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_x_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_x_ ); 
% 
%     for x=1:SpatialPoints 
%         integral_Px = zeros( SpatialPoints,1 );
%         for tildex=1:SpatialPoints
%             integral_Px( tildex ) = Weight * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
%         end
%     end
%         
%     U_x_ = ( C_x_ ./ N_x_LAG_ ).^thetaC .* ( E_x_ ./ N_x_LAG_ ).^thetaF ...
%         .* ( ( 1 - L_x_ )./N_x_LAG_ ).^thetaL ...
%         .* ( ( 1 / ( 1 + nu ) ) * Gamma^( 1 + nu ) - ( 1 / ( 1 + nu ) ) .* ( H_x_ ./ N_x_LAG_ ).^( 1 + nu ) ).^thetaH...
%         .* ( .5 * Omega ^ 2 - .5 * ( log( GN_ .* N_x_LAG_ ) ).^2 ).^ thetaN ...
%          .* ( 1 - SN_x_ ./ N_x_LAG_ ).^psi1 ...
%          .* ( dBar - SD_x_ ./ SN_x_ ).^ psi2 ...
%          .* exp( psi3 .* util_SN_x_ );
%     
%     muN = beta_ .* ( ( U_x_ .* GUTrend_ ) .^ ( 1 - varsigma ) ) .* ( 1 + ( 1 - varsigma ) .* ( ...
%             thetaH .* ( H_x_ ./ N_x_LAG_ ) .^ ( 1 + nu ) ./ ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x_ ./ N_x_LAG_ ) .^ ( 1 + nu ) ) ...
%             - thetaN .* log( N_x_ ) ./ ( 1 ./ 2 .* Omega .^ 2 - 1 ./ 2 .* log( N_x_ ) .^ 2 ) ...
%             + psi1 .* GN_ .* SN_x_ ./ ( N_x_ - GN_ .* SN_x_ ) ...
%             - ( thetaC + thetaF + thetaL + psi3 ) ...
%         ) ) ./ ( 1 - beta_ .* GmuNTrend_ .* GN_ );

%
% @#for Point1 in 1 : SpatialNumPoints
%     @#define Index1 = IndicesStringArray[Point1]
%     @#for Point2 in 1 : SpatialNumPoints
%         @#define Index2= IndicesStringArray[Point2]
%             Distance@{Index1}@{Index2}_ = getDistance( Point1 , Point2 , SpatialPointsPerDimension )
%     @#endfor
% @#endfor
% 
% @#for Point1 in 1 : SpatialNumPoints
%     @#define Index1 = IndicesStringArray[Point1]
%     E_by_F_@{Index1}_ = GetE_x_by_F_x_( Point1 );
%     N_@{Index1}_ = GetN_x_( Point1 );
%     F_@{Index1}_ = GetF_x_( Point1 );
%     K_@{Index1}_ = GetK_x_( Point1 );
%     H_@{Index1}_ = GetH_x_( Point1 );
%     Q_@{Index1}_ = GetQ_x_( Point1 );
%     SN_@{Index1}_ = GetSN_x_( Point1 , SpatialPointsPerDimension );
%     SD_@{Index1}_ = GetSD_x_( Point1 , SpatialPointsPerDimension );
%     A_@{Index1}_ = 1;
%     A_@{Index1}_LEAD_ = GA_;
%     N_LAG_@{Index1}_ = N_@{Index1}_ ./ GN_;
%     P_Over_Q_@{Index1}_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GQTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;
%     E_@{Index1}_  = E_by_F_@{Index1}_ * F_@{Index1}_;
%     L_@{Index1}_  = thetaF * gamma / ( thetaL * E_by_F_@{Index1}_  + thetaF * gamma );
%     ZF_@{Index1}_ = ( F_@{Index1}_ / L_@{Index1}_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
%     SP_@{Index1}_ = ( 1 - gamma ) * F_@{Index1}_ / ZF_@{Index1}_;
%     P_@{Index1}_ = P_Over_Q_@{Index1}_ * Q_@{Index1}_;
%     Z_@{Index1}_ = ( ( ( K_@{Index1}_ / GYTrend_ ) ^ alpha * ( A_@{Index1}_ * H_@{Index1}_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_@{Index1}_ / P_@{Index1}_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
%     SRK_@{Index1}_ = ( 1 - kappa ) * alpha * SP_@{Index1}_ * Z_@{Index1}_ / ( K_@{Index1}_ / GYTrend_ );
%     W_@{Index1}_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_@{Index1}_ * Z_@{Index1}_ / H_@{Index1}_;
%     C_@{Index1}_ = thetaC * F_@{Index1}_ / ( thetaF * P_@{Index1}_ );
%     I_@{Index1}_ = K_@{Index1}_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
%     M_@{Index1}_ = kappa * SP_@{Index1}_ * Z_@{Index1}_ / P_@{Index1}_;
%     Y_@{Index1}_ = C_@{Index1}_ + I_@{Index1}_ + M_@{Index1}_;
%     YBar_@{Index1}_ = 0
%     @#for Point2 in 1 : SpatialNumPoints
%         @#define Index2= IndicesStringArray[Point2]
%         + (1 / (SpatialPointsPerDimension^2)) * Y_@{Index2}_ * P_@{Index2}_ ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance@{Index1}@{Index2}_ )
%     @#endfor
%     ;
%     @#for Point2 in 1 : SpatialNumPoints
%         @#define Index2 = IndicesStringArray[Point2]
%         SN@{Index1}@{Index2}_ = GetSNxx_( Point1 , Point2 );
%     @#endfor
%     J_@{Index1}_ = ( Z_@{Index1}_ - ZF_@{Index1}_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_@{Index1}_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_@{Index1}_ );
%     U@{Index1}_ = ( C@{Index1}_ / N@{Index1}_LAG_ ) ^ thetaC
%     * ( E@{Index1}_ / N_LAG_@{Index1}_ ) ^ thetaF
%     * ( ( 1 - L@{Index1}_ ) / N_LAG_@{Index1}_ ) ^ thetaL
%     * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_ / N_LAG_@{Index1}_ ) ^ ( 1 + nu ) ) ^ thetaH
%     * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N_LAG_@{Index1}_ / N_LAG_ ) ^ 2 ) ^ thetaN
%     * ( 1 - SN@{Index1}_ / N_LAG_@{Index1}_ ) ^ psi1
%     * ( dBar - SD@{Index1}_ / SN@{Index1}_ ) ^ psi2
%     * exp( psi3 * ( 0
%     @#for Point2 in 1 : SpatialNumPoints
%         @#define Index2 = IndicesStringArray[Point2]
%         + (1 / (SpatialPointsPerDimension^2)) * N_LAG_@{Index1}_ / N_LAG_ * log( SN@{Index1}@{Index2}_ / N_LAG_@{Index1}_ )
%     @#endfor
%     ) );
%     U_@{Index1}_LEAD_ = U_@{Index1}_ * GUTrend_;
%     H_@{Index1}_LEAD_ = H_@{Index1}_ * GN_;
%     muN_@{Index1}_ = beta .* ( ( U_@{Index1}_ .* GUTrend_ ) ^ ( 1 - varsigma ) ) .* ( 1 + ( 1 - varsigma ) .* ( ...
%             thetaH .* ( H_@{Index1}_ ./ N_LAG_@{Index1}_ ) .^ ( 1 + nu ) ./ ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_@{Index1}_ ./ N_LAG_@{Index1}_ ) .^ ( 1 + nu ) ) ...
%             - thetaN .* log( N_@{Index1}_ ) ./ ( 1 ./ 2 .* Omega .^ 2 - 1 ./ 2 .* log( N_@{Index1}_ ) .^ 2 ) ...
%             + psi1 .* GN_ .* SN_@{Index1}_ ./ ( N_@{Index1}_ - GN_ .* SN_@{Index1}_ ) ...
%             - ( thetaC + thetaF + thetaL + psi3 ) ...
%         ) ) ./ ( 1 - beta .* GmuNTrend_ .* GN_ );
% @#endfor
% 
% @#for VariableName in AggregatedVariables
%     @{VariableName}_ = 0
%     @#for Point1 in 1 : SpatialNumPoints
%         @#define Index1 = IndicesStringArray[Point1]
%         +  (1 / (SpatialPointsPerDimension^2)) * @{VariableName}@{Index1}_
%     @#endfor
%     ;
% @#endfor
% 
%  @#include "InsertNewEndSteadyStateEquations.mod"
