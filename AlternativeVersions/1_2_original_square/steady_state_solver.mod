clear;


@#includepath "DynareTransformationEngine"

@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1

@#define SpatialDimensions = 2
@#define SpatialPointsPerDimension = 3
@#define SpatialShape = "Plane"

@#if SpatialDimensions == 1
    @#define SpatialNorm = "1"
@#else
    @#define SpatialNorm = "2"
@#endif

@#if SpatialShape[1] == "P"
    @#define SpatialShockProcesses = SpatialShockProcesses + [ "AT", "0", "Inf", "1", "0.9", "0.001", "exp(-zeta*@)#" ]
@#else
    @#define SpatialShockProcesses = SpatialShockProcesses + [ "AT", "0", "Inf", "1", "0.9", "0.001", "(exp(-zeta*@+zeta*dBar)+exp(zeta*@-zeta*dBar))/(exp(zeta*dBar)+exp(-zeta*dBar))#" ]
@#endif

@#define ShockProcesses = ShockProcesses + [ "GA", "0", "Inf", "1.005", "0.8", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "GN", "0", "Inf", "1.0025", "0.5", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "tau", "0", "Inf", "1", "0.95", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "phi", "0", "Inf", "1", "0.95", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "beta", "0", "1", "0.99", "0.95", "0.001" ]

@#define PureTrendEndoVariables = PureTrendEndoVariables + [ "AP", "GA" ]

@#include "CreateShocks.mod"

@#define EndoVariables = EndoVariables + [ "R", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GYTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GZTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GFTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GPTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GUTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GmuNTrend", "0", "Inf", "1" ]

@#for Point in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point]
    @#define EndoVariables = EndoVariables + [ "C" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "K" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "I" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "E" + Index1, "0", "Inf", "GFTrend" ]
    @#define EndoVariables = EndoVariables + [ "F" + Index1, "0", "Inf", "GFTrend" ]
    @#define EndoVariables = EndoVariables + [ "Q" + Index1, "0", "Inf", "GPTrend" ]
    @#define EndoVariables = EndoVariables + [ "J" + Index1, "0", "Inf", "GZTrend" ]
    @#define EndoVariables = EndoVariables + [ "L" + Index1, "0", "Inf", "1" ]
    @#define EndoVariables = EndoVariables + [ "H" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "N" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "SN" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "SD" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "muN" + Index1, "0", "Inf", "GmuNTrend" ]
    @#define EndoVariables = EndoVariables + [ "U" + Index1, "0", "Inf", "GUTrend" ]
@#endfor

@#define EndoVariables = EndoVariables + [ "C", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "K", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "I", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "F", "0", "Inf", "GFTrend" ]
@#define EndoVariables = EndoVariables + [ "J", "0", "Inf", "GZTrend" ]
@#define EndoVariables = EndoVariables + [ "L", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "H", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "SN", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "SD", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "U", "0", "Inf", "GUTrend" ]

@#define AggregatedVariables = [ "C", "K", "I", "F", "J", "L", "H", "SN", "SD", "U" ]

@#include "ClassifyDeclare.mod"

parameters alpha gamma kappa nu varsigma zeta lambda deltaJ deltaK Phi2 thetaC thetaF thetaL thetaH thetaN psi1 psi2 psi3 Gamma Omega dBar;
parameters UtilityParamSum;

load('param_vals.mat');

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
thetaF = 1; // food off-premises, food services + clothing = about 20% of ( PCE minus housing ) https://www.bea.gov/iTable/iTable.cfm?reqid=19&step=2#reqid=19&step=3&isuri=1&1910=x&0=-9&1921=survey&1903=65&1904=2015&1905=2017&1906=a&1911=0
thetaL = 0.25 / 0.75 * thetaF * gamma; // target of 0.75 for steady-state land use in agriculture, following data from https://www.ers.usda.gov/data-products/major-land-uses/
thetaH = 4;
thetaN = param_thetaN;
psi1 = 0.5;
psi2 = 0.5;
psi3 = psi1 * 0.02 / ( 1 - 0.02 ); // http://eyeonhousing.org/2013/01/latest-study-shows-average-buyer-expected-to-stay-in-a-home-13-years/

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
Omega = 3; // pop/km^2 for the contiguous US is 41.5. for wyoming it is 2.33 for new jersey it is 470. correspond to abs log ratios of 2.88 and 2.43 respectively.

@#if SpatialShape[1] == "P"
    dBar = ( @{SpatialDimensions} ) ^ ( 1 / ( @{SpatialNorm} ) );
@#else
    dBar = ( @{SpatialDimensions} * ( 0.5 ^ ( @{SpatialNorm} ) ) ) ^ ( 1 / ( @{SpatialNorm} ) );
@#endif


@#include "InsertNewStartSteadyStateEquations.mod"

   GYTrend_    = ( GA_ * GN_ ) ^ ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) / ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) - lambda ) );
   GZTrend_    = GYTrend_ ^ ( 1 / ( 1 + lambda ) );
   GJTrend_    = GYTrend_ ^ ( 1 / ( 1 + lambda ) );
   GFTrend_    = GYTrend_ ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
// GSRLTrend_  = GYTrend_ ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
// GWTrend_    = GYTrend_ ^ ( ( 1 - gamma ) / ( 1 + lambda ) ) / GN_;
   GPTrend_    = GYTrend_ ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
   GQTrend_    = GYTrend_ ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
   GSRKTrend_  = GYTrend_ ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
   GSPTrend_   = GYTrend_ ^ ( - gamma / ( 1 + lambda ) );
// GPiTrend_   = GYTrend_ ^ ( - gamma / ( 1 + lambda ) );
// GYBarTrend_ = GYTrend_ ^ ( - gamma / lambda );
   GUTrend_    = GYTrend_ ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN_ ^ ( thetaC + thetaF + thetaL );
   GmuNTrend_  = GYTrend_ ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN_ ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );

Xi_LEAD_ = beta_ * GN_ / GFTrend_ * GUTrend_ ^ ( 1 - varsigma );
R_ = 1 / Xi_LEAD_;
N_ = 1;
N_LAG_ = 1/GN_;

E_by_F_1_ = GetE_1_by_F_1( 1 , SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma,...
    lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_,...
    GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_ );


@#for Point1 in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point1]
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2= IndicesStringArray[Point2]
            Distance@{Index1}@{Index2}_ = getDistance( Point1 , Point2 , SpatialPointsPerDimension )
    @#endfor
@#endfor

@#for Point1 in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point1]
    E_by_F_@{Index1}_ = GetE_x_by_F_x_( Point1 );
    N_@{Index1}_ = GetN_x_( Point1 );
    F_@{Index1}_ = GetF_x_( Point1 );
    K_@{Index1}_ = GetK_x_( Point1 );
    H_@{Index1}_ = GetH_x_( Point1 );
    Q_@{Index1}_ = GetQ_x_( Point1 );
    SN_@{Index1}_ = GetSN_x_( Point1 , SpatialPointsPerDimension );
    SD_@{Index1}_ = GetSD_x_( Point1 , SpatialPointsPerDimension );
    A_@{Index1}_ = 1;
    A_@{Index1}_LEAD_ = GA_;
    N_LAG_@{Index1}_ = N_@{Index1}_ ./ GN_;
    P_Over_Q_@{Index1}_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GQTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;
    E_@{Index1}_  = E_by_F_@{Index1}_ * F_@{Index1}_;
    L_@{Index1}_  = thetaF * gamma / ( thetaL * E_by_F_@{Index1}_  + thetaF * gamma );
    ZF_@{Index1}_ = ( F_@{Index1}_ / L_@{Index1}_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
    SP_@{Index1}_ = ( 1 - gamma ) * F_@{Index1}_ / ZF_@{Index1}_;
    P_@{Index1}_ = P_Over_Q_@{Index1}_ * Q_@{Index1}_;
    Z_@{Index1}_ = ( ( ( K_@{Index1}_ / GYTrend_ ) ^ alpha * ( A_@{Index1}_ * H_@{Index1}_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_@{Index1}_ / P_@{Index1}_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
    SRK_@{Index1}_ = ( 1 - kappa ) * alpha * SP_@{Index1}_ * Z_@{Index1}_ / ( K_@{Index1}_ / GYTrend_ );
    W_@{Index1}_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_@{Index1}_ * Z_@{Index1}_ / H_@{Index1}_;
    C_@{Index1}_ = thetaC * F_@{Index1}_ / ( thetaF * P_@{Index1}_ );
    I_@{Index1}_ = K_@{Index1}_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
    M_@{Index1}_ = kappa * SP_@{Index1}_ * Z_@{Index1}_ / P_@{Index1}_;
    Y_@{Index1}_ = C_@{Index1}_ + I_@{Index1}_ + M_@{Index1}_;
    YBar_@{Index1}_ = 0
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2= IndicesStringArray[Point2]
        + (1 / (SpatialPointsPerDimension^2)) * Y_@{Index2}_ * P_@{Index2}_ ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance@{Index1}@{Index2}_ )
    @#endfor
    ;
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        SN@{Index1}@{Index2}_ = GetSNxx_( Point1 , Point2 );
    @#endfor
    J_@{Index1}_ = ( Z_@{Index1}_ - ZF_@{Index1}_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_@{Index1}_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_@{Index1}_ );
    U@{Index1}_ = ( C@{Index1}_ / N@{Index1}_LAG_ ) ^ thetaC
    * ( E@{Index1}_ / N_LAG_@{Index1}_ ) ^ thetaF
    * ( ( 1 - L@{Index1}_ ) / N_LAG_@{Index1}_ ) ^ thetaL
    * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_ / N_LAG_@{Index1}_ ) ^ ( 1 + nu ) ) ^ thetaH
    * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N_LAG_@{Index1}_ / N_LAG_ ) ^ 2 ) ^ thetaN
    * ( 1 - SN@{Index1}_ / N_LAG_@{Index1}_ ) ^ psi1
    * ( dBar - SD@{Index1}_ / SN@{Index1}_ ) ^ psi2
    * exp( psi3 * ( 0
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        + (1 / (SpatialPointsPerDimension^2)) * N_LAG_@{Index1}_ / N_LAG_ * log( SN@{Index1}@{Index2}_ / N_LAG_@{Index1}_ )
    @#endfor
    ) );
    U_@{Index1}_LEAD_ = U_@{Index1}_ * GUTrend_;
    H_@{Index1}_LEAD_ = H_@{Index1}_ * GN_;
    muN_@{Index1}_ = beta .* ( ( U_@{Index1}_ .* GUTrend_ ) ^ ( 1 - varsigma ) ) .* ( 1 + ( 1 - varsigma ) .* ( ...
            thetaH .* ( H_@{Index1}_ ./ N_LAG_@{Index1}_ ) .^ ( 1 + nu ) ./ ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_@{Index1}_ ./ N_LAG_@{Index1}_ ) .^ ( 1 + nu ) ) ...
            - thetaN .* log( N_@{Index1}_ ) ./ ( 1 ./ 2 .* Omega .^ 2 - 1 ./ 2 .* log( N_@{Index1}_ ) .^ 2 ) ...
            + psi1 .* GN_ .* SN_@{Index1}_ ./ ( N_@{Index1}_ - GN_ .* SN_@{Index1}_ ) ...
            - ( thetaC + thetaF + thetaL + psi3 ) ...
        ) ) ./ ( 1 - beta .* GmuNTrend_ .* GN_ );
@#endfor

@#for VariableName in AggregatedVariables
    @{VariableName}_ = 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        +  (1 / (SpatialPointsPerDimension^2)) * @{VariableName}@{Index1}_
    @#endfor
    ;
@#endfor

 @#include "InsertNewEndSteadyStateEquations.mod"
