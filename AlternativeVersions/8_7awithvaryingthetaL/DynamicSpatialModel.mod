@#includepath "DynareTransformationEngine"

@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1

@#define SpatialDimensions = 2
@#define SpatialPointsPerDimension = 2
@#define SpatialShape = "Torus"

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
@#define ShockProcesses = ShockProcesses + [ "MP", "0", "Inf", "1", "0", "0.01" ]

@#define PureTrendEndoVariables = PureTrendEndoVariables + [ "AP", "GA" ]

@#include "CreateShocks.mod"

@#define EndoVariables = EndoVariables + [ "R", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GYTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GZTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GFTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GPTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GUTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GSPTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GYBarTrend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GOmega1Trend", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "GmuNTrend", "0", "Inf", "1" ]

@#for Point in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point]
    @#define EndoVariables = EndoVariables + [ "Cr" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "Cs" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "K" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "I" + Index1, "0", "Inf", "GYTrend" ]
    @#define EndoVariables = EndoVariables + [ "Er" + Index1, "0", "Inf", "GFTrend" ]
    @#define EndoVariables = EndoVariables + [ "Es" + Index1, "0", "Inf", "GFTrend" ]
    @#define EndoVariables = EndoVariables + [ "F" + Index1, "0", "Inf", "GFTrend" ]
    @#define EndoVariables = EndoVariables + [ "Q" + Index1, "0", "Inf", "GPTrend" ]
    @#define EndoVariables = EndoVariables + [ "J" + Index1, "0", "Inf", "GZTrend" ]
    @#define EndoVariables = EndoVariables + [ "L" + Index1, "0", "Inf", "1" ]
    @#define EndoVariables = EndoVariables + [ "Dr" + Index1, "0", "Inf", "1" ]
    @#define EndoVariables = EndoVariables + [ "Ds" + Index1, "0", "Inf", "1" ]
    @#define EndoVariables = EndoVariables + [ "Hr" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "Hs" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "Nr" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "Ns" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "SNr" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "SNs" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "SDr" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "SDs" + Index1, "0", "Inf", "GN" ]
    @#define EndoVariables = EndoVariables + [ "muNr" + Index1, "0", "Inf", "GmuNTrend" ]
    @#define EndoVariables = EndoVariables + [ "muNs" + Index1, "0", "Inf", "GmuNTrend" ]
    @#define EndoVariables = EndoVariables + [ "Ur" + Index1, "0", "Inf", "GUTrend" ]
    @#define EndoVariables = EndoVariables + [ "Us" + Index1, "0", "Inf", "GUTrend" ]
    @#define EndoVariables = EndoVariables + [ "SPJ1" + Index1, "0", "Inf", "GSPTrend" ]
    @#define EndoVariables = EndoVariables + [ "Omega1" + Index1, "0", "Inf", "GOmega1Trend" ]
    @#define EndoVariables = EndoVariables + [ "Omega2" + Index1, "0", "Inf", "GYBarTrend" ]
    @#define EndoVariables = EndoVariables + [ "varphi1" + Index1, "0", "Inf", "GYBarTrend" ]
    @#define EndoVariables = EndoVariables + [ "varphi2" + Index1, "0", "Inf", "GOmega1Trend" ]
    @#define EndoVariables = EndoVariables + [ "varphi3" + Index1, "0", "Inf", "GSPTrend" ]
    @#define EndoVariables = EndoVariables + [ "SPJ2" + Index1, "0", "Inf", "GSPTrend" ]
@#endfor

@#define EndoVariables = EndoVariables + [ "Y", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "C", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "K", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "I", "0", "Inf", "GYTrend" ]
@#define EndoVariables = EndoVariables + [ "F", "0", "Inf", "GFTrend" ]
@#define EndoVariables = EndoVariables + [ "J", "0", "Inf", "GZTrend" ]
@#define EndoVariables = EndoVariables + [ "L", "0", "Inf", "1" ]
@#define EndoVariables = EndoVariables + [ "H", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "SN", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "SDr", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "SDs", "0", "Inf", "GN" ]
@#define EndoVariables = EndoVariables + [ "Ur", "0", "Inf", "GUTrend" ]
@#define EndoVariables = EndoVariables + [ "Us", "0", "Inf", "GUTrend" ]
@#define EndoVariables = EndoVariables + [ "Pi", "0", "Inf", "1" ]

@#define AggregatedVariables = [ "Y", "C", "K", "I", "F", "J", "L", "H", "SN", "SDr", "SDs", "Us", "Ur" ]

@#include "ClassifyDeclare.mod"

parameters alpha gamma kappa nu varsigma zeta lambda deltaJ deltaK Phi2 PhiL thetaC thetaF thetaL thetaH thetaN psi1 psi2 psi3 Gamma Omega dBar;
parameters UtilityParamSum Pi_star xi rho_R theta_y theta_pi theta beta_r_by_beta;

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
PhiL = param_PhiL;
theta = 0.5; % share of impatient households
beta_r_by_beta = 0.95; % impatient beta / patient beta

Pi_star = 1.005;
xi = param_xi;
rho_R = param_rho_R;
theta_y = param_theta_y;
theta_pi = param_theta_pi;

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

model;

    @#include "InsertNewModelEquations.mod"

   GYTrend    = ( GA * GN ) ^ ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) / ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) - lambda ) );
   GZTrend    = GYTrend ^ ( 1 / ( 1 + lambda ) );
// GJTrend    = GYTrend ^ ( 1 / ( 1 + lambda ) );
   GFTrend    = GYTrend ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
// GSRLTrend  = GYTrend ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
// GWTrend    = GYTrend ^ ( ( 1 - gamma ) / ( 1 + lambda ) ) / GN;
   GPTrend    = GYTrend ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
// GQTrend    = GYTrend ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
// GSRKTrend  = GYTrend ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
   GSPTrend   = GYTrend ^ ( - gamma / ( 1 + lambda ) );
   GYBarTrend = GYTrend ^ ( - gamma / lambda );
   GUTrend    = GYTrend ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN ^ ( thetaC + thetaF + thetaL );
   GmuNTrend  = GYTrend ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );
   GOmega1Trend   = GYTrend ^ ( - gamma / lambda - gamma / ( 1 + lambda ) );

   #betar = beta * beta_r_by_beta;

   @#for Point1 in 1 : SpatialNumPoints
       @#define Index1 = IndicesStringArray[Point1]
       #N@{Index1} = Ns@{Index1} + Nr@{Index1};
       #N@{Index1}_LAG = Ns@{Index1}_LAG + Nr@{Index1}_LAG;
   @#endfor

    #N = 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + Weight@{Index1} * N@{Index1}
    @#endfor
    ;

    #N_LAG = 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + Weight@{Index1} * N@{Index1}_LAG
    @#endfor
    ;

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        #A@{Index1} = AP * AT@{Index1};
        #A@{Index1}_LEAD = AP_LEAD * AT@{Index1}_LEAD;

        #ZF@{Index1} = ( F@{Index1} / L@{Index1}_LAG ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        #SRL@{Index1}_LEAD = gamma * F@{Index1}_LEAD / L@{Index1};
        #SP@{Index1} = ( 1 - gamma ) * F@{Index1} / ZF@{Index1};

        #ZF@{Index1}_LEAD = ( F@{Index1}_LEAD / L@{Index1} ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        #SP@{Index1}_LEAD = ( 1 - gamma ) * F@{Index1}_LEAD / ZF@{Index1}_LEAD;

        #P_star@{Index1} = ( 1 + lambda ) * Omega1@{Index1} / Omega2@{Index1};
        #P_star@{Index1}_LEAD = ( 1 + lambda ) * Omega1@{Index1}_LEAD / Omega2@{Index1}_LEAD;

        K@{Index1} = ( 1 - deltaK ) * K@{Index1}_LAG + ( 1 - Phi2 / 2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) ^ 2 ) * I@{Index1};

        #SN@{Index1} = SNr@{Index1} + SNs@{Index1};
        #SN@{Index1}_LEAD = SNr@{Index1}_LEAD + SNs@{Index1}_LEAD;
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            #SNr@{Index1}@{Index2} = psi3 * N@{Index2}_LAG / N_LAG / ( ( muNr@{Index1} - muNr@{Index2} ) / ( ( 1 - varsigma ) * Nr@{Index1}_LAG * Ur@{Index1} ^ ( 1 - varsigma ) ) + psi1 / ( Nr@{Index1}_LAG - SNr@{Index1} ) + psi2 * ( Distance@{Index1}@{Index2} * SNr@{Index1} - SDr@{Index1} ) / ( dBar * SNr@{Index1} * SNr@{Index1} - SNr@{Index1} * SDr@{Index1} ) );
            #SNs@{Index1}@{Index2} = psi3 * N@{Index2}_LAG / N_LAG / ( ( muNs@{Index1} - muNs@{Index2} ) / ( ( 1 - varsigma ) * Ns@{Index1}_LAG * Us@{Index1} ^ ( 1 - varsigma ) ) + psi1 / ( Ns@{Index1}_LAG - SNs@{Index1} ) + psi2 * ( Distance@{Index1}@{Index2} * SNs@{Index1} - SDs@{Index1} ) / ( dBar * SNs@{Index1} * SNs@{Index1} - SNs@{Index1} * SDs@{Index1} ) );
            #SN@{Index1}@{Index2} = SNr@{Index1}@{Index2} + SNs@{Index1}@{Index2};
        @#endfor

        SNr@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * SNr@{Index1}@{Index2}
        @#endfor
        ;
        SNs@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * SNs@{Index1}@{Index2}
        @#endfor
        ;

        SDr@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * Distance@{Index1}@{Index2} * SNr@{Index1}@{Index2}
        @#endfor
        ;

        SDs@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * Distance@{Index1}@{Index2} * SNs@{Index1}@{Index2}
        @#endfor
        ;

        [name = 'Ur@{Index1}'] Ur@{Index1} = ( Cr@{Index1} / Nr@{Index1}_LAG ) ^ thetaC
        * ( Er@{Index1} / Nr@{Index1}_LAG ) ^ thetaF
        * ( ( Dr@{Index1} ) / Nr@{Index1}_LAG ) ^ thetaL
        * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hr@{Index1} / Nr@{Index1}_LAG ) ^ ( 1 + nu ) ) ^ thetaH
        * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1}_LAG / N_LAG ) ^ 2 ) ^ thetaN
        * ( 1 - SNr@{Index1} / Nr@{Index1}_LAG ) ^ psi1
        * ( dBar - SDr@{Index1} / SNr@{Index1} ) ^ psi2
        * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * N@{Index2}_LAG / N_LAG * log( SNr@{Index1}@{Index2} / Nr@{Index1}_LAG )
        @#endfor
        ) );

        [name = 'Us@{Index1}'] Us@{Index1} = ( Cs@{Index1} / Ns@{Index1}_LAG ) ^ thetaC
        * ( Es@{Index1} / Ns@{Index1}_LAG ) ^ thetaF
        * ( ( Ds@{Index1} ) / Ns@{Index1}_LAG ) ^ thetaL
        * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hs@{Index1} / Ns@{Index1}_LAG ) ^ ( 1 + nu ) ) ^ thetaH
        * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1}_LAG / N_LAG ) ^ 2 ) ^ thetaN
        * ( 1 - SNs@{Index1} / Ns@{Index1}_LAG ) ^ psi1
        * ( dBar - SDs@{Index1} / SNs@{Index1} ) ^ psi2
        * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * N@{Index2}_LAG / N_LAG * log( SNs@{Index1}@{Index2} / Ns@{Index1}_LAG )
        @#endfor
        ) );

        [name = 'muNr@{Index1}'] muNr@{Index1} = betar * ( muNr@{Index1}_LEAD * GN_LEAD + Ur@{Index1}_LEAD ^ ( 1 - varsigma ) + ( 1 - varsigma ) * Ur@{Index1}_LEAD ^ ( 1 - varsigma ) * (
            thetaH * ( Hr@{Index1}_LEAD / Nr@{Index1} ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hr@{Index1}_LEAD / Nr@{Index1} ) ^ ( 1 + nu ) )
            - thetaN * ( Nr@{Index1} / N@{Index1}  ) * log( N@{Index1} / N ) / ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1} / N ) ^ 2 )
            + psi1 * SNr@{Index1}_LEAD / ( Nr@{Index1} - SNr@{Index1}_LEAD )
            - ( thetaC + thetaF + thetaL + psi3 )
        ) );

        [name = 'muNs@{Index1}'] muNs@{Index1} = beta * ( muNs@{Index1}_LEAD * GN_LEAD + Us@{Index1}_LEAD ^ ( 1 - varsigma ) + ( 1 - varsigma ) * Us@{Index1}_LEAD ^ ( 1 - varsigma ) * (
            thetaH * ( Hs@{Index1}_LEAD / Ns@{Index1} ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hs@{Index1}_LEAD / Ns@{Index1} ) ^ ( 1 + nu ) )
            - thetaN * ( Ns@{Index1} / N@{Index1}  ) * log( N@{Index1} / N ) / ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1} / N ) ^ 2 )
            + psi1 * SNs@{Index1}_LEAD / ( Ns@{Index1} - SNs@{Index1}_LEAD )
            - ( thetaC + thetaF + thetaL + psi3 )
        ) );
    @#endfor

    @#define Index1 = IndicesStringArray[1]
    #Xi_LEAD = beta * ( Ns@{Index1} / Ns@{Index1}_LAG ) * ( Es@{Index1} / Es@{Index1}_LEAD ) * ( Us@{Index1}_LEAD / Us@{Index1} ) ^ ( 1 - varsigma );

    @#for Point2 in 2 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        Er@{Index1} / Nr@{Index1}_LAG / Ur@{Index1} ^ ( 1 - varsigma ) = Er@{Index2} / Nr@{Index2}_LAG / Ur@{Index2} ^ ( 1 - varsigma );
        Es@{Index1} / Ns@{Index1}_LAG / Us@{Index1} ^ ( 1 - varsigma ) = Es@{Index2} / Ns@{Index2}_LAG / Us@{Index2} ^ ( 1 - varsigma );
    @#endfor


    @#for Point1 in 2 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        Ns@{Index1} = GN * Ns@{Index1}_LAG - SNs@{Index1}
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * SNs@{Index2}@{Index1}
        @#endfor
        ;

        Nr@{Index1} = GN * Nr@{Index1}_LAG - SNr@{Index1}
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * SNr@{Index2}@{Index1}
        @#endfor
        ;
    @#endfor

    [name = 'Nr'] theta = 0
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        + Weight@{Index2} * Nr@{Index2}
    @#endfor
    ;
    [name = 'Ns'] 1 - theta = 0
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        + Weight@{Index2} * Ns@{Index2}
    @#endfor
    ;

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        [name = 'SPJ1'] SPJ1@{Index1} = ( xi * (1-deltaJ) * ( J@{Index1}_LAG / J@{Index1} ) * ( SPJ1@{Index1}_LAG / Pi )^( -1 / lambda ) + ( 1 - xi * (1-deltaJ) * J@{Index1}_LAG / J@{Index1} ) * P_star@{Index1}^( -1 / lambda ) )^( -lambda );
        [name = 'SPJ2'] SPJ2@{Index1} = ( xi * (1-deltaJ) * ( J@{Index1}_LAG / J@{Index1} ) * ( SPJ2@{Index1}_LAG / Pi )^( -(1+lambda) / lambda ) + ( 1 - xi * (1-deltaJ) * J@{Index1}_LAG / J@{Index1} ) * P_star@{Index1}^( -(1+lambda) / lambda ) )^( -lambda/(1+lambda) );

        #P@{Index1} =  ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} *
                ( J@{Index2} * ( SPJ1@{Index2} * exp( tau * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda ) )
        @#endfor
        ) ^ ( - lambda );

        #P@{Index1}_LAG =  ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} *
                ( J@{Index2}_LAG * ( SPJ1@{Index2}_LAG * exp( tau_LAG * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda ) )
        @#endfor
        ) ^ ( - lambda );

        #P@{Index1}_LEAD =  ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} *
                ( J@{Index2}_LEAD * ( SPJ1@{Index2}_LEAD * exp( tau_LEAD * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda ) )
        @#endfor
        ) ^ ( - lambda );

        #H@{Index1} = Hs@{Index1} + Hr@{Index1};
        #H@{Index1}_LEAD = Hs@{Index1}_LEAD + Hr@{Index1}_LEAD;
        #Z@{Index1} = ( ( K@{Index1}_LAG ^ alpha * ( A@{Index1} * H@{Index1} ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1} / P@{Index1} ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        #Z@{Index1}_LEAD = ( ( K@{Index1} ^ alpha * ( A@{Index1}_LEAD * H@{Index1}_LEAD ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1}_LEAD / P@{Index1}_LEAD ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        #M@{Index1} = kappa * SP@{Index1} * Z@{Index1} / P@{Index1};
        #M@{Index1}_LEAD = kappa * SP@{Index1}_LEAD * Z@{Index1}_LEAD / P@{Index1}_LEAD;
        #SRK@{Index1} = ( 1 - kappa ) * alpha * SP@{Index1} * Z@{Index1} / K@{Index1}_LAG;
        #SRK@{Index1}_LEAD = ( 1 - kappa ) * alpha * SP@{Index1}_LEAD * Z@{Index1}_LEAD / K@{Index1};
        #W@{Index1} = ( 1 - kappa ) * ( 1 - alpha ) * SP@{Index1} * Z@{Index1} / H@{Index1};
        #SRD@{Index1} = thetaL * Es@{Index1} / ( thetaF *  Ds@{Index1} );
        #SRD@{Index1}_LEAD = thetaL * Es@{Index1}_LEAD / ( thetaF *  Ds@{Index1}_LEAD );
        #D@{Index1} = 1 - L@{Index1}_LAG;
        #D@{Index1}_LEAD = 1 - L@{Index1};
        #E@{Index1} = Es@{Index1} + Er@{Index1};
        #C@{Index1} = Cs@{Index1} + Cr@{Index1};
        #C@{Index1}_LEAD = Cs@{Index1}_LEAD + Cr@{Index1}_LEAD;

        D@{Index1} = Ds@{Index1} + Dr@{Index1};
        SRD@{Index1} = thetaL * Er@{Index1} / ( thetaF *  Dr@{Index1} );

        #PR@{Index1} = Xi_LEAD * SRD@{Index1}_LEAD;
        #PF@{Index1} = Xi_LEAD * SRL@{Index1}_LEAD;

        [name = 'Euler K'] 1 = Xi_LEAD * ( SRK@{Index1}_LEAD + Q@{Index1}_LEAD * ( 1 - deltaK ) ) / Q@{Index1};
        [name = 'Q'] P@{Index1} = Q@{Index1} * ( 1 - Phi2 / 2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) ^ 2 - Phi2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) * I@{Index1} / I@{Index1}_LAG ) + Xi_LEAD * Q@{Index1}_LEAD * Phi2 * ( I@{Index1}_LEAD / I@{Index1} - 1 ) * ( I@{Index1}_LEAD / I@{Index1} ) ^ 2;
        #Y@{Index1} = C@{Index1} + I@{Index1} + M@{Index1};
        #Y@{Index1}_LEAD = C@{Index1}_LEAD + I@{Index1}_LEAD + M@{Index1}_LEAD;
        thetaC * Er@{Index1} = thetaF * P@{Index1} * Cr@{Index1};
        thetaC * Es@{Index1} = thetaF * P@{Index1} * Cs@{Index1};
        Xi_LEAD * SRD@{Index1}_LEAD = Xi_LEAD * SRL@{Index1}_LEAD - PhiL * ( L@{Index1} / L@{Index1}_LAG - 1 ) / L@{Index1}_LAG + Xi_LEAD * PhiL * ( L@{Index1}_LEAD / L@{Index1} - 1 ) * ( L@{Index1}_LEAD / ( L@{Index1}^2 ) );
        thetaH * ( Hr@{Index1} / Nr@{Index1}_LAG ) ^ nu = thetaF * Nr@{Index1} / Er@{Index1} * W@{Index1} * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hr@{Index1} / Nr@{Index1}_LAG ) ^ ( 1 + nu ) );
        thetaH * ( Hs@{Index1} / Ns@{Index1}_LAG ) ^ nu = thetaF * Ns@{Index1} / Es@{Index1} * W@{Index1} * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hs@{Index1} / Ns@{Index1}_LAG ) ^ ( 1 + nu ) );
    @#endfor


    0 = 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + Weight@{Index1} * ( P@{Index1} * Cr@{Index1} + Er@{Index1} + SRD@{Index1} * Dr@{Index1} - W@{Index1} * Hr@{Index1})
    @#endfor
    ;

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        #YBar@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * Y@{Index2} * P@{Index2} ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance@{Index1}@{Index2} )
        @#endfor
        ;
        #YBar@{Index1}_LEAD = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * Y@{Index2}_LEAD * P@{Index2}_LEAD ^ ( ( 1 + lambda ) / lambda ) * exp( - tau_LEAD / lambda * Distance@{Index1}@{Index2} )
        @#endfor
        ;

        [name = 'Omega1'] Omega1@{Index1} = YBar@{Index1} * SP@{Index1} + xi * Xi_LEAD * ( Pi_LEAD^((1+lambda)/lambda) ) * Omega1@{Index1}_LEAD;
        [name = 'Omega2'] Omega2@{Index1} = YBar@{Index1} + xi * Xi_LEAD * ( Pi_LEAD^(1/lambda) ) * Omega2@{Index1}_LEAD;

        #pis@{Index1} = YBar@{Index1} * ( P_star@{Index1}^( - 1/lambda ) - SP@{Index1} * P_star@{Index1}^( -(1 + lambda)/lambda ) ) + varphi1@{Index1} * P_star@{Index1}^( - 1/lambda ) - varphi2@{Index1} * P_star@{Index1}^( - (1+lambda)/lambda ) + varphi3@{Index1} * ( 1 - xi ) / xi;
        #pis@{Index1}_LEAD = YBar@{Index1}_LEAD * ( P_star@{Index1}_LEAD^( - 1/lambda ) - SP@{Index1}_LEAD * P_star@{Index1}_LEAD^( -(1 + lambda)/lambda ) ) + varphi1@{Index1}_LEAD * P_star@{Index1}_LEAD^( - 1/lambda ) - varphi2@{Index1}_LEAD * P_star@{Index1}_LEAD^( - (1+lambda)/lambda ) + varphi3@{Index1}_LEAD * ( 1 - xi ) / xi;

        [name = 'varphi1'] varphi1@{Index1} = (1-deltaJ) * xi * Xi_LEAD * ( Pi_LEAD^(1/lambda) ) * ( YBar@{Index1}_LEAD +  varphi1@{Index1}_LEAD );
        [name = 'varphi2'] varphi2@{Index1} = (1-deltaJ) * xi * Xi_LEAD * ( Pi_LEAD^((1+lambda)/lambda) ) * ( YBar@{Index1}_LEAD * SP@{Index1}_LEAD +  varphi2@{Index1}_LEAD );
        [name = 'varphi3'] varphi3@{Index1} = (1-deltaJ) * xi * Xi_LEAD * ( pis@{Index1}_LEAD +  varphi3@{Index1}_LEAD );
        phi * SP@{Index1} = pis@{Index1};

        [name = 'Z demand'] Z@{Index1} = ZF@{Index1} + phi * ( J@{Index1} - ( 1 - deltaJ ) * J@{Index1}_LAG ) + J@{Index1} * SPJ2@{Index1}^(-(1+lambda)/lambda) *  YBar@{Index1};
    @#endfor

    [name = 'Euler equation'] 1 = R * Xi_LEAD / Pi_LEAD;

    % Consumption good inflation:
%    #PiC = Pi * ( 0
%    @#for Point1 in 1 : SpatialNumPoints
%        @#define Index1 = IndicesStringArray[Point1]
%        + Weight@{Index1} * ( P@{Index1} / P@{Index1}_LAG )
%    @#endfor
%    );

    % Food price inflation:
    #PiC = Pi;

    [name = 'Taylor rule'] R = R_LAG^rho_R * ( STEADY_STATE(R) * ( Y / STEADY_STATE(Y) )^theta_y * ( PiC / Pi_star )^theta_pi  )^( 1 - rho_R ) * MP ;%* AT_2_2_LAG;

    [name = 'Food market clearing']  0 = 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + Weight@{Index1} * ( E@{Index1} - F@{Index1} )
    @#endfor
    ;

    @#for VariableName in AggregatedVariables
        [name = 'Aggregete @{VariableName}'] @{VariableName} = 0
        @#for Point1 in 1 : SpatialNumPoints
            @#define Index1 = IndicesStringArray[Point1]
            + Weight@{Index1} * @{VariableName}@{Index1}
        @#endfor
        ;
    @#endfor

end;

@#if SpatialShape[1] == "P"
    @#define LoadSteadyState = 1
@#else
    @#define LoadSteadyState = 0
@#endif

@#if LoadSteadyState

    load_params_and_steady_state( 'SteadyState.txt' );

@#else

    steady_state_model;
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
           GYBarTrend_ = GYTrend_ ^ ( - gamma / lambda );
           GUTrend_    = GYTrend_ ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN_ ^ ( thetaC + thetaF + thetaL );
           GmuNTrend_  = GYTrend_ ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN_ ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );
           GOmega1Trend_ = GYTrend_ ^ ( - gamma / lambda - gamma / ( 1 + lambda ) );

        PiC_ = Pi_star;
%        Pi_ = PiC_ / GPTrend_;
        Pi_ = PiC_;
        Pi_LEAD_ = Pi_;
        Xi_LEAD_ = beta_ * GN_ / GFTrend_ * GUTrend_ ^ ( 1 - varsigma );

        R_ = Pi_LEAD_ / Xi_LEAD_;

        N_1_ = 1;
        N_ = N_1_;
        Nr_1_ = theta;
        Ns_1_ = 1-theta;

        N_1_LAG_ = 1 / GN_;
        N_LAG_ = N_1_LAG_;
        Nr_1_LAG_ = theta / GN_;
        Ns_1_LAG_ = (1-theta) / GN_;

        L_1_ = thetaF * gamma / ( thetaL + thetaF * gamma );
        betar_ = beta_ * beta_r_by_beta;

        SNr_1_ = psi3 / ( psi1 + psi3 ) * Nr_1_LAG_;
        SNs_1_ = psi3 / ( psi1 + psi3 ) * Ns_1_LAG_;
        SN_1_ = psi3 / ( psi1 + psi3 ) * N_1_LAG_;

        @#define Index1 = IndicesStringArray[1]
        SD_1_ = GetSD_1_( psi1, psi2, psi3, N_1_LAG_, SN_1_, dBar
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            , Distance@{Index1}@{Index2}_
        @#endfor
        );
        SDr_1_ = theta * SD_1_;
        SDs_1_ = (1 - theta) * SD_1_;

        A_1_ = 1;
        A_1_LEAD_ = GA_;

        AverageTransportCost_ = ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2}_ * exp( - tau_ / lambda * Distance@{Index1}@{Index2}_ )
        @#endfor
        ) ^ ( - lambda );

        P_1_Over_Q_1_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GQTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;

        F_1_ = GetF_1_( A_1_, Nr_1_, Nr_1_LAG_, Ns_1_, Ns_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GZTrend_, xi, GPTrend_, GYBarTrend_, GOmega1Trend_, Pi_);
        K_1_ = GetK_1_( 1 );
        H_1_ = GetH_1_( 1 );
        Q_1_ = GetQ_1_( 1 );
        Er_by_Es_1_ = GetEr_by_Es_1_( 1 );
        Hr_by_Hs_1_ = GetHr_by_Hs_1_( 1 );

        ZF_1_ = ( F_1_ / L_1_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        SP_1_ = ( 1 - gamma ) * F_1_ / ZF_1_;
        P_1_ = P_1_Over_Q_1_ * Q_1_;
        Z_1_ = ( ( ( K_1_ / GYTrend_ ) ^ alpha * ( A_1_ * H_1_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_1_ / P_1_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        SRK_1_ = ( 1 - kappa ) * alpha * SP_1_ * Z_1_ / ( K_1_ / GYTrend_ );
        W_1_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_1_ * Z_1_ / H_1_;
        Es_1_ = F_1_ / ( 1 + Er_by_Es_1_ );
        Er_1_ = Es_1_ * Er_by_Es_1_;
        Hs_1_ = H_1_ / ( 1 + Hr_by_Hs_1_ );
        Hr_1_ = Hs_1_ * Hr_by_Hs_1_;
        Cr_1_ = thetaC * Er_1_ / ( thetaF * P_1_ );
        Cs_1_ = thetaC * Es_1_ / ( thetaF * P_1_ );
        C_1_ = Cr_1_ + Cs_1_;
        SRD_1_ = gamma * F_1_ / L_1_;
        Dr_1_ = thetaL / thetaF * Er_1_ / SRD_1_;
        Ds_1_ = thetaL / thetaF * Es_1_ / SRD_1_;
        I_1_ = K_1_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
        M_1_ = kappa * SP_1_ * Z_1_ / P_1_;
        Y_1_ = C_1_ + I_1_ + M_1_;
        YBar_1_ = Y_1_ * P_1_ ^ ( ( 1 + lambda ) / lambda ) * AverageTransportCost_ ^ ( - 1 / lambda );
        P_star_1_ = ( 1 + lambda ) * SP_1_ * ( 1 - xi * Xi_LEAD_ * ( Pi_LEAD_^(1/lambda) ) * GYBarTrend_ ) / ( 1 - xi * Xi_LEAD_ * ( Pi_LEAD_^((1+lambda)/lambda) ) * GOmega1Trend_ );
        SPJ1_1_ = P_star_1_ * ( ( 1 - xi * (1-deltaJ) / GZTrend_ ) / ( 1 - xi * (1-deltaJ) * ( GSPTrend_ * Pi_ )^(1/lambda) / GZTrend_  ) )^(-lambda);
        J_1_ = ( SPJ1_1_ * AverageTransportCost_ / P_1_ ) ^ ( 1 / lambda );
        E_1_ = F_1_;

        Omega1_1_ = YBar_1_ * SP_1_ / ( 1 - xi * Xi_LEAD_ * ( Pi_LEAD_^((1+lambda)/lambda) ) * GOmega1Trend_ );
        Omega2_1_ = YBar_1_ / ( 1 - xi * Xi_LEAD_ * ( Pi_LEAD_^(1/lambda) ) * GYBarTrend_ );

        varphi1_1_ = (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^(1/lambda) ) * GYBarTrend_ * YBar_1_ / ( 1 - (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^(1/lambda) ) * GYBarTrend_  );
        varphi2_1_ = (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^((1+lambda)/lambda) ) * GOmega1Trend_ * YBar_1_ * SP_1_ / ( 1 - (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^((1+lambda)/lambda) ) * GOmega1Trend_ );
        pis_1_ = ( YBar_1_ * ( P_star_1_^( - 1/lambda ) - SP_1_ * P_star_1_^( -(1 + lambda)/lambda ) ) + varphi1_1_ * P_star_1_^( - 1/lambda ) - varphi2_1_ * P_star_1_^( - (1+lambda)/lambda ) ) / ( 1 - ( ( (1-deltaJ) * xi * Xi_LEAD_ * GSPTrend_ / ( 1 -  (1-deltaJ) * xi * Xi_LEAD_ * GSPTrend_ ) ) * ( 1 - xi ) / xi ) );
        varphi3_1_ = (1-deltaJ) * xi * Xi_LEAD_ * GSPTrend_ * pis_1_ / ( 1 -  (1-deltaJ) * xi * Xi_LEAD_ * GSPTrend_ );

        SPJ2_1_ = P_star_1_ * ( ( 1 - xi * (1-deltaJ) / GZTrend_ ) / ( 1 - xi * (1-deltaJ) * ( GSPTrend_ * Pi_ )^((1+lambda)/lambda) / GZTrend_  ) )^(-lambda/(1+lambda));

        @#for Point1 in 1 : SpatialNumPoints
            @#define Index1 = IndicesStringArray[Point1]
            @#for Point2 in 1 : SpatialNumPoints
                @#define Index2 = IndicesStringArray[Point2]
                SN@{Index1}@{Index2}_ = psi3 / ( psi1 / ( N_1_LAG_ - SN_1_ ) + psi2 * ( Distance@{Index1}@{Index2}_ * SN_1_ - SD_1_ ) / ( dBar * SN_1_ * SN_1_ - SN_1_ * SD_1_ ) );
                SNr@{Index1}@{Index2 }_ = theta * SN@{Index1}@{Index2}_;
                SNs@{Index1}@{Index2}_ = (1-theta) * SN@{Index1}@{Index2}_;
            @#endfor
        @#endfor

        Ur_1_ = ( Cr_1_ / Nr_1_LAG_ ) ^ thetaC
        * ( Er_1_ / Nr_1_LAG_ ) ^ thetaF
        * ( ( Dr_1_ ) / Nr_1_LAG_ ) ^ thetaL
        * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hr_1_ / Nr_1_LAG_ ) ^ ( 1 + nu ) ) ^ thetaH
         * ( 1 / 2 * Omega ^ 2 ) ^ thetaN
         * ( 1 - SN_1_ / N_1_LAG_ ) ^ psi1
         * ( dBar - SD_1_ / SN_1_ ) ^ psi2
         * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2}_ * log( SN_1_1@{Index2}_ / N_1_LAG_ )
        @#endfor
        ) );
        Us_1_ = ( Cs_1_ / Ns_1_LAG_ ) ^ thetaC
        * ( Es_1_ / Ns_1_LAG_ ) ^ thetaF
        * ( ( Ds_1_ ) / Ns_1_LAG_ ) ^ thetaL
        * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hs_1_ / Ns_1_LAG_ ) ^ ( 1 + nu ) ) ^ thetaH
         * ( 1 / 2 * Omega ^ 2 ) ^ thetaN
         * ( 1 - SN_1_ / N_1_LAG_ ) ^ psi1
         * ( dBar - SD_1_ / SN_1_ ) ^ psi2
         * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2}_ * log( SN_1_1@{Index2}_ / N_1_LAG_ )
        @#endfor
        ) );


        Ur_1_LEAD_ = Ur_1_ * GUTrend_;
        Us_1_LEAD_ = Us_1_ * GUTrend_;
        Hr_1_LEAD_ = Hr_1_ * GN_;
        Hs_1_LEAD_ = Hs_1_ * GN_;
        muNr_1_ = betar_ * ( Ur_1_LEAD_ ^ ( 1 - varsigma ) + ( 1 - varsigma ) * Ur_1_LEAD_ ^ ( 1 - varsigma ) * ( thetaH * ( Hr_1_LEAD_ / Nr_1_ ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hr_1_LEAD_ / Nr_1_ ) ^ ( 1 + nu ) ) + psi1 * SNr_1_ * GN_ / ( Nr_1_ - SNr_1_ * GN_ ) - ( thetaC + thetaF + thetaL + psi3 ) ) ) / ( 1 - betar_ * GmuNTrend_ * GN_ );
        muNs_1_ = beta_ * ( Us_1_LEAD_ ^ ( 1 - varsigma ) + ( 1 - varsigma ) * Us_1_LEAD_ ^ ( 1 - varsigma ) * ( thetaH * ( Hs_1_LEAD_ / Ns_1_ ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hs_1_LEAD_ / Ns_1_ ) ^ ( 1 + nu ) ) +  psi1 * SNs_1_ * GN_ / ( Ns_1_ - SNs_1_ * GN_ ) - ( thetaC + thetaF + thetaL + psi3 ) ) ) / ( 1 - beta_ * GmuNTrend_ * GN_ );

        @#for VariableName in AggregatedVariables
            @{VariableName}_ = @{VariableName}_1_;
        @#endfor

        @#for Point1 in 1 : SpatialNumPoints
            @#define Index1 = IndicesStringArray[Point1]
            Cr@{Index1}_ = Cr_1_;
            Cs@{Index1}_ = Cs_1_;
            Dr@{Index1}_ = Dr_1_;
            Ds@{Index1}_ = Ds_1_;
            K@{Index1}_ = K_1_;
            I@{Index1}_ = I_1_;
            Er@{Index1}_ = Er_1_;
            Es@{Index1}_ = Es_1_;
            F@{Index1}_ = F_1_;
            Q@{Index1}_ = Q_1_;
            J@{Index1}_ = J_1_;
            L@{Index1}_ = L_1_;
            Hr@{Index1}_ = Hr_1_;
            Hs@{Index1}_ = Hs_1_;
            Nr@{Index1}_ = Nr_1_;
            Ns@{Index1}_ = Ns_1_;
            SNr@{Index1}_ = SNr_1_;
            SNs@{Index1}_ = SNs_1_;
            SDr@{Index1}_ = SDr_1_;
            SDs@{Index1}_ = SDs_1_;
            muNr@{Index1}_ = muNr_1_;
            muNs@{Index1}_ = muNs_1_;
            Ur@{Index1}_ = Ur_1_;
            Us@{Index1}_ = Us_1_;
            SPJ1@{Index1}_ = SPJ1_1_;
            Omega1@{Index1}_ = Omega1_1_;
            Omega2@{Index1}_ = Omega2_1_;
            varphi1@{Index1}_ = varphi1_1_;
            varphi2@{Index1}_ = varphi2_1_;
            varphi3@{Index1}_ = varphi3_1_;
            SPJ2@{Index1}_ = SPJ2_1_;
        @#endfor

        @#include "InsertNewEndSteadyStateEquations.mod"
    end;

@#endif

@#define Deterministic = 0

shocks;
    @#if Deterministic
        // This is row 37 of sqrtm( M_.Sigma_e( 6:end, 6:end ) ) when the model is run with Deterministic = 0
        @#define ImpulseValues = [ "0.000696148115243237", "0.000371150838989358", "-8.78536948609384e-05", "-0.000422908976433854", "-0.000584004323972896", "-0.000422908976434025", "-8.78536948608953e-05", "0.000371150838989199", "0.000371150838988895", "0.00257404871124431", "0.00515834258373992", "0.00869684235182027", "0.0105438251489022", "0.00869684235181994", "0.00515834258373983", "0.00257404871124433", "-8.78536948610069e-05", "0.00515834258373976", "0.0129027156681362", "0.0274642826537827", "0.0372639907921042", "0.0274642826537831", "0.0129027156681361", "0.00515834258373973", "-0.000422908976434043", "0.00869684235182006", "0.0274642826537831", "0.0854621538241213", "0.154136302618233", "0.0854621538241214", "0.0274642826537830", "0.00869684235182020", "-0.000584004323973137", "0.0105438251489024", "0.0372639907921038", "0.154136302618233", "0.928550829279997", "0.154136302618233", "0.0372639907921042", "0.0105438251489024", "-0.000422908976434094", "0.00869684235182022", "0.0274642826537827", "0.0854621538241215", "0.154136302618233", "0.0854621538241213", "0.0274642826537830", "0.00869684235182034", "-8.78536948610440e-05", "0.00515834258374003", "0.0129027156681362", "0.0274642826537828", "0.0372639907921042", "0.0274642826537829", "0.0129027156681359", "0.00515834258373981", "0.000371150838989175", "0.00257404871124437", "0.00515834258373976", "0.00869684235182014", "0.0105438251489024", "0.00869684235182007", "0.00515834258373997", "0.00257404871124420" ]
        @#for Point1 in 1 : SpatialNumPoints
            @#define Index1 = IndicesStringArray[Point1]
            var epsilon_AT@{Index1};
            periods 1;
            values @{ImpulseValues[Point1]};
        @#endfor
    @#else
        @#include "InsertNewShockBlockLines.mod"
    @#endif
end;

options_.qz_criterium = 1 - 1e-8;

steady;
check;

@#if LoadSteadyState
    save_params_and_steady_state( 'SteadyState2.txt' );
@#else
    save_params_and_steady_state( 'SteadyState.txt' );
@#endif

@#if Deterministic
    simul( periods = 10000, maxit = 1000000, tolf = 1e-8, tolx = 1e-8, stack_solve_algo = 7, solve_algo = 0 ); // endogenous_terminal_period
@#else
    stoch_simul( order = 2, pruning, irf = 400, periods = 0, nocorr, nofunctions, nodisplay, nograph, irf_shocks = ( epsilon_AT_2_2, epsilon_GA, epsilon_GN, epsilon_tau, epsilon_phi, epsilon_beta, epsilon_MP ) ); // k_order_solver
@#endif
