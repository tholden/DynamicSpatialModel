@#includepath "DynareTransformationEngine"

@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1

@#define SpatialDimensions = 2
@#define SpatialPointsPerDimension = 7
@#define SpatialShape = "Torus"

@#if SpatialDimensions == 1
    @#define SpatialNorm = "1"
@#else
    @#define SpatialNorm = "2"
@#endif

@#if SpatialShape[1] == "P"
    @#define SpatialShockProcesses = SpatialShockProcesses + [ "AT", "0", "Inf", "1", "0.9", "0.01", "exp(-zeta*@)#" ]
    %@#define SpatialShockProcesses = SpatialShockProcesses + [ "Utilde_shock", "0", "Inf", "1", "0", "0.001", "exp(-zeta*@)#" ]
@#else
    @#define SpatialShockProcesses = SpatialShockProcesses + [ "AT", "0", "Inf", "1", "0.9", "0.047", "(exp(-zeta*@+zeta*dBar)+exp(zeta*@-zeta*dBar))/(exp(zeta*dBar)+exp(-zeta*dBar))#" ]
    %@#define SpatialShockProcesses = SpatialShockProcesses + [ "Utilde_shock", "0", "Inf", "1", "0", "0.001", "(exp(-zeta*@+zeta*dBar)+exp(zeta*@-zeta*dBar))/(exp(zeta*dBar)+exp(-zeta*dBar))#" ]

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
    @#define EndoVariables = EndoVariables + [ "Utilde" + Index1, "0", "Inf", "1" ]
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

parameters alpha gamma kappa nu varsigma zeta lambda deltaJ deltaK Phi2 PhiL thetaC thetaF thetaL thetaH thetaN psi1 psi2 psi3 Gamma dBar Omega;
parameters UtilityParamSum;


load('param_vals.mat');

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
Omega = param_Omega; // pop/km^2 for the contiguous US is 41.5. for wyoming it is 2.33 for new jersey it is 470. correspond to abs log ratios of 2.88 and 2.43 respectively.

@#if SpatialShape[1] == "P"
    dBar = ( @{SpatialDimensions} ) ^ ( 1 / ( @{SpatialNorm} ) );
@#else
    dBar = ( @{SpatialDimensions} * ( 0.5 ^ ( @{SpatialNorm} ) ) ) ^ ( 1 / ( @{SpatialNorm} ) );
@#endif


@#for Point1 in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point1]
    parameters UtildeSS@{Index1};
    UtildeSS@{Index1} =  Get_Utilde( @{SpatialPointsPerDimension} , @{Point1} , 1 );
@#endfor

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
    // GSPTrend   = GYTrend ^ ( - gamma / ( 1 + lambda ) );
    // GPiTrend   = GYTrend ^ ( - gamma / ( 1 + lambda ) );
    // GYBarTrend = GYTrend ^ ( - gamma / lambda );
       GUTrend    = GYTrend ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN ^ ( thetaC + thetaF + thetaL );
       GmuNTrend  = GYTrend ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );

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

        #ZF@{Index1} = ( F@{Index1} / L@{Index1} ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        #SRL@{Index1} = gamma * F@{Index1} / L@{Index1};
        #SP@{Index1} = ( 1 - gamma ) * F@{Index1} / ZF@{Index1};

        #ZF@{Index1}_LEAD = ( F@{Index1}_LEAD / L@{Index1}_LEAD ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        #SP@{Index1}_LEAD = ( 1 - gamma ) * F@{Index1}_LEAD / ZF@{Index1}_LEAD;

        K@{Index1} = ( 1 - deltaK ) * K@{Index1}_LAG + ( 1 - Phi2 / 2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) ^ 2 ) * I@{Index1};

        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            #SN@{Index1}@{Index2} = psi3 * N@{Index2}_LAG / N_LAG / ( ( muN@{Index1} - muN@{Index2} ) / ( ( 1 - varsigma ) * N@{Index1}_LAG * U@{Index1} ^ ( 1 - varsigma ) ) + psi1 / ( N@{Index1}_LAG - SN@{Index1} ) + psi2 * ( Distance@{Index1}@{Index2} * SN@{Index1} - SD@{Index1} ) / ( dBar * SN@{Index1} * SN@{Index1} - SN@{Index1} * SD@{Index1} ) );
        @#endfor

[name = 'Migration into @{Index1} (point @{Point1})']
        SN@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * SN@{Index1}@{Index2}
        @#endfor
        ;

        SD@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * Distance@{Index1}@{Index2} * SN@{Index1}@{Index2}
        @#endfor
        ;

[name = 'Utility @{Index1} (point @{Point1})']
        U@{Index1} = ( C@{Index1} / N@{Index1}_LAG ) ^ thetaC
        * ( E@{Index1} / N@{Index1}_LAG ) ^ thetaF
        * ( ( 1 - L@{Index1} ) / N@{Index1}_LAG ) ^ thetaL
        * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1} / N@{Index1}_LAG ) ^ ( 1 + nu ) ) ^ thetaH
        * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1}_LAG / N_LAG ) ^ 2 ) ^ thetaN
        * ( 1 - SN@{Index1} / N@{Index1}_LAG ) ^ psi1
        * ( dBar - SD@{Index1} / SN@{Index1} ) ^ psi2
        * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * N@{Index2}_LAG / N_LAG * log( SN@{Index1}@{Index2} / N@{Index1}_LAG )
        @#endfor
        ) ) * Utilde@{Index1};

        %Utilde@{Index1} = Utilde@{Index1}_LAG * Utilde_shock@{Index1};
        Utilde@{Index1} = UtildeSS@{Index1};

        muN@{Index1} = beta * ( muN@{Index1}_LEAD * GN_LEAD + U@{Index1}_LEAD ^ ( 1 - varsigma ) + ( 1 - varsigma ) * U@{Index1}_LEAD ^ ( 1 - varsigma ) * (
            thetaH * ( H@{Index1}_LEAD / N@{Index1} ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_LEAD / N@{Index1} ) ^ ( 1 + nu ) )
            - thetaN * log( N@{Index1} / N ) / ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1} / N ) ^ 2 )
            + psi1 * SN@{Index1}_LEAD / ( N@{Index1} - SN@{Index1}_LEAD )
            - ( thetaC + thetaF + thetaL + psi3 )
        ) );
    @#endfor

    @#define Index1 = IndicesStringArray[1]
    #Xi_LEAD = beta * ( N@{Index1} / N@{Index1}_LAG ) * ( E@{Index1} / E@{Index1}_LEAD ) * ( U@{Index1}_LEAD / U@{Index1} ) ^ ( 1 - varsigma );

    @#for Point2 in 2 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]

[name = 'Marginal Utility equalization between  @{Index1} and @{Index2} (points @{Point1},@{Point2})']
        E@{Index1} / N@{Index1}_LAG / U@{Index1} ^ ( 1 - varsigma ) = E@{Index2} / N@{Index2}_LAG / U@{Index2} ^ ( 1 - varsigma );
    @#endfor

    @#for Point1 in 2 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]


[name = 'Population law of motion  @{Index1} (point @{Point1})']
        N@{Index1} = GN * N@{Index1}_LAG - SN@{Index1}
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * SN@{Index2}@{Index1}
        @#endfor
        ;
    @#endfor

[name = 'total population']
    1 = 0
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        + Weight@{Index2} * N@{Index2}
    @#endfor
    ;

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        #P@{Index1} = ( 1 + lambda ) * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * J@{Index2} * ( SP@{Index2} * exp( tau * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda )
        @#endfor
        ) ^ ( - lambda );

        #P@{Index1}_LEAD = ( 1 + lambda ) * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * J@{Index2}_LEAD * ( SP@{Index2}_LEAD * exp( tau * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda )
        @#endfor
        ) ^ ( - lambda );

        #Z@{Index1} = ( ( K@{Index1}_LAG ^ alpha * ( A@{Index1} * H@{Index1} ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1} / P@{Index1} ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        #Z@{Index1}_LEAD = ( ( K@{Index1} ^ alpha * ( A@{Index1}_LEAD * H@{Index1}_LEAD ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1}_LEAD / P@{Index1}_LEAD ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        #M@{Index1} = kappa * SP@{Index1} * Z@{Index1} / P@{Index1};
        #SRK@{Index1} = ( 1 - kappa ) * alpha * SP@{Index1} * Z@{Index1} / K@{Index1}_LAG;
        #SRK@{Index1}_LEAD = ( 1 - kappa ) * alpha * SP@{Index1}_LEAD * Z@{Index1}_LEAD / K@{Index1};
        #W@{Index1} = ( 1 - kappa ) * ( 1 - alpha ) * SP@{Index1} * Z@{Index1} / H@{Index1};

[name = 'Capital decision (Euler) @{Index1} (point @{Point1})']
        1 = Xi_LEAD * ( SRK@{Index1}_LEAD + Q@{Index1}_LEAD * ( 1 - deltaK ) ) / Q@{Index1};
[name = 'Q - price of capital @{Index1} (point @{Point1})']
        P@{Index1} = Q@{Index1} * ( 1 - Phi2 / 2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) ^ 2 - Phi2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) * I@{Index1} / I@{Index1}_LAG ) + Xi_LEAD * Q@{Index1}_LEAD * Phi2 * ( I@{Index1}_LEAD / I@{Index1} - 1 ) * ( I@{Index1}_LEAD / I@{Index1} ) ^ 2;
        #Y@{Index1} = C@{Index1} + I@{Index1} + M@{Index1};
        thetaC * E@{Index1} = thetaF * P@{Index1} * C@{Index1};
        thetaL * E@{Index1} = thetaF * ( SRL@{Index1} - PhiL * ( L@{Index1} / L@{Index1}_LAG - 1 ) / L@{Index1}_LAG + Xi_LEAD * PhiL * ( L@{Index1}_LEAD / L@{Index1} - 1 ) * ( L@{Index1}_LEAD / ( L@{Index1}^2 ) ) ) * ( 1 - L@{Index1} );
[name = 'Labour supply @{Index1} (point @{Point1})']
        thetaH * ( H@{Index1} / N@{Index1}_LAG ) ^ nu = thetaF * N@{Index1}_LAG / E@{Index1} * W@{Index1} * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1} / N@{Index1}_LAG ) ^ ( 1 + nu ) );
    @#endfor

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        #YBar@{Index1} = 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Weight@{Index2} * Y@{Index2} * P@{Index2} ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance@{Index1}@{Index2} )
        @#endfor
        ;

        #Pi@{Index1} = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP@{Index1} ^ ( - 1 / lambda ) * YBar@{Index1};

[name = 'Firm entry @{Index1} (point @{Point1})']
        phi * SP@{Index1} = Pi@{Index1} + ( 1 - deltaJ ) * Xi_LEAD * phi_LEAD * SP@{Index1}_LEAD;

 [name = 'Raw good market clearing @{Index1} (point @{Point1})']
       Z@{Index1} = ZF@{Index1} + phi * ( J@{Index1} - ( 1 - deltaJ ) * J@{Index1}_LAG ) + J@{Index1} * ( ( 1 + lambda ) * SP@{Index1} ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar@{Index1};
    @#endfor

    1 = R * Xi_LEAD;

    0 = 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + Weight@{Index1} * ( E@{Index1} - F@{Index1} )
    @#endfor
    ;

    @#for VariableName in AggregatedVariables
        @{VariableName} = 0
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
        // GYBarTrend_ = GYTrend_ ^ ( - gamma / lambda );
           GUTrend_    = GYTrend_ ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN_ ^ ( thetaC + thetaF + thetaL );
           GmuNTrend_  = GYTrend_ ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN_ ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );

           Xi_LEAD_ = beta_ * GN_ / GFTrend_ * GUTrend_ ^ ( 1 - varsigma );
           R_ = 1 / Xi_LEAD_;
           N_ = 1;
           N_LAG_ = 1/GN_;

           @#for Point1 in 1 : SpatialNumPoints
               @#define Index1 = IndicesStringArray[Point1]
                Utilde@{Index1}_ = UtildeSS@{Index1};
           @#endfor

           //E_by_F_1_ = GetE_1_by_F_1_homotopy( 1 , @{SpatialPointsPerDimension},  GN_ , nu, gamma, Gamma, Omega, 
           //    lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_,
           //    GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1,
           //    psi2, psi3, tau_, dBar, beta_, varsigma );
           E_by_F_1_ = GetE_1_by_F_1( 1 , @{SpatialPointsPerDimension},  GN_ , nu, gamma, Gamma, Omega, 
               lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_,
               GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1,
               psi2, psi3, tau_, dBar, beta_, varsigma );


           @#for Point1 in 1 : SpatialNumPoints
               @#define Index1 = IndicesStringArray[Point1]
               E_by_F_@{Index1}_ = GetE_x_by_F_x( @{Point1} );
               N@{Index1}_ = GetN_x(  @{Point1} );
               F@{Index1}_ = GetF_x(  @{Point1} );
               K@{Index1}_ = GetK_x(  @{Point1} );
               H@{Index1}_ = GetH_x(  @{Point1} );
               Q@{Index1}_ = GetQ_x(  @{Point1} );
               SN@{Index1}_ = GetSN_x_(  @{Point1} , @{SpatialPointsPerDimension} );
               SD@{Index1}_ = GetSD_x_(  @{Point1} , @{SpatialPointsPerDimension} );
               A_@{Index1}_ = 1;
               A_@{Index1}_LEAD_ = GA_;
               N_LAG_@{Index1}_ = N@{Index1}_ / GN_;
               N_LEAD_@{Index1}_ = N@{Index1}_ * GN_;
               P_Over_Q_@{Index1}_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GQTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;
               E@{Index1}_  = E_by_F_@{Index1}_ * F@{Index1}_;
               L@{Index1}_  = thetaF * gamma / ( thetaL * E_by_F_@{Index1}_  + thetaF * gamma );
               ZF_@{Index1}_ = ( F@{Index1}_ / L@{Index1}_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
               SP_@{Index1}_ = ( 1 - gamma ) * F@{Index1}_ / ZF_@{Index1}_;
               P_@{Index1}_ = P_Over_Q_@{Index1}_ * Q@{Index1}_;
               Z_@{Index1}_ = ( ( ( K@{Index1}_ / GYTrend_ ) ^ alpha * ( A_@{Index1}_ * H@{Index1}_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_@{Index1}_ / P_@{Index1}_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
               SRK_@{Index1}_ = ( 1 - kappa ) * alpha * SP_@{Index1}_ * Z_@{Index1}_ / ( K@{Index1}_ / GYTrend_ );
               W_@{Index1}_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_@{Index1}_ * Z_@{Index1}_ / H@{Index1}_;
               C@{Index1}_ = thetaC * E@{Index1}_ / ( thetaF * P_@{Index1}_ );
               I@{Index1}_ = K@{Index1}_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
               M_@{Index1}_ = kappa * SP_@{Index1}_ * Z_@{Index1}_ / P_@{Index1}_;
               Y_@{Index1}_ = C@{Index1}_ + I@{Index1}_ + M_@{Index1}_;
       @#endfor
        weight__ = (1 / (@{SpatialPointsPerDimension}^2));
        @#for Point1 in 1 : SpatialNumPoints
              @#define Index1 = IndicesStringArray[Point1]
               YBar_@{Index1}_ = 0
               @#for Point2 in 1 : SpatialNumPoints
                   @#define Index2= IndicesStringArray[Point2]
                   + weight__ * Y_@{Index2}_ * P_@{Index2}_ ^ ( ( 1 + lambda ) / lambda ) * exp( - tau_ / lambda * Distance@{Index1}@{Index2}_ )
               @#endfor
               ;
               @#for Point2 in 1 : SpatialNumPoints
                   @#define Index2 = IndicesStringArray[Point2]
                   SN@{Index1}@{Index2}_ = GetSN_xx(  @{Point1} ,  @{Point2} );
               @#endfor
               J@{Index1}_ = ( Z_@{Index1}_ - ZF_@{Index1}_ ) / ( phi_ * ( 1 - ( 1 - deltaJ ) / GJTrend_ ) + ( ( 1 + lambda ) * SP_@{Index1}_ ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar_@{Index1}_ );
               U@{Index1}_ = ( C@{Index1}_ / N_LAG_@{Index1}_ ) ^ thetaC
               * ( E@{Index1}_ / N_LAG_@{Index1}_ ) ^ thetaF
               * ( ( 1 - L@{Index1}_ ) / N_LAG_@{Index1}_ ) ^ thetaL
               * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_ / N_LAG_@{Index1}_ ) ^ ( 1 + nu ) ) ^ thetaH
               * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N_LAG_@{Index1}_ / N_LAG_ ) ^ 2 ) ^ thetaN
               * ( 1 - SN@{Index1}_ / N_LAG_@{Index1}_ ) ^ psi1
               * ( dBar - SD@{Index1}_ / SN@{Index1}_ ) ^ psi2
               * exp( psi3 * ( 0
               @#for Point2 in 1 : SpatialNumPoints
                   @#define Index2 = IndicesStringArray[Point2]
                   + weight__ * N_LAG_@{Index2}_ / N_LAG_ * log( SN@{Index1}@{Index2}_ / N_LAG_@{Index1}_ )
               @#endfor
               ) ) * Utilde@{Index1}_;
               U_@{Index1}_LEAD_ = U@{Index1}_ * GUTrend_;
               H@{Index1}_LEAD_ = H@{Index1}_ * GN_;
               muN@{Index1}_ = beta_ * ( ( U@{Index1}_ * GUTrend_ ) ^ ( 1 - varsigma ) ) * ( 1 + ( 1 - varsigma ) * (
                       thetaH * ( H@{Index1}_ / N_LAG_@{Index1}_ ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_ / N_LAG_@{Index1}_ ) ^ ( 1 + nu ) )
                       - thetaN * log( N@{Index1}_ ) / ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N@{Index1}_ ) ^ 2 )
                       + psi1 * GN_ * SN@{Index1}_ / ( N@{Index1}_ - GN_ * SN@{Index1}_ )
                       - ( thetaC + thetaF + thetaL + psi3 )
                   ) ) / ( 1 - beta_ * GmuNTrend_ * GN_ );
           @#endfor

           @#for VariableName in AggregatedVariables
               @{VariableName}_ = 0
               @#for Point1 in 1 : SpatialNumPoints
                   @#define Index1 = IndicesStringArray[Point1]
                   + weight__ * @{VariableName}@{Index1}_
               @#endfor
               ;
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
    stoch_simul( order = 1, irf = 400, periods = 0, nocorr, nofunctions, nodisplay, nograph, irf_shocks = ( epsilon_AT_1_1, epsilon_AT_4_4, epsilon_GA, epsilon_GN, epsilon_tau, epsilon_phi, epsilon_beta ) ); // k_order_solver
    %stoch_simul( order = 1, irf = 0, periods = 10000, nocorr, nofunctions, nodisplay, nograph ); // k_order_solver
@#endif
