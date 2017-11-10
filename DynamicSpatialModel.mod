@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1

@#define SpatialDimensions = 1
@#define SpatialPointsPerDimension = 10
@#define SpatialShape = "Torus"
@#define SpatialNorm = "1"

@#define SpatialShockProcesses = SpatialShockProcesses + [ "AT", "0", "Inf", "1", "0.9", "0.001", "(exp(-zeta*@+zeta/2)+exp(zeta*@-zeta/2))/(exp(zeta/2)+exp(-zeta/2))#" ]

@#define ShockProcesses = ShockProcesses + [ "GA", "0", "Inf", "1.005", "0.8", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "GN", "0", "Inf", "1.0025", "0.5", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "tau", "0", "Inf", "2", "0.95", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "phi", "0", "Inf", "1", "0.95", "0.001" ]
@#define ShockProcesses = ShockProcesses + [ "beta", "0", "1", "0.99", "0.95", "0.001" ]

@#define PureTrendEndoVariables = PureTrendEndoVariables + [ "AP", "GA" ]

@#include "CreateShocks.mod"

@#define EndoVariables = EndoVariables + [ "R", "0", "Inf", 1 ]

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

@#include "ClassifyDeclare.mod"

parameters alpha gamma kappa nu varsigma zeta lambda deltaJ deltaK Phi2 thetaC thetaF thetaL thetaH psi1 psi2 psi3 Gamma;

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
thetaC = 0.15;
thetaF = 0.15;
thetaL = 0.3;
thetaH = 0.3;
psi1 = 0.05;
psi2 = 0.025;
psi3 = 0.025;
Gamma = 1;

model;
    @#include "InsertNewModelEquations.mod"
    
    #MaxDistance = 0.5;

    #  GYTrend    = ( GA * GN ) ^ ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) / ( ( 1 - alpha ) * ( 1 - kappa ) * ( 1 + lambda ) - lambda ) );
    #  GZTrend    = GYTrend ^ ( 1 / ( 1 + lambda ) );
    // GJTrend    = GYTrend ^ ( 1 / ( 1 + lambda ) );
    #  GFTrend    = GYTrend ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
    // GSRLTrend  = GYTrend ^ ( ( 1 - gamma ) / ( 1 + lambda ) );
    // GWTrend    = GYTrend ^ ( ( 1 - gamma ) / ( 1 + lambda ) ) / GN;
    #  GPTrend    = GYTrend ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
    // GQTrend    = GYTrend ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
    // GSRKTrend  = GYTrend ^ ( - ( gamma + lambda ) / ( 1 + lambda ) );
    // GSPTrend   = GYTrend ^ ( - gamma / ( 1 + lambda ) );
    // GPiTrend   = GYTrend ^ ( - gamma / ( 1 + lambda ) );
    // GYBarTrend = GYTrend ^ ( - gamma / lambda );
    #  GUTrend    = GYTrend ^ ( thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) / GN ^ ( thetaC + thetaF + thetaL );
    #  GmuNTrend  = GYTrend ^ ( (  thetaC + thetaF * ( 1 - gamma ) / ( 1 + lambda ) ) * ( 1 - varsigma ) ) / GN ^ ( ( thetaC + thetaF + thetaL ) * ( 1 - varsigma ) );
    
    #N = ( 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + N@{Index1}
    @#endfor
    ) / @{SpatialNumPoints};

    #N_LAG = ( 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + N@{Index1}_LAG
    @#endfor
    ) / @{SpatialNumPoints};
    
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        
        #A@{Index1} = AP * AT@{Index1};
        
        #ZF@{Index1} = ( F@{Index1} / L@{Index1} ^ gamma ) ^ ( 1 / ( 1 - gamma );
        #SRL@{Index1} = gamma * F@{Index1} / L@{Index1};
        #SP@{Index1} = ( 1 - gamma ) * F@{Index1} / ZF@{Index1};

        #ZF@{Index1}_LEAD = ( F@{Index1}_LEAD / L@{Index1}_LEAD ^ gamma ) ^ ( 1 / ( 1 - gamma );
        #SP@{Index1}_LEAD = ( 1 - gamma ) * F@{Index1}_LEAD / ZF@{Index1}_LEAD;
        
        K@{Index1} = ( 1 - deltaK ) * K@{Index1}_LAG + ( 1 - Phi2 / 2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) ^ 2 ) * I@{Index1};
        
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            #SN@{Index1}@{Index2} = psi3 * N@{Index2}_LAG / N_LAG / ( ( muN@{Index1} - muN@{Index2} ) / ( ( 1 - varsigma ) * N@{Index1}_LAG * U@{Index1} ^ ( 1 - varsigma ) ) + psi1 / ( N@{Index1}_LAG - SN@{Index1} ) + psi2 * ( Distance@{Index1}@{Index2} * SN@{Index1} - SD@{Index1} ) / ( MaxDistance * SN@{Index1} * SN@{Index1} - SN@{Index1} * SD@{Index1} ) );
        @#endfor
        
        SN@{Index1} = ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + SN@{Index1}@{Index2}
        @#endfor
        ) / @{SpatialNumPoints};

        SD@{Index1} = ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Distance@{Index1}@{Index2} * SN@{Index1}@{Index2}
        @#endfor
        ) / @{SpatialNumPoints};
        
        U@{Index1} = ( C@{Index1} / N@{Index1}_LAG ) ^ thetaC * ( E@{Index1} / N@{Index1}_LAG ) ^ thetaF * ( ( 1 - L@{Index1} ) / N@{Index1}_LAG ) ^ thetaL * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1} / N@{Index1}_LAG ) ^ ( 1 + nu ) ) ^ thetaH * ( 1 - SN@{Index1} / N@{Index1}_LAG ) ^ psi1 * ( MaxDistance - SD@{Index1} / SN@{Index1} ) ^ psi2 * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + N@{Index2}_LAG / N_LAG * log( SN@{Index1}@{Index2} / N@{Index1}_LAG )
        @#endfor
        ) / @{SpatialNumPoints} );
        
        muN@{Index1} = beta * ( muN@{Index1}_LEAD * GN@{Index1}_LEAD + U@{Index1}_LEAD ^ ( 1 - varsigma ) + ( 1 - varsigma ) * U@{Index1}_LEAD ^ ( 1 - varsigma ) * ( thetaH * ( H@{Index1}_LEAD / N@{Index1} ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_LEAD / N@{Index1} ) ^ ( 1 + nu ) ) + psi1 * SN@{Index1}_LEAD / ( N@{Index1} - SN@{Index1}_LEAD ) - ( thetaC + thetaF + thetaL + psi3 ) ) );
    @#endfor
    
    @#define Index1 = IndicesStringArray[1]
    #Xi_LEAD = beta * ( N@{Index1} / N@{Index1}_LAG ) * ( E@{Index1} / E@{Index1}_LEAD ) * ( U@{Index1}_LEAD / U@{Index1} ) ^ ( 1 - varsigma );

    @#for Point2 in 2 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]

        E@{Index1} / N@{Index1}_LAG / U@{Index1} ^ ( 1 - varsigma ) = E@{Index2} / N@{Index2}_LAG / U@{Index2} ^ ( 1 - varsigma );
    @#endfor

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        N@{Index1} = GN * N@{Index1}_LAG - SN@{Index1} + ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + SN@{Index2}@{Index1}
        @#endfor
        ) / @{SpatialNumPoints};
        
        #P@{Index1} = ( 1 + lambda ) * ( ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + J@{Index2} * ( SP@{Index2} * exp( tau * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda )
        @#endfor
        ) / @{SpatialNumPoints} ) ^ ( - lambda );
        
        #P@{Index1}_LEAD = ( 1 + lambda ) * ( ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + J@{Index2}_LEAD * ( SP@{Index2}_LEAD * exp( tau * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda )
        @#endfor
        ) / @{SpatialNumPoints} ) ^ ( - lambda );
        
        #Z@{Index1} = ( ( K@{Index1}_LAG ^ alpha * ( A@{Index1} * H@{Index1} ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1} / P@{Index1} ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        #Z@{Index1}_LEAD = ( ( K@{Index1} ^ alpha * ( A@{Index1}_LEAD * H@{Index1}_LEAD ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1}_LEAD / P@{Index1}_LEAD ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        #M@{Index1} = kappa * SP@{Index1} * Z@{Index1} / P@{Index1};
        #SRK@{Index1} = ( 1 - kappa ) * alpha * SP@{Index1} * Z@{Index1} / K@{Index1}_LAG;
        #SRK@{Index1}_LEAD = ( 1 - kappa ) * alpha * SP@{Index1}_LEAD * Z@{Index1}_LEAD / K@{Index1};
        #W@{Index1} = ( 1 - kappa ) * ( 1 - alpha ) * SP@{Index1} * Z@{Index1} / H@{Index1};
        
        1 = Xi_LEAD * ( SRK@{Index1}_LEAD + Q@{Index1}_LEAD * ( 1 - deltaK ) ) / Q@{Index1};
        P@{Index1} = Q@{Index1} * ( 1 - Phi2 / 2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) ^ 2 - Phi2 * ( I@{Index1} / I@{Index1}_LAG - 1 ) * I@{Index1} / I@{Index1}_LAG ) + Xi_LEAD * Q@{Index1}_LEAD * Phi2 * ( I@{Index1}_LEAD / I@{Index1} - 1 ) * ( I@{Index1}_LEAD / I@{Index1} ) ^ 2;
        #Y@{Index1} = C@{Index1} + I@{Index1} + M@{Index1};
        thetaC * E@{Index1} = thetaF * P@{Index1} * C@{Index1};
        thetaL * E@{Index1} = thetaF * SRL@{Index1} * ( 1 - L@{Index1} );
        thetaH * ( H@{Index1} / N@{Index1}_LAG ) ^ nu = thetaF * N@{Index1} / E@{Index1} * W@{Index1} * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1} / N@{Index1}_LAG ) ^ ( 1 + nu ) );
    @#endfor

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        
        #YBar@{Index1} = ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + Y@{Index2} * P@{Index2} ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance@{Index1}@{Index2} )
        @#endfor
        ) / @{SpatialNumPoints};
        
        #Pi@{Index1} = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP@{Index1} ^ ( - 1 / lambda ) * YBar@{Index1};
        
        phi * SP@{Index1} = Pi@{Index1} + ( 1 - deltaJ ) * Xi_LEAD * phi_LEAD * SP@{Index1}_LEAD;
        
        Z@{Index1} = ZF@{Index1} + phi * ( J@{Index1} - ( 1 - deltaJ ) * J@{Index1}_LAG ) + J@{Index1} * ( ( 1 + lambda ) * SP@{Index1} ) ^ ( - ( 1 + lambda ) / lambda ) *  YBarSP@{Index1};
    @#endfor
    
    1 = R * Xi_LEAD;
    
    0 = ( 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + E@{Index1} - F@{Index1}
    @#endfor
    ) / @{SpatialNumPoints};
    
end;

shocks;
    @#include "InsertNewShockBlockLines.mod"
end;

//steady_state_model;
//    R_ = 1 / beta;
//    @#for Point1 in 1 : SpatialNumPoints
//        @#define Index1 = IndicesStringArray[Point1]
//        C@{Index1}_ = ( 1 - alpha ) ^ ( ( 1 - alpha ) / ( 1 + nu ) );
//    @#endfor
//    @#include "InsertNewSteadyStateEquations.mod"
//end;

steady;
//check;

stoch_simul( order = 1, irf = 0, periods = 1000, nocorr, nofunctions ) log_R;
