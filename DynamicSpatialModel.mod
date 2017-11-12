@#includepath "DynareTransformationEngine"

@#include "Initialize.mod"
@#define UsingGrowthSyntax = 1

@#define SpatialDimensions = 1
@#define SpatialPointsPerDimension = 6
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
psi3 = psi1 * 0.02 / ( 1 - 0.02 ); // http://eyeonhousing.org/2013/01/latest-study-shows-average-buyer-expected-to-stay-in-a-home-13-years/
psi2 = 0.1 - psi1 - psi3;
Gamma = 1;

model;
    @#include "InsertNewModelEquations.mod"
    
    #MaxDistance = 0.5;

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
        #A@{Index1}_LEAD = AP_LEAD * AT@{Index1}_LEAD;
        
        #ZF@{Index1} = ( F@{Index1} / L@{Index1} ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        #SRL@{Index1} = gamma * F@{Index1} / L@{Index1};
        #SP@{Index1} = ( 1 - gamma ) * F@{Index1} / ZF@{Index1};

        #ZF@{Index1}_LEAD = ( F@{Index1}_LEAD / L@{Index1}_LEAD ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
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
        
        muN@{Index1} = beta * ( muN@{Index1}_LEAD * GN_LEAD + U@{Index1}_LEAD ^ ( 1 - varsigma ) + ( 1 - varsigma ) * U@{Index1}_LEAD ^ ( 1 - varsigma ) * ( thetaH * ( H@{Index1}_LEAD / N@{Index1} ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1}_LEAD / N@{Index1} ) ^ ( 1 + nu ) ) + psi1 * SN@{Index1}_LEAD / ( N@{Index1} - SN@{Index1}_LEAD ) - ( thetaC + thetaF + thetaL + psi3 ) ) );
    @#endfor
    
    @#define Index1 = IndicesStringArray[1]
    #Xi_LEAD = beta * ( N@{Index1} / N@{Index1}_LAG ) * ( E@{Index1} / E@{Index1}_LEAD ) * ( U@{Index1}_LEAD / U@{Index1} ) ^ ( 1 - varsigma );

    @#for Point2 in 2 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]

        E@{Index1} / N@{Index1}_LAG / U@{Index1} ^ ( 1 - varsigma ) = E@{Index2} / N@{Index2}_LAG / U@{Index2} ^ ( 1 - varsigma );
    @#endfor

    @#for Point1 in 2 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]

        N@{Index1} = GN * N@{Index1}_LAG - SN@{Index1} + ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + SN@{Index2}@{Index1}
        @#endfor
        ) / @{SpatialNumPoints};
    @#endfor
    
    1 = ( 0
    @#for Point2 in 1 : SpatialNumPoints
        @#define Index2 = IndicesStringArray[Point2]
        + N@{Index2}
    @#endfor
    ) / @{SpatialNumPoints};    

    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        
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
        
        Z@{Index1} = ZF@{Index1} + phi * ( J@{Index1} - ( 1 - deltaJ ) * J@{Index1}_LAG ) + J@{Index1} * ( ( 1 + lambda ) * SP@{Index1} ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar@{Index1};
    @#endfor
    
    1 = R * Xi_LEAD;
    
    0 = ( 0
    @#for Point1 in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point1]
        + E@{Index1} - F@{Index1}
    @#endfor
    ) / @{SpatialNumPoints};
    
end;

@#define Deterministic = 0

shocks;
    @#if Deterministic
        var epsilon_AT_1;
        periods 1:10;
        values 1;
    @#else
        @#include "InsertNewShockBlockLines.mod"
    @#endif
end;

@#define LoadSteadyState = 0

@#if LoadSteadyState

    load_params_and_steady_state( 'SteadyState.txt' );

@#else

    steady_state_model;
        @#include "InsertNewStartSteadyStateEquations.mod"
        
        MaxDistance_ = 0.5;

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

        L_1_ = thetaF * gamma / ( thetaL + thetaF * gamma );

        R_ = 1 / Xi_LEAD_;

        N_1_ = 1;
        N_ = N_1_;

        N_1_LAG_ = 1 / GN_;
        N_LAG_ = N_1_LAG_;

        SN_1_ = psi3 / ( psi1 + psi3 ) * N_1_LAG_;

        SD_1_ = GetSD_1_( psi1, psi2, psi3, N_1_LAG_, SN_1_, MaxDistance_
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            , Distance_1@{Index2}_ 
        @#endfor
        );    

        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            SN_1_@{Index2}_ = psi3 / ( psi1 / ( N_1_LAG_ - SN_1_ ) + psi2 * ( Distance_1@{Index2}_ * SN_1_ - SD_1_ ) / ( MaxDistance_ * SN_1_ * SN_1_ - SN_1_ * SD_1_ ) );
        @#endfor

        A_1_ = 1;
        A_1_LEAD_ = GA_;

        AverageTransportCost_ = ( ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + exp( - tau_ / lambda * Distance_1@{Index2}_ )
        @#endfor
        ) / @{SpatialNumPoints} ) ^ ( - lambda );

        P_1_Over_Q_1_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GQTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;

        F_1_ = GetF_1_( A_1_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ );
        K_1_ = GetK_1_( 1 );
        H_1_ = GetH_1_( 1 );
        Q_1_ = GetQ_1_( 1 );
        
        ZF_1_ = ( F_1_ / L_1_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
        SP_1_ = ( 1 - gamma ) * F_1_ / ZF_1_;
        P_1_ = P_1_Over_Q_1_ * Q_1_;
        Z_1_ = ( ( ( K_1_ / GYTrend_ ) ^ alpha * ( A_1_ * H_1_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_1_ / P_1_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
        SRK_1_ = ( 1 - kappa ) * alpha * SP_1_ * Z_1_ / ( K_1_ / GYTrend_ );
        W_1_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_1_ * Z_1_ / H_1_;
        C_1_ = thetaC * F_1_ / ( thetaF * P_1_ );
        I_1_ = K_1_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
        M_1_ = kappa * SP_1_ * Z_1_ / P_1_;
        Y_1_ = C_1_ + I_1_ + M_1_;
        YBar_1_ = Y_1_ * P_1_ ^ ( ( 1 + lambda ) / lambda ) * AverageTransportCost_ ^ ( - 1 / lambda );
        Pi_1_ = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP_1_ ^ ( - 1 / lambda ) * YBar_1_;
        J_1_ = ( ( 1 + lambda ) * SP_1_ * AverageTransportCost_ / P_1_ ) ^ ( 1 / lambda );
        E_1_ = F_1_;

        U_1_ = ( C_1_ / N_1_LAG_ ) ^ thetaC * ( F_1_ / N_1_LAG_ ) ^ thetaF * ( ( 1 - L_1_ ) / N_1_LAG_ ) ^ thetaL * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H_1_ / N_1_LAG_ ) ^ ( 1 + nu ) ) ^ thetaH * ( 1 - SN_1_ / N_1_LAG_ ) ^ psi1 * ( MaxDistance_ - SD_1_ / SN_1_ ) ^ psi2 * exp( psi3 * ( 0
        @#for Point2 in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + log( SN_1_@{Index2}_ / N_1_LAG_ )
        @#endfor
        ) / @{SpatialNumPoints} );

        U_1_LEAD_ = U_1_ * GUTrend_;
        H_1_LEAD_ = H_1_ * GN_;
        muN_1_ = beta_ * ( U_1_LEAD_ ^ ( 1 - varsigma ) + ( 1 - varsigma ) * U_1_LEAD_ ^ ( 1 - varsigma ) * ( thetaH * ( H_1_LEAD_ / N_1_ ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H_1_LEAD_ / N_1_ ) ^ ( 1 + nu ) ) + psi1 * SN_1_ * GN_ / ( N_1_ - SN_1_ * GN_ ) - ( thetaC + thetaF + thetaL + psi3 ) ) ) / ( 1 - beta_ * GmuNTrend_ * GN_ );
        
        @#for Point1 in 2 : SpatialNumPoints
            @#define Index1 = IndicesStringArray[Point1]
            C@{Index1}_ = C_1_;
            K@{Index1}_ = K_1_;
            I@{Index1}_ = I_1_;
            E@{Index1}_ = E_1_;
            F@{Index1}_ = F_1_;
            Q@{Index1}_ = Q_1_;
            J@{Index1}_ = J_1_;
            L@{Index1}_ = L_1_;
            H@{Index1}_ = H_1_;
            N@{Index1}_ = N_1_;
            SN@{Index1}_ = SN_1_;
            SD@{Index1}_ = SD_1_;
            muN@{Index1}_ = muN_1_;
            U@{Index1}_ = U_1_;
        @#endfor
        
        @#include "InsertNewEndSteadyStateEquations.mod"
    end;

@#endif

options_.qz_criterium = 1 - 1e-8;

steady;
// check;

@#if LoadSteadyState
    save_params_and_steady_state( 'SteadyState2.txt' );
@#else
    save_params_and_steady_state( 'SteadyState.txt' );
@#endif

@#if Deterministic
    simul( periods = 1000, maxit = 1000000, tolf = 1e-8, tolx = 1e-8 ); // endogenous_terminal_period
@#else
    stoch_simul( order = 1, irf = 0, periods = 10000, nocorr, nofunctions ); // k_order_solver
@#endif
