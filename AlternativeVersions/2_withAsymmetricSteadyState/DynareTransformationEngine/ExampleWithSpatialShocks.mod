@#include "Initialize.mod"

// Options determining the spatial topology
@#define SpatialDimensions = 1
@#define SpatialShape = "Torus"
// The other option is "Plane"

@#if SpatialDimensions == 1
    @#define SpatialNorm = "1"
    @#define SpatialPointsPerDimension = 100
@#else
    @#define SpatialNorm = "2"
    @#define SpatialPointsPerDimension = 10
@#endif

// This is an option specific to this file. Setting it to 1 turns on spatial diffusion of technology here.
@#define WithDiffusion = 1

@#if WithDiffusion
    @#if SpatialShape[1] == "P"
        @#define SpatialDiffusionShockProcesses = SpatialDiffusionShockProcesses + [ "A", "0", "Inf", "1", "rho", "sigma", "chi", "exp(-eta*@)#", "exp(-zeta*@)#" ]
        // In order, these are the variable name, its minimum, its maximum, its steady-state, its persistence, its standard deviation, the amount of diffusion, the function governing diffusion (with "@" representing the input distance, and "#" at the end of the string) and the function governing correlation in the shock (with "@" representing the input distance, and "#" at the end of the string)
    @#else
        @#define SpatialDiffusionShockProcesses = SpatialDiffusionShockProcesses + [ "A", "0", "Inf", "1", "rho", "sigma", "chi", "(exp(-eta*@+eta*dBar)+exp(eta*@-eta*dBar))/(exp(eta*dBar)+exp(-eta*dBar))#", "(exp(-zeta*@+zeta*dBar)+exp(zeta*@-zeta*dBar))/(exp(zeta*dBar)+exp(-zeta*dBar))#" ]
        // In order, these are the variable name, its minimum, its maximum, its steady-state, its persistence, its standard deviation, the amount of diffusion, the function governing diffusion (with "@" representing the input distance, and "#" at the end of the string) and the function governing correlation in the shock (with "@" representing the input distance, and "#" at the end of the string)
    @#endif
@#else
    @#if SpatialShape[1] == "P"
        @#define SpatialShockProcesses = SpatialShockProcesses + [ "A", "0", "Inf", "1", "rho", "sigma", "exp(-zeta*@)#" ]
        // In order, these are the variable name, its minimum, its maximum, its steady-state, its persistence, its standard deviation, and the function governing correlation in the shock (with "@" representing the input distance, and "#" at the end of the string)
    @#else
        @#define SpatialShockProcesses = SpatialShockProcesses + [ "A", "0", "Inf", "1", "rho", "sigma", "(exp(-zeta*@+zeta*dBar)+exp(zeta*@-zeta*dBar))/(exp(zeta*dBar)+exp(-zeta*dBar))#" ]
        // In order, these are the variable name, its minimum, its maximum, its steady-state, its persistence, its standard deviation, and the function governing correlation in the shock (with "@" representing the input distance, and "#" at the end of the string)
    @#endif
@#endif

@#include "CreateShocks.mod"

@#define EndoVariables = EndoVariables + [ "R", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "A", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "C", "0", "Inf" ]

@#for Point in 1 : SpatialNumPoints
    @#define Index = IndicesStringArray[Point]
    @#define EndoVariables = EndoVariables + [ "C" + Index, "0", "Inf" ]
    @#define EndoVariables = EndoVariables + [ "B" + Index, "-Inf", "Inf" ]
@#endfor

@#include "ClassifyDeclare.mod"

parameters alpha beta nu rho chi eta zeta sigma phi dBar;

alpha = 0.3;
beta = 0.99;
nu = 2;
rho = 0.95;
chi = 0.5;
eta = 8;
zeta = 4;
sigma = 0.02;
phi = 1e-6;

@#if SpatialShape[1] == "P"
    dBar = ( @{SpatialDimensions} ) ^ ( 1 / ( @{SpatialNorm} ) );
@#else
    dBar = ( @{SpatialDimensions} * ( 0.5 ^ ( @{SpatialNorm} ) ) ) ^ ( 1 / ( @{SpatialNorm} ) );
@#endif

model;
    @#include "InsertNewModelEquations.mod"

    @#for Point in 1 : SpatialNumPoints
        @#define Index = IndicesStringArray[Point]
        #L@{Index} = ( ( 1 - alpha ) * A@{Index} ^ alpha / C@{Index} ) ^ ( 1 / ( alpha + nu ) );
        #Y@{Index} = A@{Index} ^ alpha * L@{Index} ^ ( 1 - alpha );
        1 + phi * B@{Index} = beta * R * C@{Index} / C@{Index}_LEAD;
        C@{Index} + B@{Index} + phi / 2 * B@{Index} ^ 2 = Y@{Index} + R_LAG * B@{Index}_LAG;
    @#endfor

    #B =
    @#for Point in 1 : SpatialNumPoints
        @#define Index = IndicesStringArray[Point]
        + Weight@{Index} * B@{Index}
    @#endfor
    ;

    B = 0;

    log( A ) =
    @#for Point in 1 : SpatialNumPoints
        @#define Index = IndicesStringArray[Point]
        + Weight@{Index} * log( A@{Index} )
    @#endfor
    ;

    C =
    @#for Point in 1 : SpatialNumPoints
        @#define Index = IndicesStringArray[Point]
        + Weight@{Index} * C@{Index}
    @#endfor
    ;
end;

shocks;
    @#include "InsertNewShockBlockLines.mod"
end;

steady_state_model;
    R_ = 1 / beta;
    C_ = ( 1 - alpha ) ^ ( ( 1 - alpha ) / ( 1 + nu ) );
    A_ = 1;
    
    @#for Point in 1 : SpatialNumPoints
        @#define Index = IndicesStringArray[Point]
        C@{Index}_ = C_;
        B@{Index}_ = 0;
    @#endfor
    
    @#include "InsertNewSteadyStateEquations.mod"
end;

options_.qz_criterium = 1 - 1e-8;

steady;
check;

stoch_simul( order = 2, irf = 0, periods = 10000, nocorr, nofunctions, nodecomposition, nomoments, nodisplay, nograph, k_order_solver );
