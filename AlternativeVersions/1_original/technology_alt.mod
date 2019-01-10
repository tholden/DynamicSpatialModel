@#includepath "DynareTransformationEngine"
@#include "Initialize.mod"
@#define SpatialNumPoints = 49
@#define SpatialDimensions = 1
@#define SpatialPointsPerDimension = SpatialNumPoints / SpatialDimensions
@#define Indices = EmptyNumericArray
@#define IndicesArray = EmptyNumericArray
@#define IndicesStringArray = EmptyStringArray
@#define Depth = 0
@#include "CreateIndicesArray.mod"
@#include "DefineNumbers.mod"


@#define estimation_toggle = 1
@#define stoch_simul_toggle = 0

@#for Point in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point]
    var a@{Index1} da@{Index1};
    varexo eps_a@{Index1};
@#endfor

var a;
varexo eps_a;

parameters rho_a rho_ax zetta sigma_a sigma_ax;

rho_a = 0.7912;
rho_ax  = 0.8703;
zetta = 19;
sigma_a = 0.0072;
sigma_ax = 0.0043;

% shut down aggregate
//zetta = 5.4;
//sigma_a = 0.000072;
//sigma_ax = 0.0052;
//rho_ax  = 0.8632;

@#for Point1 in 1 : SpatialNumPoints
  @#for Point2  in 1 : SpatialNumPoints
          @#define Index1 = IndicesStringArray[Point1]
          @#define Index2 = IndicesStringArray[Point2]
          @#define Point1_num = Numbers[Point1+1]
          @#define Point2_num = Numbers[Point2+1]
          parameters d@{Index1}@{Index2} w@{Index1}@{Index2} ;
          d@{Index1}@{Index2} = getDistance(@{Point1_num},@{Point2_num});
          w@{Index1}@{Index2} = exp( -zetta * d@{Index1}@{Index2} );
  @#endfor
@#endfor

model;

a = rho_a * a(-1) + sigma_a * eps_a;

@#for Point1 in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point1]
    #u@{Index1} = sigma_ax * ( eps_a@{Index1}
    @#for Point2  in 1 : SpatialNumPoints
            @#define Index2 = IndicesStringArray[Point2]
            + exp( -zetta * d@{Index1}@{Index2} ) * eps_a@{Index2}
    @#endfor
    );
    a@{Index1} = a + rho_ax * ( a@{Index1}(-1) - a(-1) ) + u@{Index1};
    da@{Index1} = a@{Index1} - a@{Index1}(-1);
@#endfor

end;

shocks;
var eps_a; stderr 1;
@#for Point in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point]
    var eps_a@{Index1}; stderr 1;
@#endfor
end;


@#if estimation_toggle

estimated_params;
%---% Bayesian:
//sigma_a,          INV_GAMMA_PDF,0.004, 2;       //technology
//sigma_ax,         INV_GAMMA_PDF,0.004, 2;       //state-technology
//rho_a,            BETA_PDF,     0.5, 0.2;     //AR1 technology
//rho_ax,           BETA_PDF,     0.5, 0.2;     //AR1 state-technology
//zetta,            INV_GAMMA_PDF,   10,   20;       //Spatial correlation
%---% Maximum likelihood:
sigma_a,    0.005;
sigma_ax,   0.005;
rho_a,      0.5;
rho_ax,     0.5;
zetta,      10, 0 , 100;
end;


varobs
@#for Point in 1 : SpatialNumPoints
    @#define Index1 = IndicesStringArray[Point]
     da@{Index1}
@#endfor
;

estimation(datafile='employment_and_productivity.xlsx',xls_sheet=dynare_2,xls_range = B1:AX28, mh_jscale=1
 ,mode_compute=9   
)
    @#for Point in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point]
        da@{Index1} 
    @#endfor
;

@#endif

@#if stoch_simul_toggle
stoch_simul(order=1,irf_shocks = ( eps_a_26 ), nograph)
    @#for Point in 1 : SpatialNumPoints
        @#define Index1 = IndicesStringArray[Point]
        a@{Index1} 
    @#endfor
; % 26 = Nebraska
@#endif