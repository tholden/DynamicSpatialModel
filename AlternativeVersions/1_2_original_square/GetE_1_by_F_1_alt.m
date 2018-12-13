function E_1_by_F_1_ = GetE_1_by_F_1_alt( Index, SpatialPointsPerDimension , psi3 , psi1 , GN_ )
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    E_by_F_0 = 1 * ones(SpatialPoints,1); 
    N_0 = 1 * ones(SpatialPoints,1); 
    F_0 = 0.0883 * ones(SpatialPoints,1); 
    K_0 = 0.5574 * ones(SpatialPoints,1); 
    H_0 = 0.6248 * ones(SpatialPoints,1); 
    Q_0 = 5.5109 * ones(SpatialPoints,1);  
    SN_0 = ( psi3 / ( psi1 + psi3 ) / GN_ ) * ones( SpatialPoints*SpatialPoints , 1);   
    

initial_values = log( [ E_by_F_0' , N_0' , F_0' , K_0' , H_0' , Q_0' , SN_0' ] );
H0 = 1e-2*eye(length(initial_values));  %Initial Hessian 
crit = 1e-12;                            %Tolerance
nit = 100000;                             %Number of iterations
E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ = exp( csminwel( @GetResidual ,initial_values,H0,[],crit,nit )); 

global E_by_F_x N_x F_x K_x H_x Q_x SNxx_
    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints )';
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints )';
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints )';
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints )';
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints )';
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints )';
    SNxx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end )';
    
    SNxx_ = reshape(SNxx_,[SpatialPoints,SpatialPoints]);
    
    E_1_by_F_1_ = E_by_F_x( Index );
end

function Residual = GetResidual( E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ )
   
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
thetaF = 1; % food off-premises, food services + clothing = about 20% of ( PCE minus housing ) https://www.bea.gov/iTable/iTable.cfm?reqid=19&step=2#reqid=19&step=3&isuri=1&1910=x&0=-9&1921=survey&1903=65&1904=2015&1905=2017&1906=a&1911=0
thetaL = 0.25 / 0.75 * thetaF * gamma; % target of 0.75 for steady-state land use in agriculture, following data from https://www.ers.usda.gov/data-products/major-land-uses/
thetaH = 4;
thetaN = param_thetaN;
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

    d = getDistanceMatrix( shape, SpatialPointsPerDimension ); % create distance 

    Weight = 1 / SpatialPoints;
    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints )';
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints )';
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints )';
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints )';
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints )';
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints )';
    SNxx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end )';
    SNxx_       = reshape(SNxx_,[SpatialPoints,SpatialPoints]);
    
    A_x = 1;
    N_x_LAG_ = N_x ./ GN_; 
    Px_Over_Q_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GSRKTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;
     
    util_SN = zeros( SpatialPoints,1 );
    SD_x = zeros( SpatialPoints,1 );
    SN_x_out = zeros( SpatialPoints,1 );
    SN_x_in = zeros( SpatialPoints,1 );
    integral_N = zeros( SpatialPoints,1 );
    for x=1:SpatialPoints 
        integral_util_SN = zeros( SpatialPoints,1 );
        integral_SD = zeros( SpatialPoints,1 );
        integral_SN_out = zeros( SpatialPoints,1 );
        integral_SN_in = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_util_SN( tildex ) = Weight * GN_ * N_x_LAG_( tildex ) * log( SNxx_( x , tildex ) / N_x_LAG_( x ) );
            integral_SD( tildex ) = Weight * d( x , tildex ) * SNxx_( x , tildex  );
            integral_SN_out( tildex ) = Weight * ( SNxx_( x , tildex ) );
            integral_SN_in( tildex ) = Weight * ( SNxx_( tildex , x ) );
        end
        util_SN( x ) = sum( integral_util_SN );
        SD_x( x ) = sum( integral_SD );
        SN_x_out( x ) = sum( integral_SN_out );
        SN_x_in( x ) = sum( integral_SN_in );
        integral_N( x ) = Weight * ( N_x( x ) );
    end 
    
    SN_x_ = SN_x_out;
    E_x_ = E_by_F_x .* F_x;    
    L_x_ = thetaF .* gamma ./ ( thetaL .* E_by_F_x + thetaF .* gamma );
    ZF_x_ = ( F_x ./ L_x_ .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) );
    SP_x_ = ( 1 - gamma ) .* F_x ./ ZF_x_;
    P_x_ = Px_Over_Q_ .* Q_x;
    Z_x_ = ( ( ( K_x ./ GYTrend_ ) .^ alpha .* ( A_x .* H_x ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* SP_x_ ./ P_x_ ) .^ kappa ) .^ ( 1 / ( 1 - kappa ) );
    SRK_x_ = ( 1 - kappa ) .* alpha .* SP_x_ .* Z_x_ ./ ( K_x ./ GYTrend_ );
    W_x_ = ( 1 - kappa ) .* ( 1 - alpha ) .* SP_x_ .* Z_x_ ./ H_x;
    C_x_ = thetaC .* F_x ./ ( thetaF .* P_x_ );
    I_x_ = K_x .* ( 1 - ( 1 - deltaK ) ./ GYTrend_ ) ./ ( 1 - Phi2 ./ 2 .* ( GYTrend_ - 1 ) .^ 2 );
    M_x_ = kappa .* SP_x_ .* Z_x_ ./ P_x_;
    Y_x_ = C_x_ + I_x_ + M_x_;
    
    YBar_x_ = zeros( SpatialPoints,1 );
    for x=1:SpatialPoints 
        integral_YBar = zeros( 1,SpatialPoints );
        for tildex=1:SpatialPoints
            integral_YBar( tildex ) = Weight * Y_x_( tildex ) * ( P_x_( tildex ) )^((1+lambda)/lambda) * exp( -tau_ * d( x , tildex ) / lambda );
        end
        YBar_x_( x ) = sum( integral_YBar );
    end
    
    J_x_ = ( Z_x_ - ZF_x_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_x_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_x_ ); 

    for x=1:SpatialPoints 
        integral_Px = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_Px( tildex ) = Weight * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
        end
    end
        
    U_x_ = ( C_x_ ./ N_x_LAG_ ).^thetaC .* ( E_x_ ./ N_x_LAG_ ).^thetaF ...
        .* ( ( 1 - L_x_ )./N_x_LAG_ ).^thetaL ...
        .* ( ( 1 / ( 1 + nu ) ) * Gamma^( 1 + nu ) - ( 1 / ( 1 + nu ) ) .* ( H_x ./ N_x_LAG_ ).^( 1 + nu ) ).^thetaH...
        .* ( .5 * Omega ^ 2 - .5 * ( log( GN_ .* N_x_LAG_ ) ).^2 ).^ thetaN ...
         .* ( 1 - SN_x_ ./ N_x_LAG_ ).^psi1 ...
         .* ( dBar - SD_x ./ SN_x_ ).^ psi2 ...
         .* exp( psi3 .* util_SN );
    
    muN_x_ = beta_ .* ( ( U_x_ .* GUTrend_ ) .^ ( 1 - varsigma ) ) .* ( 1 + ( 1 - varsigma ) .* ( ...
            thetaH .* ( H_x ./ N_x_LAG_ ) .^ ( 1 + nu ) ./ ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x ./ N_x_LAG_ ) .^ ( 1 + nu ) ) ...
            - thetaN .* log( N_x ) ./ ( 1 ./ 2 .* Omega .^ 2 - 1 ./ 2 .* log( N_x ) .^ 2 ) ...
            + psi1 .* GN_ .* SN_x_ ./ ( N_x - GN_ .* SN_x_ ) ...
            - ( thetaC + thetaF + thetaL + psi3 ) ...
        ) ) ./ ( 1 - beta_ .* GmuNTrend_ .* GN_ );
    
    Residual = zeros( SpatialPoints , 6 );
    
    Residual( 1 , 1 ) = sum( E_x_ ) - sum( F_x );
    Residual( 2:end , 1 ) = E_x_( 1:end-1 ) .* N_x_LAG_( 2:end ) .* ( U_x_( 2:end ) ).^( 1 - varsigma ) - E_x_( 2:end ) .* N_x_LAG_( 1:end-1 ) .* ( U_x_( 1:end-1 ) ).^( 1 - varsigma );
    Residual( 1 , 2 ) = sum( integral_N ) - 1;
    Residual( 2:end , 2 ) = SN_x_out( 2:end ) - SN_x_in( 2:end);
    Residual( : , 3 ) = Xi_LEAD_ .* GSRKTrend_ .* SRK_x_ ./ ( 1 - Xi_LEAD_ .* GSRKTrend_ .* ( 1 - deltaK ) ) - Q_x;
    Residual( : , 4 ) = thetaF .* N_x_LAG_ ./ F_x .* W_x_ .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x ./ N_x_LAG_ ).^ ( 1 + nu ) ) - thetaH .* ( H_x ./ N_x_LAG_ ).^ nu;
    Residual( : , 5 ) = lambda ./ ( 1 + lambda ) .* ( 1 + lambda ) .^ ( - 1 ./ lambda ) .* SP_x_.^ ( - 1 ./ lambda ) .* YBar_x_ + ( 1 - deltaJ ).* Xi_LEAD_.* phi_ .* SP_x_ .* GSPTrend_ - phi_.* SP_x_;
    Residual( : , 6 ) = P_x_ - ( 1 + lambda ) .* sum( integral_Px ).^(-lambda);
        
    Residual = reshape(Residual,[1,numel(Residual)]);

    Residual_SN = zeros( SpatialPoints , SpatialPoints); 
    for x=1:SpatialPoints
        for tildex=1:SpatialPoints
            Residual_SN( x , tildex ) = muN_x_( x ) - muN_x_( tildex ) - ( 1 - varsigma ) .* N_x_LAG_( x ) .* ( U_x_( x ) ).^( 1 - varsigma ) .* ( GN_ .* psi3 .* N_x_LAG_( tildex ) ./ SNxx_( x , tildex ) - psi1 ./ ( N_x_LAG_( x ) - SN_x_( x ) ) - psi2 .* ( d( ( x ), tildex ) .* ( SN_x_( x ) ) - SD_x( x ) ) ./ ( dBar .* ( SN_x_( x ) ).^2 - SN_x_( x ) .* SD_x( x ) ) );                              
        end
    end
    
     Residual = [ Residual reshape(Residual_SN,[1,numel(Residual_SN)]) ];
end


        
        
        
        