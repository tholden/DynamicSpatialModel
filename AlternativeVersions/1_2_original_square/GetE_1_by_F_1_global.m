function E_1_by_F_1_ = GetE_1_by_F_1_global( shape , Index, SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma)
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    E_by_F_0 = 1 * ones(SpatialPoints,1); 
    N_0 = 1 * ones(SpatialPoints,1); 
    F_0 = 0.0883 * ones(SpatialPoints,1); 
    K_0 = 0.5574 * ones(SpatialPoints,1); 
    H_0 = 0.6248 * ones(SpatialPoints,1); 
    Q_0 = 5.5109 * ones(SpatialPoints,1);  
    SN_0 = 0.001 * ones( SpatialPoints*SpatialPoints , 1);   
    d = getDistanceMatrix( shape, SpatialPointsPerDimension ); % create distance
%'HessianApproximation', {'lbfgs',30}, 
opts = optimoptions( @fmincon, 'Display', 'iter', 'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false);
problem = createOptimProblem('fmincon',...
    'objective', @( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ) GetResidual( exp( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ), shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma),...
    'x0',log( [ E_by_F_0' , N_0' , F_0' , K_0' , H_0' , Q_0' , SN_0' ]' ),'options',opts );
problem.nonlcon = @( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ) constraints( exp( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ), shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma);
gs = GlobalSearch(MultiStart('FunctionTolerance',1e-8,'UseParallel',true));
E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ = exp( run(gs,problem) );

global E_by_F_x N_x F_x K_x H_x Q_x SNxx_
    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints );
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints );
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints );
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints );
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints );
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints );
    SNxx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end );
    
    SNxx_ = reshape(SNxx_,[SpatialPoints,SpatialPoints]);
    
    E_1_by_F_1_ = E_by_F_x( Index );
end

function Residual = GetResidual( E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_, shape, SpatialPointsPerDimension, d, GN_, nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_, varsigma )
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    if strcmp('P',shape)
        Weight = ones(SpatialPointsPerDimension , SpatialPointsPerDimension);
        for i=1:SpatialPointsPerDimension
            for j=1:SpatialPointsPerDimension
                if i==1 || i==SpatialPointsPerDimension
                    if j==1 || j== SpatialPointsPerDimension
                        Weight( i , j ) = (.5 * (1 / ( SpatialPointsPerDimension - 1 )) ) ^2;
                    else
                        Weight( i , j ) = (.5 * (1 / ( SpatialPointsPerDimension - 1 )) ) * (1 / ( SpatialPointsPerDimension - 1 ));
                    end
                else
                    if j==1 || j== SpatialPointsPerDimension
                        Weight( i , j ) = (.5 * (1 / ( SpatialPointsPerDimension - 1 )) ) * (1 / ( SpatialPointsPerDimension - 1 ));
                    else
                        Weight( i , j ) = (1 / ( SpatialPointsPerDimension - 1 ))^2;
                    end
                end
            end
        end
        Weight = reshape(Weight,[SpatialPoints,1]);
    else
        Weight = ( 1 / ( SpatialPointsPerDimension * SpatialPointsPerDimension ) ) * ones( SpatialPoints , 1);
    end
    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints );
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints );
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints );
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints );
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints );
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints );
    SNxx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end );
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
            integral_util_SN( tildex ) = Weight( tildex ) * GN_ * N_x_LAG_( tildex ) * log( SNxx_( x , tildex ) / N_x_LAG_( x ) );
            integral_SD( tildex ) = Weight( tildex ) * d( x , tildex ) * SNxx_( x , tildex  );
            integral_SN_out( tildex ) = Weight( tildex ) * ( SNxx_( x , tildex ) );
            integral_SN_in( tildex ) = Weight( tildex ) * ( SNxx_( tildex , x ) );
        end
        util_SN( x ) = sum( integral_util_SN );
        SD_x( x ) = sum( integral_SD );
        SN_x_out( x ) = sum( integral_SN_out );
        SN_x_in( x ) = sum( integral_SN_in );
        integral_N( x ) = Weight( x ) * ( N_x( x ) );
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
            integral_YBar( tildex ) = Weight( tildex ) * Y_x_( tildex ) * ( P_x_( tildex ) )^((1+lambda)/lambda) * exp( -tau_ * d( x , tildex ) / lambda );
        end
        YBar_x_( x ) = sum( integral_YBar );
    end
    
    J_x_ = ( Z_x_ - ZF_x_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_x_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_x_ ); 

    for x=1:SpatialPoints 
        integral_Px = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_Px( tildex ) = Weight( tildex ) * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
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
     Residual = sum( Residual.^2 );
end


function [c,ceq] = constraints( E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_, shape, SpatialPointsPerDimension, d, GN_, nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_, varsigma )

    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:9 );
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 9+1:2*9 );
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*9+1:3*9 );
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*9+1:4*9 );
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*9+1:5*9 );
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*9+1:6*9 );
    SN_xx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*9+1:end );
    SN_xx_       = reshape(SN_xx_,[9,9]);
  
  % Population higher closer to the centre / population equal everywhere
%   %c   =   []; % e.g. where x1 - 1 >= x2, c(1) = x2 + 1 - x1
  c( 1 )    =   N_x( 2 ) - N_x( 5 ); 
  c( 2 )    =   N_x( 1 ) - N_x( 2 );
  c( 3 )    =   K_x( 2 ) - K_x( 5 ); 
  c( 4 )    =   K_x( 1 ) - K_x( 2 );
  c( 5 )    =   F_x( 5 ) - F_x( 2 ); 
  c( 6 )    =   F_x( 2 ) - F_x( 1 );
  c( 7 )    =   E_by_F_x( 2 ) - E_by_F_x( 5 );
  c( 8 )    =   E_by_F_x( 1 ) - E_by_F_x( 2 );
  c( 9 )    =   H_x( 2 ) - H_x( 5 ); 
  c( 10 )   =   H_x( 1 ) - H_x( 2 );
  ceq( 1 )  =   N_x( 4 ) - N_x( 2 );
  ceq( 2 )  =   N_x( 6 ) - N_x( 2 );
  ceq( 3 )  =   N_x( 8 ) - N_x( 2 );
  ceq( 4 )  =   N_x( 3 ) - N_x( 1 );
  ceq( 5 )  =   N_x( 7 ) - N_x( 1 );
  ceq( 6 )  =   N_x( 9 ) - N_x( 1 );
  ceq( 7 )  =   K_x( 4 ) - K_x( 2 );
  ceq( 8 )  =   K_x( 6 ) - K_x( 2 );
  ceq( 9 )  =   K_x( 8 ) - K_x( 2 );
  ceq( 10 ) =   K_x( 3 ) - K_x( 1 );
  ceq( 11 ) =   K_x( 7 ) - K_x( 1 );
  ceq( 12 ) =   K_x( 9 ) - K_x( 1 );
  ceq( 13 ) =   F_x( 4 ) - F_x( 2 );
  ceq( 14 ) =   F_x( 6 ) - F_x( 2 );
  ceq( 15 ) =   F_x( 8 ) - F_x( 2 );
  ceq( 16 ) =   F_x( 3 ) - F_x( 1 );
  ceq( 17 ) =   F_x( 7 ) - F_x( 1 );
  ceq( 18 ) =   F_x( 9 ) - F_x( 1 );
  ceq( 19 ) =   E_by_F_x( 4 ) - E_by_F_x( 2 );
  ceq( 20 ) =   E_by_F_x( 6 ) - E_by_F_x( 2 );
  ceq( 21 ) =   E_by_F_x( 8 ) - E_by_F_x( 2 );
  ceq( 22 ) =   E_by_F_x( 3 ) - E_by_F_x( 1 );
  ceq( 23 ) =   E_by_F_x( 7 ) - E_by_F_x( 1 );
  ceq( 24 ) =   E_by_F_x( 9 ) - E_by_F_x( 1 );
  ceq( 25 ) =   H_x( 4 ) - H_x( 2 );
  ceq( 26 ) =   H_x( 6 ) - H_x( 2 );
  ceq( 27 ) =   H_x( 8 ) - H_x( 2 );
  ceq( 28 ) =   H_x( 3 ) - H_x( 1 );
  ceq( 29 ) =   H_x( 7 ) - H_x( 1 );
  ceq( 30 ) =   H_x( 9 ) - H_x( 1 );
  
  % SN_xx - migration higher to closer locations
  % SN_xx - migration across equivalent distances higher for more highly
  % populated locations

  % Square:
  % From location x = 1
%   c( 11 )   =   SN_xx_( 1 , 2 ) - SN_xx_( 1 , 1 );
%   c( 12 )   =   SN_xx_( 1 , 5 ) - SN_xx_( 1 , 2 );
%   c( 13 )   =   SN_xx_( 1 , 3 ) - SN_xx_( 1 , 5 );
%   c( 14 )   =   SN_xx_( 1 , 6 ) - SN_xx_( 1 , 3 );
%   c( 15 )   =   SN_xx_( 1 , 9 ) - SN_xx_( 1 , 6 );
%   ceq( 31 ) =   SN_xx_( 1 , 2 ) - SN_xx_( 1 , 4 );
%   ceq( 32 ) =   SN_xx_( 1 , 3 ) - SN_xx_( 1 , 7 );
%   ceq( 33 ) =   SN_xx_( 1 , 6 ) - SN_xx_( 1 , 8 );
%   
%   % Symmetric location x = 1 , 7 , 3 , 9
%   ceq( 34 ) =   SN_xx_( 1 , 1 ) - SN_xx_( 3 , 3 );
%   ceq( 35 ) =   SN_xx_( 1 , 2 ) - SN_xx_( 3 , 2 );
%   ceq( 36 ) =   SN_xx_( 1 , 3 ) - SN_xx_( 3 , 1 );
%   ceq( 37 ) =   SN_xx_( 1 , 4 ) - SN_xx_( 3 , 6 );
%   ceq( 38 ) =   SN_xx_( 1 , 5 ) - SN_xx_( 3 , 5 );
%   ceq( 39 ) =   SN_xx_( 1 , 6 ) - SN_xx_( 3 , 4 );
%   ceq( 40 ) =   SN_xx_( 1 , 7 ) - SN_xx_( 3 , 9 );
%   ceq( 41 ) =   SN_xx_( 1 , 8 ) - SN_xx_( 3 , 8 );
%   ceq( 42 ) =   SN_xx_( 1 , 9 ) - SN_xx_( 3 , 7 );
%   ceq( 43 ) =   SN_xx_( 1 , 1 ) - SN_xx_( 7 , 7 );
%   ceq( 44 ) =   SN_xx_( 1 , 2 ) - SN_xx_( 7 , 8 );
%   ceq( 45 ) =   SN_xx_( 1 , 3 ) - SN_xx_( 7 , 9 );
%   ceq( 46 ) =   SN_xx_( 1 , 4 ) - SN_xx_( 7 , 4 );
%   ceq( 47 ) =   SN_xx_( 1 , 5 ) - SN_xx_( 7 , 5 );
%   ceq( 48 ) =   SN_xx_( 1 , 6 ) - SN_xx_( 7 , 6 );
%   ceq( 49 ) =   SN_xx_( 1 , 7 ) - SN_xx_( 7 , 1 );
%   ceq( 50 ) =   SN_xx_( 1 , 8 ) - SN_xx_( 7 , 2 );
%   ceq( 51 ) =   SN_xx_( 1 , 9 ) - SN_xx_( 7 , 3 );
%   ceq( 52 ) =   SN_xx_( 1 , 1 ) - SN_xx_( 9 , 9 );
%   ceq( 53 ) =   SN_xx_( 1 , 2 ) - SN_xx_( 9 , 8 );
%   ceq( 54 ) =   SN_xx_( 1 , 3 ) - SN_xx_( 9 , 7 );
%   ceq( 55 ) =   SN_xx_( 1 , 4 ) - SN_xx_( 9 , 6 );
%   ceq( 56 ) =   SN_xx_( 1 , 5 ) - SN_xx_( 9 , 5 );
%   ceq( 57 ) =   SN_xx_( 1 , 6 ) - SN_xx_( 9 , 4 );
%   ceq( 58 ) =   SN_xx_( 1 , 7 ) - SN_xx_( 9 , 3 );
%   ceq( 59 ) =   SN_xx_( 1 , 8 ) - SN_xx_( 9 , 2 );
%   ceq( 60 ) =   SN_xx_( 1 , 9 ) - SN_xx_( 9 , 1 );
%   
%   % From location x = 2
%   c( 16 )   =   SN_xx_( 2 , 1 ) - SN_xx_( 1 , 2 );
%   c( 17 )   =   SN_xx_( 2 , 4 ) - SN_xx_( 1 , 1 );
%   c( 18 )   =   SN_xx_( 2 , 8 ) - SN_xx_( 1 , 4 );
%   c( 19 )   =   SN_xx_( 2 , 7 ) - SN_xx_( 1 , 8 );
%   ceq( 61 ) =   SN_xx_( 2 , 1 ) - SN_xx_( 1 , 3 );
%   ceq( 62 ) =   SN_xx_( 2 , 6 ) - SN_xx_( 1 , 4 );
%   ceq( 63 ) =   SN_xx_( 2 , 9 ) - SN_xx_( 1 , 7 );
%   
%   % Symmetric location x = 2 , 4 , 6 , 8
%   ceq( 64 ) =   SN_xx_( 2 , 1 ) - SN_xx_( 4 , 1 );
%   ceq( 65 ) =   SN_xx_( 2 , 2 ) - SN_xx_( 4 , 4 );
%   ceq( 66 ) =   SN_xx_( 2 , 3 ) - SN_xx_( 4 , 7 );
%   ceq( 67 ) =   SN_xx_( 2 , 4 ) - SN_xx_( 4 , 2 );
%   ceq( 68 ) =   SN_xx_( 2 , 5 ) - SN_xx_( 4 , 5 );
%   ceq( 69 ) =   SN_xx_( 2 , 6 ) - SN_xx_( 4 , 8 );
%   ceq( 70 ) =   SN_xx_( 2 , 7 ) - SN_xx_( 4 , 3 );
%   ceq( 71 ) =   SN_xx_( 2 , 8 ) - SN_xx_( 4 , 6 );
%   ceq( 72 ) =   SN_xx_( 2 , 9 ) - SN_xx_( 4 , 9 );
%   ceq( 73 ) =   SN_xx_( 2 , 1 ) - SN_xx_( 6 , 9 );
%   ceq( 74 ) =   SN_xx_( 2 , 2 ) - SN_xx_( 6 , 6 );
%   ceq( 75 ) =   SN_xx_( 2 , 3 ) - SN_xx_( 6 , 3 );
%   ceq( 76 ) =   SN_xx_( 2 , 4 ) - SN_xx_( 6 , 8 );
%   ceq( 77 ) =   SN_xx_( 2 , 5 ) - SN_xx_( 6 , 5 );
%   ceq( 78 ) =   SN_xx_( 2 , 6 ) - SN_xx_( 6 , 2 );
%   ceq( 79 ) =   SN_xx_( 2 , 7 ) - SN_xx_( 6 , 7 );
%   ceq( 80 ) =   SN_xx_( 2 , 8 ) - SN_xx_( 6 , 4 );
%   ceq( 81 ) =   SN_xx_( 2 , 9 ) - SN_xx_( 6 , 1 );
%   ceq( 82 ) =   SN_xx_( 2 , 1 ) - SN_xx_( 8 , 7 );
%   ceq( 83 ) =   SN_xx_( 2 , 2 ) - SN_xx_( 8 , 8 );
%   ceq( 84 ) =   SN_xx_( 2 , 3 ) - SN_xx_( 8 , 9 );
%   ceq( 85 ) =   SN_xx_( 2 , 4 ) - SN_xx_( 8 , 4 );
%   ceq( 86 ) =   SN_xx_( 2 , 5 ) - SN_xx_( 8 , 5 );
%   ceq( 87 ) =   SN_xx_( 2 , 6 ) - SN_xx_( 8 , 6 );
%   ceq( 88 ) =   SN_xx_( 2 , 7 ) - SN_xx_( 8 , 1 );
%   ceq( 89 ) =   SN_xx_( 2 , 8 ) - SN_xx_( 8 , 2 );
%   ceq( 90 ) =   SN_xx_( 2 , 9 ) - SN_xx_( 8 , 3 );



 % Food consumption / production ratio:
  c( 20 )    =   E_by_F_x( 1 ) - 1;
  c( 21 )    =   1 - E_by_F_x( 5 );
 
  
  % population equal everywhere
 ceq = [ceq c];
 c = []; 
 
SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    Weight = getWeightMatrix( shape , SpatialPointsPerDimension );
    Weight = reshape(Weight,[SpatialPoints,1]);
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
            integral_util_SN( tildex ) = Weight( tildex ) * GN_ * N_x_LAG_( tildex ) * log( SN_xx_( x , tildex ) / N_x_LAG_( x ) );
            integral_SD( tildex ) = Weight( tildex ) * d( x , tildex ) * SN_xx_( x , tildex  );
            integral_SN_out( tildex ) = Weight( tildex ) * ( SN_xx_( x , tildex ) );
            integral_SN_in( tildex ) = Weight( tildex ) * ( SN_xx_( tildex , x ) );
        end
        util_SN( x ) = sum( integral_util_SN );
        SD_x( x ) = sum( integral_SD );
        SN_x_out( x ) = sum( integral_SN_out );
        SN_x_in( x ) = sum( integral_SN_in );
        integral_N( x ) = Weight( x ) * ( N_x( x ) );
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
            integral_YBar( tildex ) = Weight( tildex ) * Y_x_( tildex ) * ( P_x_( tildex ) )^((1+lambda)/lambda) * exp( -tau_ * d( x , tildex ) / lambda );
        end
        YBar_x_( x ) = sum( integral_YBar );
    end
    
    J_x_ = ( Z_x_ - ZF_x_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_x_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_x_ ); 

    for x=1:SpatialPoints 
        integral_Px = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_Px( tildex ) = Weight( tildex ) * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
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
            Residual_SN( x , tildex ) = muN_x_( x ) - muN_x_( tildex ) - ( 1 - varsigma ) .* N_x_LAG_( x ) .* ( U_x_( x ) ).^( 1 - varsigma ) .* ( GN_ .* psi3 .* N_x_LAG_( tildex ) ./ SN_xx_( x , tildex ) - psi1 ./ ( N_x_LAG_( x ) - SN_x_( x ) ) - psi2 .* ( d( ( x ), tildex ) .* ( SN_x_( x ) ) - SD_x( x ) ) ./ ( dBar .* ( SN_x_( x ) ).^2 - SN_x_( x ) .* SD_x( x ) ) );                              
        end
    end
    
     Residual = [ Residual reshape(Residual_SN,[1,numel(Residual_SN)]) ];
     
    % ceq = [ceq Residual ];
 
end        
        
        