function F_x = GetF_x_(  Index , SpatialPointsPerDimension, A_x_, N_x_, N_x_LAG_, nu, gamma, Omega, Gamma, L_x_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaN, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_x_Over_Q_x_, Xi_LEAD_, deltaK, Phi2, GmuNTrend_, GUTrend_, GSPTrend_, tau_, GN_,  psi1, psi2, psi3, dBar, beta_, varsigma )
    global F_x_K_x_H_x_Q_x_SN_xx_ 
    d = getDistanceMatrix( 'T', SpatialPointsPerDimension );
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    F_0_ = 0.0819 * ones( SpatialPoints , 1 );
    K_0_ = 0.4425 * ones( SpatialPoints , 1 );
    H_0_ = 0.6252 * ones( SpatialPoints , 1 );
    Q_0_ = 6.4416 * ones( SpatialPoints , 1 );
    SN_xx_0_ = 0.001 * ones( SpatialPoints , SpatialPoints );   
    F_x_K_x_H_x_Q_x_SN_xx_ = exp( fsolve( @( log_F_x_K_x_H_x_Q_x_SN_xx_ ) GetResidual( exp( log_F_x_K_x_H_x_Q_x_SN_xx_ ), SpatialPointsPerDimension, d, A_x_, N_x_, N_x_LAG_, nu, gamma, Omega, Gamma, L_x_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaN, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_x_Over_Q_x_, Xi_LEAD_, deltaK, Phi2, GmuNTrend_, GUTrend_, GSPTrend_, tau_, GN_,  psi1, psi2, psi3, dBar, beta_, varsigma  ), ...
        log( [ F_0_ K_0_ H_0_ Q_0_ SN_xx_0_ ] ), ...
        optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
    global F_x_ K_x_ H_x_ Q_x_ SN_xx_
    F_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 1 );
    K_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 2 );
    H_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 3 );
    Q_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 4 );
    SN_xx_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 5:end );
    disp( F_x_K_x_H_x_Q_x_SN_xx_ );
    
    F_x = F_x_( Index );
end

function Residual = GetResidual( F_x_K_x_H_x_Q_x_SN_xx_, SpatialPointsPerDimension, d, A_x_, N_x_, N_x_LAG_, nu, gamma, Omega, Gamma, L_x_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaN, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_x_Over_Q_x_, Xi_LEAD_, deltaK, Phi2, GmuNTrend_, GUTrend_, GSPTrend_ , tau_, GN_, psi1, psi2, psi3, dBar, beta_,varsigma )
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    F_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 1 );
    K_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 2 );
    H_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 3 );
    Q_x_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 4 );
    SN_xx_ = F_x_K_x_H_x_Q_x_SN_xx_( : , 5:end );
    
    N_x_LAG_ = N_x_LAG_ * ones( SpatialPoints,1 );
    N_x = N_x_ * ones( SpatialPoints,1 );
    
    Weight = 1 / SpatialPoints;
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
            integral_util_SN( tildex ) = Weight * GN_ * N_x_LAG_( tildex ) * log( SN_xx_( x , tildex ) / N_x_LAG_( x ) );
            integral_SD( tildex ) = Weight * d( x , tildex ) * SN_xx_( x , tildex  );
            integral_SN_out( tildex ) = Weight * ( SN_xx_( x , tildex ) );
            integral_SN_in( tildex ) = Weight * ( SN_xx_( tildex , x ) );
        end
        util_SN( x ) = sum( integral_util_SN );
        SD_x( x ) = sum( integral_SD );
        SN_x_out( x ) = sum( integral_SN_out );
        SN_x_in( x ) = sum( integral_SN_in );
        integral_N( x ) = Weight * ( N_x( x ) );
    end 
    
    SN_x_ = SN_x_out;
    E_x_ = F_x_;
    ZF_x_ = ( F_x_ ./ L_x_ .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) );
    SP_x_ = ( 1 - gamma ) .* F_x_ ./ ZF_x_;
    P_x_ = P_x_Over_Q_x_ .* Q_x_;
    Z_x_ = ( ( ( K_x_ ./ GYTrend_ ) .^ alpha .* ( A_x_ .* H_x_ ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* SP_x_ ./ P_x_ ) .^ kappa ) .^ ( 1 ./ ( 1 - kappa ) );
    SRK_x_ = ( 1 - kappa ) .* alpha .* SP_x_ .* Z_x_ ./ ( K_x_ ./ GYTrend_ );
    W_x_ = ( 1 - kappa ) .* ( 1 - alpha ) .* SP_x_ .* Z_x_ ./ H_x_;
    C_x_ = thetaC .* F_x_ ./ ( thetaF .* P_x_ );
    I_x_ = K_x_ .* ( 1 - ( 1 - deltaK ) ./ GYTrend_ ) ./ ( 1 - Phi2 ./ 2 .* ( GYTrend_ - 1 ) .^ 2 );
    M_x_ = kappa .* SP_x_ .* Z_x_ ./ P_x_;
    Y_x_ = C_x_ + I_x_ + M_x_;
    %YBar_x_ = Y_x_ .* P_x_ .^ ( ( 1 + lambda ) / lambda ) .* AverageTransportCost_ .^ ( - 1 ./ lambda );
    
    YBar_x_ = zeros( SpatialPoints,1 );
    for x=1:SpatialPoints 
        integral_YBar = zeros( 1,SpatialPoints );
        for tildex=1:9
            integral_YBar( tildex ) = Weight * Y_x_( tildex ) * ( P_x_( tildex ) )^((1+lambda)/lambda) * exp( -tau_ * d( x , tildex ) / lambda );
        end
        YBar_x_( x ) = sum( integral_YBar );
    end
    
    Pi_x_ = lambda ./ ( 1 + lambda ) .* ( 1 + lambda ) .^ ( - 1 ./ lambda ) .* SP_x_ .^ ( - 1 ./ lambda ) .* YBar_x_;    
    J_x_ = ( Z_x_ - ZF_x_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_x_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_x_ ); 

    for x=1:SpatialPoints 
        integral_Px = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_Px( tildex ) = Weight * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
        end
    end
    
    
    U_x_ = ( C_x_ ./ N_x_LAG_ ).^thetaC .* ( E_x_ ./ N_x_LAG_ ).^thetaF ...
        .* ( ( 1 - L_x_ )./N_x_LAG_ ).^thetaL ...
        .* ( ( 1 / ( 1 + nu ) ) * Gamma^( 1 + nu ) - ( 1 / ( 1 + nu ) ) .* ( H_x_ ./ N_x_LAG_ ).^( 1 + nu ) ).^thetaH...
        .* ( .5 * Omega ^ 2 - .5 * ( log( GN_ .* N_x_LAG_ ) ).^2 ).^ thetaN ...
         .* ( 1 - SN_x_ ./ N_x_LAG_ ).^psi1 ...
         .* ( dBar - SD_x ./ SN_x_ ).^ psi2 ...
         .* exp( psi3 .* util_SN );
    
    muN_x_ = beta_ .* ( ( U_x_ .* GUTrend_ ) .^ ( 1 - varsigma ) ) .* ( 1 + ( 1 - varsigma ) .* ( ...
            thetaH .* ( H_x_ ./ N_x_LAG_ ) .^ ( 1 + nu ) ./ ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x_ ./ N_x_LAG_ ) .^ ( 1 + nu ) ) ...
            - thetaN .* log( N_x ) ./ ( 1 ./ 2 .* Omega .^ 2 - 1 ./ 2 .* log( N_x ) .^ 2 ) ...
            + psi1 .* GN_ .* SN_x_ ./ ( N_x - GN_ .* SN_x_ ) ...
            - ( thetaC + thetaF + thetaL + psi3 ) ...
        ) ) ./ ( 1 - beta_ .* GmuNTrend_ .* GN_ );

    Residual = zeros( SpatialPoints , 4+SpatialPoints );
    
    Residual( :, 1 ) = Xi_LEAD_ .* GSRKTrend_ .* SRK_x_ ./ ( 1 - Xi_LEAD_ .* GSRKTrend_ .* ( 1 - deltaK ) ) - Q_x_;
    Residual( :, 2 ) = thetaF .* N_x_ ./ F_x_ .* W_x_ .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x_ ./ N_x_LAG_ ) .^ ( 1 + nu ) ) - thetaH .* ( H_x_ ./ N_x_LAG_ ) .^ nu;
    Residual( :, 3 ) = Pi_x_ + ( 1 - deltaJ ) .* Xi_LEAD_ .* phi_ .* SP_x_ .* GSPTrend_ - phi_ .* SP_x_;
    Residual( :, 4 ) = P_x_ - ( 1 + lambda ) .* sum( integral_Px ).^(-lambda);
    
    
    for x=1:SpatialPoints
        for tildex=1:SpatialPoints
            Residual( x , 4+tildex ) = muN_x_( x ) - muN_x_( tildex ) - ( 1 - varsigma ) .* N_x_LAG_( x ) .* ( U_x_( x ) ).^( 1 - varsigma ) .* ( GN_ .* psi3 .* N_x_LAG_( tildex ) ./ SN_xx_( x , tildex ) - psi1 ./ ( N_x_LAG_( x ) - SN_x_( x ) ) - psi2 .* ( d( ( x ), tildex ) .* ( SN_x_( x ) ) - SD_x( x ) ) ./ ( dBar .* ( SN_x_( x ) ).^2 - SN_x_( x ) .* SD_x( x ) ) );                              
        end
    end

end
