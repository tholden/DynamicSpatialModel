function N_x_SN_x_SD_x_util_SN_x = GetN_x_SN_x_SD_x_util_SN_x( N_0_SN_0_SD_0_util_SN_0 , SpatialPoints,  GN_ , GUTrend_, GmuNTrend_, thetaC, thetaF, thetaL, thetaN, thetaH, nu, Gamma, Omega, dBar, d, Weight, psi1, psi2, psi3, beta_, varsigma, C_x_, E_x_, L_x_, H_x ) 
    N_x_SN_x_SD_x_util_SN_x = exp( fsolve( @( log_N_x_SN_x_SD_x_util_SN_x ) GetResidual( exp( log_N_x_SN_x_SD_x_util_SN_x ) , SpatialPoints,  GN_ , GUTrend_, GmuNTrend_, thetaC, thetaF, thetaL, thetaN, thetaH, nu, Gamma, Omega, dBar, d, Weight, psi1, psi2, psi3, beta_, varsigma, C_x_, E_x_, L_x_, H_x), ...
        log( N_0_SN_0_SD_0_util_SN_0 ), ...
        optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
    disp( N_x_SN_x_SD_x_util_SN_x );
end

function Residual = GetResidual( N_x_SN_x_SD_x_util_SN_x , SpatialPoints,  GN_ , GUTrend_, GmuNTrend_, thetaC, thetaF, thetaL, thetaN, thetaH, nu, Gamma, Omega, dBar, d, Weight, psi1, psi2, psi3, beta_, varsigma, C_x_, E_x_, L_x_, H_x )

    N_x_        = N_x_SN_x_SD_x_util_SN_x( : , 1 );
    SN_x_       = N_x_SN_x_SD_x_util_SN_x( : , 2 );
    SD_x_       = N_x_SN_x_SD_x_util_SN_x( : , 3 );
    util_SN_x_  = N_x_SN_x_SD_x_util_SN_x( : , 4 );

    if SpatialPoints==9
        N_x_ = [ N_x_( 1 ) N_x_( 2 ) N_x_( 1 ) N_x_( 2 ) N_x_( 3 ) N_x_( 2 ) N_x_( 1 ) N_x_( 2 ) N_x_( 1 )   ]';
        SN_x_ = [ SN_x_( 1 ) SN_x_( 2 ) SN_x_( 1 ) SN_x_( 2 ) SN_x_( 3 ) SN_x_( 2 ) SN_x_( 1 ) SN_x_( 2 ) SN_x_( 1 )   ]';
        SD_x_ = [ SD_x_( 1 ) SD_x_( 2 ) SD_x_( 1 ) SD_x_( 2 ) SD_x_( 3 ) SD_x_( 2 ) SD_x_( 1 ) SD_x_( 2 ) SD_x_( 1 )   ]';
        util_SN_x_ = [ util_SN_x_( 1 ) util_SN_x_( 2 ) util_SN_x_( 1 ) util_SN_x_( 2 ) util_SN_x_( 3 ) util_SN_x_( 2 ) util_SN_x_( 1 ) util_SN_x_( 2 ) util_SN_x_( 1 )   ]';
    end
    
    N_x_LAG_ = N_x_ ./ GN_;

    U_x_ = ( C_x_ ./ N_x_LAG_ ).^thetaC .* ( E_x_ ./ N_x_LAG_ ).^thetaF ...
            .* ( ( 1 - L_x_ )./N_x_LAG_ ).^thetaL ...
            .* ( ( 1 / ( 1 + nu ) ) * Gamma^( 1 + nu ) - ( 1 / ( 1 + nu ) ) .* ( H_x ./ N_x_LAG_ ).^( 1 + nu ) ).^thetaH...
            .* ( .5 * Omega ^ 2 - .5 * ( log( GN_ .* N_x_LAG_ ) ).^2 ).^ thetaN ...
             .* ( 1 - SN_x_ ./ N_x_LAG_ ).^psi1 ...
             .* ( dBar - SD_x_ ./ SN_x_ ).^ psi2 ...
             .* exp( psi3 .* util_SN_x_ );

    muN_x_ = beta_ .* ( ( U_x_ .* GUTrend_ ) .^ ( 1 - varsigma ) ) .* ( 1 + ( 1 - varsigma ) .* ( ...
                thetaH .* ( H_x ./ N_x_LAG_ ) .^ ( 1 + nu ) ./ ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x ./ N_x_LAG_ ) .^ ( 1 + nu ) ) ...
                - thetaN .* log( N_x_ ) ./ ( 1 ./ 2 .* Omega .^ 2 - 1 ./ 2 .* log( N_x_ ) .^ 2 ) ...
                + psi1 .* GN_ .* SN_x_ ./ ( N_x_ - GN_ .* SN_x_ ) ...
                - ( thetaC + thetaF + thetaL + psi3 ) ...
            ) ) ./ ( 1 - beta_ .* GmuNTrend_ .* GN_ );

    SN_xx_ = zeros( SpatialPoints , SpatialPoints );
    for xtilde=1:SpatialPoints
        SN_xx_( : , xtilde ) = psi3 .* GN_  .* N_x_LAG_ ./ ( ( muN_x_ - muN_x_( xtilde ) ) ./ ( ( 1 - varsigma ) .* N_x_LAG_ .* U_x_.^ ( 1 - varsigma ) ) + psi1 ./ ( N_x_LAG_ - SN_x_ ) + psi2 .* ( d( : , xtilde ) .* SN_x_ - SD_x_ ) ./ ( dBar .* SN_x_ .* SN_x_ - SN_x_ .* SD_x_ ) );
    end

    SN_x_out = zeros( SpatialPoints,1 );
    SN_x_in = zeros( SpatialPoints,1 );
    integral_Nx = zeros( SpatialPoints,1 );
    SD = zeros( SpatialPoints,1 );
    util_SN = zeros( SpatialPoints,1 );
    for x=1:SpatialPoints 
        integral_SN_out = zeros( SpatialPoints,1 );
        integral_SN_in = zeros( SpatialPoints,1 );
        integral_SD = zeros( SpatialPoints,1 );
        integral_util_SN = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_SN_out( tildex ) = Weight( tildex ) * ( SN_xx_( x , tildex ) );
            integral_SN_in( tildex ) = Weight( tildex ) * ( SN_xx_( tildex , x ) );
            integral_SD( tildex ) = Weight( tildex ) * d( x , tildex ) * SN_xx_( x , tildex  );
            integral_util_SN( tildex ) = Weight( tildex ) * GN_ * N_x_LAG_( tildex ) * log( SN_xx_( x , tildex ) / N_x_LAG_( x ) );
        end
        SN_x_out( x ) = sum( integral_SN_out );
        SN_x_in( x ) = sum( integral_SN_in );
        integral_Nx( x ) = Weight( x ) * ( N_x_( x ) );
        SD( x ) = sum( integral_SD );
        util_SN( x ) = sum( integral_util_SN );
    end
    
    if SpatialPoints==9
        Residual = zeros( 3 , 4 );
    end
    
    Residual( 1 , 1 ) = 1 - sum( integral_Nx );
    Residual( 2 , 1 ) = SN_x_out( 2 ) - SN_x_in( 2 );
    Residual( 3 , 1 ) = SN_x_out( 3 ) - SN_x_in( 3 );
    Residual( 1 , 2 ) = SN_x_out( 1 ) - SN_x_( 1 );
    Residual( 2 , 2 ) = SN_x_out( 2 ) - SN_x_( 2 );
    Residual( 3 , 2 ) = SN_x_out( 5 ) - SN_x_( 5 );
    Residual( 1 , 3 ) = SD_x_( 1 ) - SD( 1 );
    Residual( 2 , 3 ) = SD_x_( 2 ) - SD( 2 );
    Residual( 3 , 3 ) = SD_x_( 5 ) - SD( 5 );
    Residual( 1 , 4 ) = util_SN_x_( 1 ) - util_SN( 1 );
    Residual( 2 , 4 ) = util_SN_x_( 2 ) - util_SN( 2 );
    Residual( 3 , 4 ) = util_SN_x_( 5 ) - util_SN( 5 );

end