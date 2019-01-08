function Residual = GetResidual_fsolve( E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_, shape, SpatialPointsPerDimension, d, GN_, nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_, varsigma )
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    Weight = getWeightMatrix( shape , SpatialPointsPerDimension );
    Weight = reshape(Weight,[SpatialPoints,1]);
    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints );
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints );
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints );
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints );
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints );
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints );
    SNxx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end );
    SNxx_       = reshape(SNxx_,[SpatialPoints,SpatialPoints]);
    
    load('Omega_vec.mat','Omega_vec');
    Omega = reshape(Omega_vec,[SpatialPoints,1]);
    
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
    C_x_ = thetaC .* E_x_ ./ ( thetaF .* P_x_ );
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
        .* ( .5 * Omega .^ 2 - .5 * ( log( GN_ .* N_x_LAG_ ) ).^2 ).^ thetaN ...
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
    Residual( 2:end , 1 ) = E_x_( 1 ) .* N_x_LAG_( 2:end ) .* ( U_x_( 2:end ) ).^( 1 - varsigma ) - E_x_( 2:end ) .* N_x_LAG_( 1 ) .* ( U_x_( 1 ) ).^( 1 - varsigma );
    Residual( 1 , 2 ) = sum( integral_N ) - 1;
    Residual( 2:end , 2 ) = SN_x_out( 2:end ) - SN_x_in( 2:end);
    Residual( : , 3 ) = Xi_LEAD_ .* GSRKTrend_ .* ( SRK_x_ + ( 1 - deltaK ) .* Q_x ) - Q_x;
    Residual( : , 4 ) = thetaF .* N_x_LAG_ ./ E_x_ .* W_x_ .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x ./ N_x_LAG_ ).^ ( 1 + nu ) ) - thetaH .* ( H_x ./ N_x_LAG_ ).^ nu;
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
