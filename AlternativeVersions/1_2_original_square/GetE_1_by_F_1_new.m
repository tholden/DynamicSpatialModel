function E_1_by_F_1_ = GetE_1_by_F_1_new( solver, shape , Index, SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma)
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    if SpatialPoints==9
        SpatialPointsToCompute = 3;
    end
    E_by_F_0 = 1 * ones(SpatialPointsToCompute,1); 
    F_0 = 0.0883 * ones(SpatialPointsToCompute,1); 
    K_0 = 0.5574 * ones(SpatialPointsToCompute,1); 
    H_0 = 0.6248 * ones(SpatialPointsToCompute,1); 
    Q_0 = 5.5109 * ones(SpatialPointsToCompute,1); 
    N_0 =  ones(SpatialPointsToCompute,1);
    SN_0 = .2 * ones(SpatialPointsToCompute,1);
    SD_0 = .01 * ones(SpatialPointsToCompute,1);
    exp_util_SN_0 = exp( -2.3026 * ones(SpatialPointsToCompute,1) );
    d = getDistanceMatrix( shape, SpatialPointsPerDimension ); % create distance

    if strcmp(solver,'fmincon')
        % Asymmetric
        %A = zeros(5 , (6+SpatialPoints)*SpatialPoints);
        %A(1,1) = 1; A(2,3) = 1; A(3,7) = 1; A(4,9) = 1; A(5,5) = -1;
        %b = ones(5,1); b(5,1) = -1;
        % Symmetric
        A = []; b = []; lb = [];
        E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x = exp( fmincon( @( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ) GetResidual( exp( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ), solver, shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma), ...
          log( [ E_by_F_0' , F_0' , K_0' , H_0' , Q_0' , N_0' , SN_0' , SD_0' , exp_util_SN_0' ]' ), ...
          A,b,[],[],lb,[],@( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ) constraints( exp( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ), SpatialPointsPerDimension ),...
          ...optimoptions( @fmincon, 'HessianApproximation', {'lbfgs',30}, 'Display', 'iter', 'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12) ) );
          optimoptions( @fmincon,'Display', 'iter', 'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false) ) );
    elseif strcmp(solver,'fsolve')
        E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x = exp( fsolve( @( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ) GetResidual( exp( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ), solver, shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma), ...
        log( [ E_by_F_0' , F_0' , K_0' , H_0' , Q_0' , N_0' , SN_0' , SD_0' , exp_util_SN_0' ]' ), ...
        optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
    elseif strcmp(solver,'fmincon_global')
        %'HessianApproximation', {'lbfgs',30}, 
        opts = optimoptions( @fmincon, 'Display', 'iter', 'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false);
        problem = createOptimProblem('fmincon',...
        'objective', @( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ) GetResidual( exp( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ),  solver, shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma),...
        'x0',log( [ E_by_F_0' , F_0' , K_0' , H_0' , Q_0' , N_0' , SN_0' , SD_0' , exp_util_SN_0' ]' ),'options',opts );
        problem.nonlcon = @( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ) constraints( exp( log_E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x ), SpatialPointsPerDimension);
        gs = GlobalSearch(MultiStart('FunctionTolerance',1e-8,'UseParallel',true));
        E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x = exp( run(gs,problem) );
    end

global E_by_F_x F_x K_x H_x Q_x N_x_ SN_x_ SD_x_ util_SN_x_
    E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x = reshape(E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x,[SpatialPointsToCompute,9]);
    E_by_F_x    = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 1 );
    F_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 2 );
    K_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 3 );
    H_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 4 );
    Q_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 5 );
    N_x_        = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 6 );
    SN_x_       = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 7 );
    SD_x_       = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 8 );
    util_SN_x_  = log( E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 9 ) );
    
    E_1_by_F_1_ = E_by_F_x( Index );
end

function Residual = GetResidual( E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x,  solver, shape, SpatialPointsPerDimension, d, GN_, nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_, varsigma )
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    if SpatialPoints==9
        SpatialPointsToCompute = 3;
    end
    Weight = getWeightMatrix( shape , SpatialPointsPerDimension );
    Weight = reshape(Weight,[SpatialPoints,1]);
    E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x = reshape(E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x,[SpatialPointsToCompute,9]);
    E_by_F_x    = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 1 );
    F_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 2 );
    K_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 3 );
    H_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 4 );
    Q_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 5 );
    N_x_        = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 6 );
    SN_x_       = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 7 );
    SD_x_       = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 8 );
    util_SN_x_  = log( E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 9 ) );
    
    if SpatialPoints==9
        E_by_F_x = [ E_by_F_x( 1 ) E_by_F_x( 2 ) E_by_F_x( 1 ) E_by_F_x( 2 ) E_by_F_x( 3 ) E_by_F_x( 2 ) E_by_F_x( 1 ) E_by_F_x( 2 ) E_by_F_x( 1 )   ]';
        F_x = [ F_x( 1 ) F_x( 2 ) F_x( 1 ) F_x( 2 ) F_x( 3 ) F_x( 2 ) F_x( 1 ) F_x( 2 ) F_x( 1 )   ]';
        K_x = [ K_x( 1 ) K_x( 2 ) K_x( 1 ) K_x( 2 ) K_x( 3 ) K_x( 2 ) K_x( 1 ) K_x( 2 ) K_x( 1 )   ]';
        H_x = [ H_x( 1 ) H_x( 2 ) H_x( 1 ) H_x( 2 ) H_x( 3 ) H_x( 2 ) H_x( 1 ) H_x( 2 ) H_x( 1 )   ]';
        Q_x = [ Q_x( 1 ) Q_x( 2 ) Q_x( 1 ) Q_x( 2 ) Q_x( 3 ) Q_x( 2 ) Q_x( 1 ) Q_x( 2 ) Q_x( 1 )   ]';
        N_x_ = [ N_x_( 1 ) N_x_( 2 ) N_x_( 1 ) N_x_( 2 ) N_x_( 3 ) N_x_( 2 ) N_x_( 1 ) N_x_( 2 ) N_x_( 1 )   ]';
        SN_x_ = [ SN_x_( 1 ) SN_x_( 2 ) SN_x_( 1 ) SN_x_( 2 ) SN_x_( 3 ) SN_x_( 2 ) SN_x_( 1 ) SN_x_( 2 ) SN_x_( 1 )   ]';
        SD_x_ = [ SD_x_( 1 ) SD_x_( 2 ) SD_x_( 1 ) SD_x_( 2 ) SD_x_( 3 ) SD_x_( 2 ) SD_x_( 1 ) SD_x_( 2 ) SD_x_( 1 )   ]';
        util_SN_x_ = [ util_SN_x_( 1 ) util_SN_x_( 2 ) util_SN_x_( 1 ) util_SN_x_( 2 ) util_SN_x_( 3 ) util_SN_x_( 2 ) util_SN_x_( 1 ) util_SN_x_( 2 ) util_SN_x_( 1 )   ]';
    end
    
    A_x = 1;
    Px_Over_Q_ = ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 - Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ) + Xi_LEAD_ * GSRKTrend_ * Phi2 * ( GYTrend_ - 1 ) * GYTrend_ ^ 2;
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

    integral_Px = zeros( SpatialPoints,1 );
    for x=1:SpatialPoints 
        integral_Px_term = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_Px_term( tildex ) = Weight( tildex ) * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
        end
        integral_Px( x ) = sum( integral_Px_term );
    end
                
    N_x_LAG_ = N_x_ ./ GN_;

    U_x_ = ( C_x_ ./ N_x_LAG_ ).^thetaC .* ( E_x_ ./ N_x_LAG_ ).^thetaF ...
            .* ( ( 1 - L_x_ )./N_x_LAG_ ).^thetaL ...
            .* ( ( 1 / ( 1 + nu ) ) * Gamma^( 1 + nu ) - ( 1 / ( 1 + nu ) ) .* ( H_x ./ N_x_LAG_ ).^( 1 + nu ) ).^thetaH...
            .* ( .5 * Omega ^ 2 - .5 * log( GN_ .* N_x_LAG_ ).^2 ).^ thetaN ...
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
        Residual = zeros( 3 , 9 );
    end 
    
    
    Residual( 1 , 1 )   = sum( E_x_ ) - sum( F_x );
    Residual( 2 , 1 )   = E_x_( 1 ) .* N_x_LAG_( 2 ) .* ( U_x_( 2 ) ).^( 1 - varsigma ) - E_x_( 2 ) .* N_x_LAG_( 1 ) .* ( U_x_( 1 ) ).^( 1 - varsigma );
    Residual( 3 , 1 )   = E_x_( 1 ) .* N_x_LAG_( 5 ) .* ( U_x_( 5 ) ).^( 1 - varsigma ) - E_x_( 5 ) .* N_x_LAG_( 1 ) .* ( U_x_( 1 ) ).^( 1 - varsigma );
    Residual( 1 , 2 ) = Xi_LEAD_ .* GSRKTrend_ .* ( SRK_x_( 1 ) +  ( 1 - deltaK ) .* Q_x( 1 ) ) - Q_x( 1 );
    Residual( 2 , 2 ) = Xi_LEAD_ .* GSRKTrend_ .* ( SRK_x_( 2 ) +  ( 1 - deltaK ) .* Q_x( 2 ) ) - Q_x( 2 );
    Residual( 3 , 2 ) = Xi_LEAD_ .* GSRKTrend_ .* ( SRK_x_( 5 ) +  ( 1 - deltaK ) .* Q_x( 5 ) ) - Q_x( 5 );
    Residual( 1 , 3 ) = thetaF .* N_x_LAG_( 1 ) ./ E_x_( 1 ) .* W_x_( 1 ) .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x( 1 ) ./ N_x_LAG_( 1 ) ).^ ( 1 + nu ) ) - thetaH .* ( H_x( 1 ) ./ N_x_LAG_( 1 ) ).^ nu;
    Residual( 2 , 3 ) = thetaF .* N_x_LAG_( 2 ) ./ E_x_( 1 ) .* W_x_( 2 ) .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x( 2 ) ./ N_x_LAG_( 2 ) ).^ ( 1 + nu ) ) - thetaH .* ( H_x( 2 ) ./ N_x_LAG_( 2 ) ).^ nu;
    Residual( 3 , 3 ) = thetaF .* N_x_LAG_( 5 ) ./ E_x_( 1 ) .* W_x_( 5 ) .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x( 5 ) ./ N_x_LAG_( 5 ) ).^ ( 1 + nu ) ) - thetaH .* ( H_x( 5 ) ./ N_x_LAG_( 5 ) ).^ nu;
    Residual( 1 , 4 ) = lambda ./ ( 1 + lambda ) .* ( 1 + lambda ) .^ ( - 1 ./ lambda ) .* SP_x_( 1 ).^ ( - 1 ./ lambda ) .* YBar_x_( 1 ) + ( 1 - deltaJ ).* Xi_LEAD_.* phi_ .* SP_x_( 1 ) .* GSPTrend_ - phi_.* SP_x_( 1 );
    Residual( 2 , 4 ) = lambda ./ ( 1 + lambda ) .* ( 1 + lambda ) .^ ( - 1 ./ lambda ) .* SP_x_( 2 ).^ ( - 1 ./ lambda ) .* YBar_x_( 2 ) + ( 1 - deltaJ ).* Xi_LEAD_.* phi_ .* SP_x_( 2 ) .* GSPTrend_ - phi_.* SP_x_( 2 );
    Residual( 3 , 4 ) = lambda ./ ( 1 + lambda ) .* ( 1 + lambda ) .^ ( - 1 ./ lambda ) .* SP_x_( 5 ).^ ( - 1 ./ lambda ) .* YBar_x_( 5 ) + ( 1 - deltaJ ).* Xi_LEAD_.* phi_ .* SP_x_( 5 ) .* GSPTrend_ - phi_.* SP_x_( 5 );
    Residual( 1 , 5 ) = P_x_( 1 ) - ( 1 + lambda ) .* ( integral_Px( 1 ) ).^(-lambda);
    Residual( 2 , 5 ) = P_x_( 2 ) - ( 1 + lambda ) .* ( integral_Px( 2 ) ).^(-lambda);
    Residual( 3 , 5 ) = P_x_( 5 ) - ( 1 + lambda ) .* ( integral_Px( 5 ) ).^(-lambda);
    Residual( 1 , 6 ) = 1 - sum( integral_Nx );
    Residual( 2 , 6 ) = SN_x_out( 1 ) - SN_x_in( 1 );
    Residual( 3 , 6)  = SN_x_out( 2 ) - SN_x_in( 2 );
    Residual( 1 , 7 ) = SN_x_out( 1 ) - SN_x_( 1 );
    Residual( 2 , 7 ) = SN_x_out( 2 ) - SN_x_( 2 );
    Residual( 3 , 7 ) = SN_x_out( 5 ) - SN_x_( 5 );
    Residual( 1 , 8 ) = SD_x_( 1 ) - SD( 1 );
    Residual( 2 , 8 ) = SD_x_( 2 ) - SD( 2 );
    Residual( 3 , 8 ) = SD_x_( 5 ) - SD( 5 );
    Residual( 1 , 9 ) = util_SN_x_( 1 ) - util_SN( 1 );
    Residual( 2 , 9 ) = util_SN_x_( 2 ) - util_SN( 2 );
    Residual( 3 , 9 ) = util_SN_x_( 5 ) - util_SN( 5 );
    
    
    if strcmp(solver,'fmincon') || strcmp(solver,'fmincon_global')
        Residual = reshape(Residual,[1,numel(Residual)]);
        Residual = sum( Residual.^2 );
    end
end


function [c,ceq] = constraints( E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x, SpatialPointsPerDimension )

    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    if SpatialPoints==9
        SpatialPointsToCompute = 3;
    end
    E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x = reshape(E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x,[SpatialPointsToCompute,9]);
    E_by_F_x    = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 1 );
    F_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 2 );
    K_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 3 );
    H_x         = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 4 );
    N_x_        = E_x_by_F_x_F_x_K_x_H_x_Q_x_N_x_SN_x_SD_x_util_SN_x( : , 6 );
  
  % Population higher closer to the centre / population equal everywhere
%   %c   =   []; % e.g. where x1 - 1 >= x2, c(1) = x2 + 1 - x1

% for SpatialPoints==9, pt 1 = 1, pt 2 = 2, pt 3 = 5
  ceq = [];
  c( 1 )    =   N_x_( 2 ) - N_x_( 3 ); 
  c( 2 )    =   N_x_( 1 ) - N_x_( 2 );
  c( 3 )    =   K_x( 2 ) - K_x( 3 ); 
  c( 4 )    =   K_x( 1 ) - K_x( 2 );
  c( 5 )    =   F_x( 3 ) - F_x( 2 ); 
  c( 6 )    =   F_x( 2 ) - F_x( 1 );
  c( 7 )    =   H_x( 2 ) - H_x( 3 ); 
  c( 8 )    =   H_x( 1 ) - H_x( 2 );
  c( 9 )    =   E_by_F_x( 2 ) - E_by_F_x( 3 );
  c( 10 )   =   E_by_F_x( 1 ) - E_by_F_x( 2 );
  
  % for torus 
%   ceq = c;
%   c = [];
  
end        
        
        