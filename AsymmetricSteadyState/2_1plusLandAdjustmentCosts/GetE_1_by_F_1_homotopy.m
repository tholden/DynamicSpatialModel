function E_1_by_F_1_ = GetE_1_by_F_1_homotopy( shape , SpatialPointsPerDimension,  GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma)
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ first_run
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    d = getDistanceMatrix( shape, SpatialPointsPerDimension ); % create distance

    if first_run==1
        % First iteration (symmetric equilibrium)
        E_by_F_0 = 1 * ones(SpatialPoints,1);
        N_0 = 1 * ones(SpatialPoints,1);
        F_0 = 0.0787 * ones(SpatialPoints,1);
        K_0 = 0.3911 * ones(SpatialPoints,1);
        H_0 = 0.6252 * ones(SpatialPoints,1);
        Q_0 = 7.0009 * ones(SpatialPoints,1);  
        SN_0_ = psi3 / ( psi1 + psi3 ) * N_0;
        N_LAG = 1/GN_;
        N_0_LAG = N_0./GN_;
        SD_0 = zeros(SpatialPoints,1);
        for i=1:SpatialPoints
            SD_0( i ) = SolveInitialSD( psi1, psi2, psi3, N_0_LAG( i ), N_LAG,  SN_0_( i ), dBar, d( i, : ) ); 
        end
        for i=1:SpatialPoints
            for j=1:SpatialPoints
                SN_xx_0( i , j) = ( N_0_LAG( j ) / N_LAG ) * psi3 / ( psi1 / ( N_0_LAG( i ) - SN_0_( i ) ) + psi2 * ( d( i , j ) * SN_0_( i ) - SD_0( i ) ) / ( dBar * SN_0_( i ) * SN_0_( i ) - SN_0_( i ) * SD_0( i ) ) );
            end
        end
        SN_xx_0 = reshape(SN_xx_0,[SpatialPoints*SpatialPoints , 1]);
        
        disp('Solving symmetric equilibrium case...')
        homotopy = 0;
        steadystate_lsqnonlin.options =  optimoptions( 'lsqnonlin' , 'Display', 'iter', 'FunctionTolerance', 1e-8, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);
        steadystate_lsqnonlin.solver = 'lsqnonlin';
        steadystate_lsqnonlin.objective = @( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ) GetResidual( exp( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ), shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma, homotopy);
        steadystate_lsqnonlin.x0 = log( [ E_by_F_0' , N_0' , F_0' , K_0' , H_0' , Q_0' , SN_xx_0' ]' );
        E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ = exp( lsqnonlin(steadystate_lsqnonlin) );

        crit=0;
        step=0.1;
        while crit<1
            E_by_F_0    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints );
            N_0         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints );
            F_0         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints );
            K_0         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints );
            H_0         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints );
            Q_0         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints );
            SN_xx_0       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end );
            %lh=homotopy;homotopy=min(1,homotopy+(1.1-homotopy)/2);
            lh=homotopy;homotopy=min(1,lh+step);
            count=0;error_count=0;
            while count<=error_count   
                try  
                    disp(['Solving asymmetric equilibrium via homotopy. Scaling factor: ',num2str(homotopy)])
%                     steadystate_lsqnonlin.options =  optimoptions( 'lsqnonlin' , 'Display', 'iter','MaxIterations', 25);%,'algorithm','levenberg-marquardt');
%                     steadystate_lsqnonlin.solver = 'lsqnonlin';
%                     steadystate_lsqnonlin.objective = @( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ) GetResidual( exp( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ), shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma, homotopy);
%                     steadystate_lsqnonlin.x0 = log( [ E_by_F_0' , N_0' , F_0' , K_0' , H_0' , Q_0' , SN_xx_0' ]' );
%                     [ x_temp,resnorm,res,exitflag,output ] = lsqnonlin(steadystate_lsqnonlin);
%                     if resnorm>1e-8
%                         error('Equilbrium not solved')
%                     end
%                     if exitflag==0
%                         error('Max iterations reached')
%                     end
                    
                    steadystate_fsolve.options =  optimoptions( 'fsolve' , 'Display', 'iter','MaxIterations', 40);
                    steadystate_fsolve.solver = 'fsolve';
                    steadystate_fsolve.objective = @( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ) GetResidual( exp( log_E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ ), shape, SpatialPointsPerDimension, d, GN_ , nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_ , varsigma, homotopy);
                    steadystate_fsolve.x0 = log( [ E_by_F_0' , N_0' , F_0' , K_0' , H_0' , Q_0' , SN_xx_0' ]' );
                    [ x_temp,fval,exitflag,output ] = fsolve(steadystate_fsolve);
                    if ~exitflag==1
                        error('Equilbrium not solved')
                    end
                    
                    disp(['Solved equilibrium at homotopy = ',num2str(homotopy)])
                    E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_ = exp( x_temp );
                    step=step*2;
                    crit = homotopy;
                catch
                    step = step/2;
                    disp(['Solver failed. Reducing scaling factor from ',num2str(homotopy),' to ',num2str(lh+step),' and halving step to ',num2str(step)])
                    homotopy=min(1,lh+step); 
                    error_count=error_count+1;
                end
                count = count+1;
            end
        end
    end
     
     figure;
        subplot(3,2,1);surf(reshape(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints ),[SpatialPointsPerDimension,SpatialPointsPerDimension])); title('E_by_F');
        subplot(3,2,2);surf(reshape(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints ),[SpatialPointsPerDimension,SpatialPointsPerDimension])); title('N');
        subplot(3,2,3);surf(reshape(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints ),[SpatialPointsPerDimension,SpatialPointsPerDimension])); title('F');
        subplot(3,2,4);surf(reshape(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints ),[SpatialPointsPerDimension,SpatialPointsPerDimension])); title('K');
        subplot(3,2,5);surf(reshape(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints ),[SpatialPointsPerDimension,SpatialPointsPerDimension])); title('H');
        subplot(3,2,6);surf(reshape(E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints ),[SpatialPointsPerDimension,SpatialPointsPerDimension])); title('Q');
                    
    
    first_run=0;
    
    global E_by_F_x N_x F_x K_x H_x Q_x SNxx_
    E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints );
    N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints );
    F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints );
    K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints );
    H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints );
    Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints );
    SNxx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end );
    SNxx_ = reshape(SNxx_,[SpatialPoints,SpatialPoints]);
    
    E_1_by_F_1_ = E_by_F_x( 1 );
    
end

function Residual = GetResidual( E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_, shape, SpatialPointsPerDimension, d, GN_, nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_, varsigma, homotopy )
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

    integral_Px = zeros( SpatialPoints,1 );
    for x=1:SpatialPoints 
        integral_Px_temp = zeros( SpatialPoints,1 );
        for tildex=1:SpatialPoints
            integral_Px_temp( tildex ) = Weight( tildex ) * J_x_( tildex ) * ( SP_x_( tildex ) * exp( tau_ * d( x , tildex ) ) )^(-1/lambda);
        end
        integral_Px( x,1 ) = sum( integral_Px_temp );
    end
    
    Utilde_x_ = ones(SpatialPoints,1);
    for x=1:SpatialPoints
        Utilde_x_( x ) = Get_Utilde(  SpatialPointsPerDimension, x, homotopy ); 
    end
        
    U_x_ = ( C_x_ ./ N_x_LAG_ ).^thetaC .* ( E_x_ ./ N_x_LAG_ ).^thetaF ...
        .* ( ( 1 - L_x_ )./N_x_LAG_ ).^thetaL ...
        .* ( ( 1 / ( 1 + nu ) ) * Gamma^( 1 + nu ) - ( 1 / ( 1 + nu ) ) .* ( H_x ./ N_x_LAG_ ).^( 1 + nu ) ).^thetaH...
        .* ( .5 * Omega .^ 2 - .5 * ( log( GN_ .* N_x_LAG_ ) ).^2 ).^ thetaN ...
         .* ( 1 - SN_x_ ./ N_x_LAG_ ).^psi1 ...
         .* ( dBar - SD_x ./ SN_x_ ).^ psi2 ...
         .* exp( psi3 .* util_SN ) .* Utilde_x_;
    
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
    Residual( : , 6 ) = P_x_ - ( 1 + lambda ) .* integral_Px.^(-lambda);
        
    Residual = reshape(Residual,[1,numel(Residual)]);

    Residual_SN = zeros( SpatialPoints , SpatialPoints); 
    for x=1:SpatialPoints
        for tildex=1:SpatialPoints
            Residual_SN( x , tildex ) = muN_x_( x ) - muN_x_( tildex ) - ( 1 - varsigma ) .* N_x_LAG_( x ) .* ( U_x_( x ) ).^( 1 - varsigma ) .* ( GN_ .* psi3 .* N_x_LAG_( tildex ) ./ SNxx_( x , tildex ) - psi1 ./ ( N_x_LAG_( x ) - SN_x_( x ) ) - psi2 .* ( d( ( x ), tildex ) .* ( SN_x_( x ) ) - SD_x( x ) ) ./ ( dBar .* ( SN_x_( x ) ).^2 - SN_x_( x ) .* SD_x( x ) ) );                              
        end
    end
    
     Residual = [ Residual reshape(Residual_SN,[1,numel(Residual_SN)]) ];
end




        
        