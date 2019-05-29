function F_1_ = GetF_1_( A_1_, GWTrend_, GFTrend_, GN_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, PhiD, phiD, deltaD, thetaN )
    global F_1_K_1_H_1_Q_1_
    F_1_K_1_H_1_Q_1_ = exp( fsolve( @( log_F_1_K_1_H_1_Q_1_ ) GetResidual( exp( log_F_1_K_1_H_1_Q_1_ ), A_1_, GWTrend_, GFTrend_, GN_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, PhiD, phiD, deltaD, thetaN ), ...
        [ -3.17866 -2.30363 -1.5714 2.67453 ], ...log( [ 0.0756 0.3423 0.6247 7.6171 ] ), ...
        optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
    F_1_ = F_1_K_1_H_1_Q_1_( 1 );
    disp( F_1_K_1_H_1_Q_1_ );
end

function Residual = GetResidual( F_1_K_1_H_1_Q_1_, A_1_, GWTrend_, GFTrend_, GN_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, PhiD, phiD, deltaD, thetaN )
    F_1_ = F_1_K_1_H_1_Q_1_( 1 );
    K_1_ = F_1_K_1_H_1_Q_1_( 2 );
    H_1_ = F_1_K_1_H_1_Q_1_( 3 );
    Q_1_ = F_1_K_1_H_1_Q_1_( 4 );
     
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
    QD_1_ = ( 1 - PhiD * ( 1 - GN_ )^2 / 2 + PhiD * ( 1 - GN_ ) * GN_ - Xi_LEAD_ * GWTrend_ * PhiD * ( 1 - GN_ ) * GN_^2 )^(-1);
    Omega_1_  = N_1_  + phiD * Xi_LEAD_ * GFTrend_ * ( F_1_ / W_1_ ) * ( thetaN / thetaF ) / ( QD_1_ - Xi_LEAD_ * GWTrend_ * ( 1 - deltaD ) * QD_1_ ); 
    HD_1_ = Omega_1_ * ( 1 - (1-deltaD)/GN_ ) / ( phiD * ( 1 - (PhiD/2) * ( 1 - GN_ )^2 ) ); 
    HT_1_ = H_1_ + HD_1_;
    
    Residual = zeros( 4, 1 );
    
    Residual( 1 ) = Xi_LEAD_ * GSRKTrend_ * SRK_1_ / ( 1 - Xi_LEAD_ * GSRKTrend_ * ( 1 - deltaK ) ) - Q_1_;
    Residual( 2 ) = thetaF * N_1_ / F_1_ * W_1_ * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( HT_1_ / N_1_LAG_ ) ^ ( 1 + nu ) ) - thetaH * ( HT_1_ / N_1_LAG_ ) ^ nu;
    Residual( 3 ) = Pi_1_ + ( 1 - deltaJ ) * Xi_LEAD_ * phi_ * SP_1_ * GSPTrend_ - phi_ * SP_1_;
    Residual( 4 ) = ZF_1_ + phi_ * ( J_1_ - ( 1 - deltaJ ) * J_1_ / GJTrend_ ) + J_1_ * ( ( 1 + lambda ) * SP_1_ ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar_1_ - Z_1_;
end
