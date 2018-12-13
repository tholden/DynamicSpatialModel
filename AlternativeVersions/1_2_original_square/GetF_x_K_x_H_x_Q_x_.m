function [ F_x_ , K_x_ , H_x_ , Q_x_ ] = GetF_x_K_x_H_x_Q_x_( init_values, A_x_, N_x_LAG_, nu, gamma, Gamma, L_x_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_x_Over_Q_x_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ )
    global F_x_K_x_H_x_Q_x_
    F_x_K_x_H_x_Q_x_ = exp( fsolve( @( log_F_x_K_x_H_x_Q_x_ ) GetResidual( exp( log_F_x_K_x_H_x_Q_x_ ), A_x_, N_x_LAG_, nu, gamma, Gamma, L_x_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_x_Over_Q_x_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ ), ...
        log( init_values ), ...
        optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
        F_x_ = F_x_K_x_H_x_Q_x_( 1 );
        K_x_ = F_x_K_x_H_x_Q_x_( 2 );
        H_x_ = F_x_K_x_H_x_Q_x_( 3 );
        Q_x_ = F_x_K_x_H_x_Q_x_( 4 );
end

function Residual = GetResidual( F_x_K_x_H_x_Q_x_, A_x_, N_x_LAG_, nu, gamma, Gamma, L_x_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_x_Over_Q_x_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ )
    F_x_ = F_x_K_x_H_x_Q_x_( 1 );
    K_x_ = F_x_K_x_H_x_Q_x_( 2 );
    H_x_ = F_x_K_x_H_x_Q_x_( 3 );
    Q_x_ = F_x_K_x_H_x_Q_x_( 4 );
    
    ZF_x_ = ( F_x_ / L_x_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
    SP_x_ = ( 1 - gamma ) * F_x_ / ZF_x_;
    P_x_ = P_x_Over_Q_x_ * Q_x_;
    Z_x_ = ( ( ( K_x_ / GYTrend_ ) ^ alpha * ( A_x_ * H_x_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_x_ / P_x_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
    SRK_x_ = ( 1 - kappa ) * alpha * SP_x_ * Z_x_ / ( K_x_ / GYTrend_ );
    W_x_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_x_ * Z_x_ / H_x_;
    C_x_ = thetaC * F_x_ / ( thetaF * P_x_ );
    I_x_ = K_x_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
    M_x_ = kappa * SP_x_ * Z_x_ / P_x_;
    Y_x_ = C_x_ + I_x_ + M_x_;
    YBar_x_ = Y_x_ * P_x_ ^ ( ( 1 + lambda ) / lambda ) * AverageTransportCost_ ^ ( - 1 / lambda );
    Pi_x_ = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP_x_ ^ ( - 1 / lambda ) * YBar_x_;
    J_x_ = ( ( 1 + lambda ) * SP_x_ * AverageTransportCost_ / P_x_ ) ^ ( 1 / lambda );
    
    Residual = zeros( 4, 1 );
    
    Residual( 1 ) = Xi_LEAD_ * GSRKTrend_ * SRK_x_ / ( 1 - Xi_LEAD_ * GSRKTrend_ * ( 1 - deltaK ) ) - Q_x_;
    Residual( 2 ) = thetaF * N_x_LAG_ / F_x_ * W_x_ * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H_x_ / N_x_LAG_ ) ^ ( 1 + nu ) ) - thetaH * ( H_x_ / N_x_LAG_ ) ^ nu;
    Residual( 3 ) = Pi_x_ + ( 1 - deltaJ ) * Xi_LEAD_ * phi_ * SP_x_ * GSPTrend_ - phi_ * SP_x_;
    Residual( 4 ) = ZF_x_ + phi_ * ( J_x_ - ( 1 - deltaJ ) * J_x_ / GJTrend_ ) + J_x_ * ( ( 1 + lambda ) * SP_x_ ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar_x_ - Z_x_;
end
