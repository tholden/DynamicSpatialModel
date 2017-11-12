function F_1_ = GetF_1_( A_1_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, GQTrend_, deltaK, Phi2, GSPTrend_ )
    global F_1_Z_1_YBar_1_
    F_1_Z_1_YBar_1_ = exp( fsolve( @( log_F_1_Z_1_YBar_1_ ) GetResidual( exp( log_F_1_Z_1_YBar_1_ ), A_1_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, GQTrend_, deltaK, Phi2, GSPTrend_ ), ...
        zeros( 3, 1 ), optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
    F_1_ = F_1_Z_1_YBar_1_( 1 );
end

function Residual = GetResidual( F_1_Z_1_YBar_1_, A_1_, N_1_, N_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, GQTrend_, deltaK, Phi2, GSPTrend_ )
    F_1_ = F_1_Z_1_YBar_1_( 1 );
    Z_1_ = F_1_Z_1_YBar_1_( 2 );
    YBar_1_ = F_1_Z_1_YBar_1_( 3 );
    
    SP_1_ = ( 1 - gamma ) * F_1_ / ( F_1_ / L_1_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );

    P_1_ = ( 1 + lambda ) * ( ( Z_1_ - ( F_1_ / L_1_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) ) ) / ( phi_ * ( 1 - ( 1 - deltaJ ) / GJTrend_ ) + ( ( 1 + lambda ) * SP_1_ ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar_1_ ) ) ^ ( - lambda ) * SP_1_ * AverageTransportCost_;

    Z_1_Alt = A_1_ * ( thetaF * N_1_ / F_1_ * ( 1 - kappa ) * ( 1 - alpha ) * SP_1_ * Z_1_ / ( 1 + nu ) * Gamma ^ ( 1 + nu ) / ( thetaH / N_1_LAG_ ^ nu + thetaF * N_1_ / F_1_ * ( 1 - kappa ) * ( 1 - alpha ) * SP_1_ * Z_1_ / ( 1 + nu ) / N_1_LAG_ ^ ( 1 + nu ) ) ) ^ ( 1 / ( 1 + nu ) ) * ( ( ( 1 - kappa ) * alpha * SP_1_ * GSRKTrend_ * P_1_Over_Q_1_ / P_1_ / ( 1 / Xi_LEAD_ - GQTrend_ * ( 1 - deltaK ) ) ) ^ alpha * ( kappa * SP_1_ / P_1_ ) ^ ( kappa / ( 1 - kappa ) ) ) ^ ( 1 / ( 1 - alpha ) );
    YBar_1_Alt = ( thetaC * F_1_ / ( thetaF * P_1_ ) + ( 1 - kappa ) * alpha * SP_1_ * Z_1_ * GYTrend_ * GSRKTrend_ * P_1_Over_Q_1_ / P_1_ / ( 1 / Xi_LEAD_ - GQTrend_ * ( 1 - deltaK ) ) * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 ) + kappa * SP_1_ * Z_1_ / P_1_ ) * P_1_ ^ ( ( 1 + lambda ) / lambda ) * AverageTransportCost_ ^ ( - 1 / lambda );
    SP_1_Alt = ( ( phi_ - ( 1 - deltaJ ) * Xi_LEAD_ * phi_ * GSPTrend_ ) / ( lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) / YBar_1_ ) ) ^ ( - lambda / ( 1 + lambda ) );
    
    Residual = [ Z_1_ - Z_1_Alt, YBar_1_ - YBar_1_Alt, SP_1_ - SP_1_Alt ];
end
