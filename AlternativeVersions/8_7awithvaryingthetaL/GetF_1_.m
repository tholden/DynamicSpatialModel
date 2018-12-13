function F_1_ = GetF_1_( A_1_, Nr_1_, Nr_1_LAG_, Ns_1_, Ns_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GZTrend_, xi, GPTrend_, GYBarTrend_, GOmega1Trend_, Pi_)
    global F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_
    F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_ = exp( fsolve( @( log_F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_ ) GetResidual( exp( log_F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_ ), A_1_, Nr_1_, Nr_1_LAG_, Ns_1_, Ns_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GZTrend_, xi, GPTrend_, GYBarTrend_, GOmega1Trend_, Pi_), ...
        log( [ 0.0882 0.5575 0.6255 5.5178 1 1 ] ), ...
        optimoptions( @fsolve, 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', false ) ) );
    F_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 1 );
    disp( F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_ );
end

function Residual = GetResidual( F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_, A_1_, Nr_1_, Nr_1_LAG_, Ns_1_, Ns_1_LAG_, nu, gamma, Gamma, L_1_, lambda, phi_, deltaJ, GJTrend_, AverageTransportCost_, thetaL, thetaC, thetaF, thetaH, kappa, alpha, GYTrend_, GSRKTrend_, P_1_Over_Q_1_, Xi_LEAD_, deltaK, Phi2, GSPTrend_ , GZTrend_, xi, GPTrend_ , GYBarTrend_, GOmega1Trend_, Pi_)
    F_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 1 );
    K_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 2 );
    H_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 3 );
    Q_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 4 );
    Er_by_Es_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 5 );
    Hr_by_Hs_1_ = F_1_K_1_H_1_Q_1_Er_by_Es_1_Hr_by_Hs_1_( 6 );

    ZF_1_ = ( F_1_ / L_1_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) );
    SP_1_ = ( 1 - gamma ) * F_1_ / ZF_1_;
    P_1_ = P_1_Over_Q_1_ * Q_1_;
    Z_1_ = ( ( ( K_1_ / GYTrend_ ) ^ alpha * ( A_1_ * H_1_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_1_ / P_1_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
    SRK_1_ = ( 1 - kappa ) * alpha * SP_1_ * Z_1_ / ( K_1_ / GYTrend_ );
    W_1_ = ( 1 - kappa ) * ( 1 - alpha ) * SP_1_ * Z_1_ / H_1_;
    Es_1_ = F_1_ / ( 1 + Er_by_Es_1_ );
    Er_1_ = Es_1_ * Er_by_Es_1_;
    Hs_1_ = H_1_ / ( 1 + Hr_by_Hs_1_ );
    Hr_1_ = Hs_1_ * Hr_by_Hs_1_;
    Cr_1_ = thetaC * Er_1_ / ( thetaF * P_1_ );
    Cs_1_ = thetaC * Es_1_ / ( thetaF * P_1_ );
    C_1_ = Cr_1_ + Cs_1_;
    SRD_1_ = gamma * F_1_ / L_1_;
    Dr_1_ = thetaL / thetaF * Er_1_ / SRD_1_;
    I_1_ = K_1_ * ( 1 - ( 1 - deltaK ) / GYTrend_ ) / ( 1 - Phi2 / 2 * ( GYTrend_ - 1 ) ^ 2 );
    M_1_ = kappa * SP_1_ * Z_1_ / P_1_;
    Y_1_ = C_1_ + I_1_ + M_1_;
    YBar_1_ = Y_1_ * P_1_ ^ ( ( 1 + lambda ) / lambda ) * AverageTransportCost_ ^ ( - 1 / lambda );
    P_star_1_ = ( 1 + lambda ) * SP_1_ * ( 1 - xi * Xi_LEAD_ * ( Pi_^(1/lambda) ) * GYBarTrend_ ) / ( 1 - xi * Xi_LEAD_ * ( Pi_^((1+lambda)/lambda) ) * GOmega1Trend_ );

    SPJ1_1_ = P_star_1_ * ( ( 1 - xi * (1-deltaJ) / GZTrend_ ) / ( 1 - xi * (1-deltaJ) * ( GSPTrend_ * Pi_ )^(1/lambda) / GZTrend_  ) )^(-lambda);

    J_1_ = ( SPJ1_1_ * AverageTransportCost_ / P_1_ ) ^ ( 1 / lambda );

    SPJ2_1_ = P_star_1_ * ( ( 1 - xi * (1-deltaJ) / GZTrend_ ) / ( 1 - xi * (1-deltaJ) * ( GSPTrend_ * Pi_ )^((1+lambda)/lambda) / GZTrend_  ) )^(-lambda/(1+lambda));
    varphi1_1_ = (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^(1/lambda) ) * GYBarTrend_ * YBar_1_ / ( 1 - (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^(1/lambda) ) * GYBarTrend_  );
    varphi2_1_ = (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^((1+lambda)/lambda) ) * GOmega1Trend_ * YBar_1_ * SP_1_ / ( 1 - (1-deltaJ) * xi * Xi_LEAD_ * ( Pi_^((1+lambda)/lambda) ) * GOmega1Trend_ );
    pis_1_ = ( YBar_1_ * ( P_star_1_^( - 1/lambda ) - SP_1_ * P_star_1_^( -(1 + lambda)/lambda ) ) + varphi1_1_ * P_star_1_^( - 1/lambda ) - varphi2_1_ * P_star_1_^( - (1+lambda)/lambda ) ) / ( 1 - ( ( (1-deltaJ) * xi * Xi_LEAD_ * GSPTrend_ / ( 1 -  (1-deltaJ) * xi * Xi_LEAD_ * GSPTrend_ ) ) * ( 1 - xi ) / xi ) );

    Residual = zeros( 6, 1 );
    Residual( 1 ) = Xi_LEAD_ * GSRKTrend_ * SRK_1_ / ( 1 - Xi_LEAD_ * GSRKTrend_ * ( 1 - deltaK ) ) - Q_1_;
    Residual( 2 ) = thetaF * Nr_1_ / Er_1_ * W_1_ * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hr_1_ / Nr_1_LAG_ ) ^ ( 1 + nu ) ) - thetaH * ( Hr_1_ / Nr_1_LAG_ ) ^ nu;
    Residual( 3 ) = pis_1_ - phi_ * SP_1_;
    Residual( 4 ) = ZF_1_ + phi_ * ( J_1_ - ( 1 - deltaJ ) * J_1_ / GJTrend_ ) +  J_1_ * SPJ2_1_^(-(1+lambda)/lambda) * YBar_1_  - Z_1_;
    Residual( 5 ) = P_1_ * Cr_1_ + Er_1_ + SRD_1_ * Dr_1_ - W_1_ * Hr_1_;
    Residual( 6 ) = thetaF * Ns_1_ / Es_1_ * W_1_ * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( Hs_1_ / Ns_1_LAG_ ) ^ ( 1 + nu ) ) - thetaH * ( Hs_1_ / Ns_1_LAG_ ) ^ nu;
end
