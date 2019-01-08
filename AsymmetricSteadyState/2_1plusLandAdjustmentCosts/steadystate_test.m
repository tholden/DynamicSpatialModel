


YBar(x) = sum( w(xtilde) * Y(xtilde) * P(xtilde) ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance(x,xtilde) ) );
Pi(x)   = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP(x) ^ ( - 1 / lambda ) * YBar(x);

% Firm entry - nonlinear solver
phi * SP(x) = Pi(x) + ( 1 - deltaJ ) * Xi_LEAD * phi_LEAD * SP_LEAD(x);

% Raw good market clearing 
Z(x) = ZF(x) + phi * ( J(x) - ( 1 - deltaJ ) * J_LAG(x) ) + J(x) * ( ( 1 + lambda ) * SP(x) ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar(x);








% Firm entry 
0 = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP(x) ^ ( - 1 / lambda ) * YBar(x) + ( 1 - deltaJ ) * Xi_LEAD * phi_LEAD * SP_LEAD(x) - phi * SP(x);
0 = lambda ./ ( 1 + lambda ) .* ( 1 + lambda ) .^ ( - 1 ./ lambda ) .* SP_x_.^ ( - 1 ./ lambda ) .* YBar_x_ + ( 1 - deltaJ ).* Xi_LEAD_.* phi_ .* SP_x_ .* GSPTrend_ - phi_.* SP_x_;
% this is correct in non-linear solver, so must be an error in definition
% of one of: YBar_x_ (#), SP_x_ (#), Xi_LEAD_ (#)

% Raw good market clearing 
%mod
J(x) = ( Z(x) - ZF(x) ) /  ( phi * ( 1 - ( 1 - deltaJ ) / GJTrend_  ) + ( ( 1 + lambda ) * SP(x) ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar(x) );
%mod SS
J@{Index1}_ = ( Z_@{Index1}_ - ZF_@{Index1}_ ) / ( phi_ * ( 1 - ( 1 - deltaJ ) / GJTrend_ ) + ( ( 1 + lambda ) * SP_@{Index1}_ ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar_@{Index1}_ );
%SS solver
J_x_ = ( Z_x_ - ZF_x_ ) ./ ( phi_ .* ( 1 - ( 1 - deltaJ ) ./ GJTrend_ ) + ( ( 1 + lambda ) .* SP_x_ ) .^ ( - ( 1 + lambda ) ./ lambda ) .*  YBar_x_ ); 
% This equation matches across files, so must be an error in the definition
% of one of: Z(x) (#), ZF(x) (#), SP_x_ (#), YBar_x_ (#)



YBar(x) = sum( w(xtilde) *             Y(xtilde)         * P(xtilde) ^     ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance(x,xtilde) ) );
YBar_x_( x ) = sum( Weight( tildex ) * Y_x_( tildex ) * ( P_x_( tildex ) )^((1+lambda)/lambda) * exp( -tau_ * d( x , tildex ) / lambda ) );
YBar_@{Index1}_ = sum( weight__ * Y_@{Index2}_ * P_@{Index2}_ ^ ( ( 1 + lambda ) / lambda ) * exp( - tau_ / lambda * Distance@{Index1}@{Index2}_ ) );
% Definition of Ybar matches across files, so either correct, or must be an
% error in defintion of one of: Y(x) (#), P(x) (#)


%mod
SP@{Index1} = ( 1 - gamma ) * F@{Index1} / ZF@{Index1}
%mod SS
SP_@{Index1}_ = ( 1 - gamma ) * F@{Index1}_ / ZF_@{Index1}_;
%SS solver
SP_x_ = ( 1 - gamma ) .* F_x ./ ZF_x_;
% Definition of SP matches across files, so either correct, or must be an
% error in definition of one of: F(x), ZF(x) (#)


%mod
Xi_LEAD = beta * ( N(1) / N_LAG(1) ) * ( E(1) / E_LEAD(1) ) * ( U_LEAD(1) / U(1) ) ^ ( 1 - varsigma );
%mod SS
Xi_LEAD_ = beta_ * GN_ / GFTrend_ * GUTrend_ ^ ( 1 - varsigma );
% This looks fine


%mod
ZF@{Index1} = ( F@{Index1} / L@{Index1} ^ gamma ) ^ ( 1 / ( 1 - gamma ) )
%mod SS
ZF_@{Index1}_ = ( F@{Index1}_ / L@{Index1}_ ^ gamma ) ^ ( 1 / ( 1 - gamma ) )
% SS solver
ZF_x_ = ( F_x ./ L_x_ .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) )
% Definition of ZF matches across files, so either correct, or must be an
% error in definition of one of: F(x)*, L(x)

%mod
Z@{Index1} = ( ( K@{Index1}_LAG ^ alpha * ( A@{Index1} * H@{Index1} ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP@{Index1} / P@{Index1} ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) )
%mod SS
Z_@{Index1}_ = ( ( ( K@{Index1}_ / GYTrend_ ) ^ alpha * ( A_@{Index1}_ * H@{Index1}_ ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP_@{Index1}_ / P_@{Index1}_ ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) )
%SS solver
Z_x_ = ( ( ( K_x ./ GYTrend_ ) .^ alpha .* ( A_x .* H_x ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* SP_x_ ./ P_x_ ) .^ kappa ) .^ ( 1 / ( 1 - kappa ) )
% Definition of Z matches across files, so either correct, or must be an
% error in definition of one of: K(x)*, A(x) (#), SP(x) (#), P(x) (#)


%mod
P@{Index1} = ( 1 + lambda ) * ( sum ( Weight@{Index2} * J@{Index2} * ( SP@{Index2} * exp( tau * Distance@{Index1}@{Index2} ) ) ^ ( - 1 / lambda ) ) ) ^ ( - lambda );
% SS solver (non-linear solver equation)
P_x_ = ( 1 + lambda ) .* sum( integral_Px ).^(-lambda)
