    N = sum( Weight( 1i ) * N( 1i ) );

% loop over i
    A( 1i )  = AP * AT( 1i );
    ZF( 1i ) = ( F( 1i ) /  L(1i)^ gamma ) ^ ( 1 / ( 1 - gamma ) );
SRL( 1i ) = gamma * F( 1i ) / L( 1i );
    SP( 1i ) = ( 1 - gamma ) * F( 1i ) / ZF( 1i );
    K( 1i ) *( 1 - ( 1 - deltaK )/GYTrend ) = ( 1 - Phi2 / 2 * ( GYTrend - 1 ) ^ 2 ) * I( 1i );

% loop over i,j
    SN( 1i , 1j ) = psi3 * N_LAG( 1j ) / N_LAG / ( ( muN( 1i ) - muN( 1j ) ) / ( ( 1 - varsigma ) * N_LAG( 1i ) * U( 1i ) ^ ( 1 - varsigma ) ) + psi1 / ( N_LAG( 1i ) - SN( 1i ) ) + psi2 * ( Distance( 1i , 1j ) * SN( 1i ) - SD( 1i ) ) / ( dBar * SN( 1i ) * SN( 1i ) - SN( 1i ) * SD( 1i ) ) );
      
    
% loop over i
    SN( 1i ) = sum( Weight( 1j ) * SN( 1i , 1j ) );
    SD( 1i ) = sum( Weight( 1j ) * Distance( 1i , 1j ) * SN( 1i , 1j ) );
    U( 1i ) = ( C( 1i ) / N_LAG( 1i )) ^ thetaC * ( E( 1i ) / N_LAG( 1i ) ) ^ thetaF ...
        * ( ( 1 - L( 1i ) ) / N_LAG( 1i ) ) ^ thetaL  * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H( 1i ) / N_LAG( 1i ) ) ^ ( 1 + nu ) ) ^ thetaH ...
        * ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N_LAG( 1i ) / N_LAG ) ^ 2 ) ^ thetaN * ( 1 - SN( 1i ) / N_LAG( 1i ) ) ^ psi1 * ( dBar - SD( 1i ) / SN( 1i ) ) ^ psi2 ...
        * exp( psi3 * sum( Weight( 1j ) * N_LAG( 1j ) / N_LAG * log( SN( 1i , 1j ) / N_LAG( 1i ) ) ) );
    muN( 1i ) = beta * ( muN_LEAD( 1i ) * GN_LEAD + U_LEAD( 1i ) ^ ( 1 - varsigma ) ...
        + ( 1 - varsigma ) * U_LEAD( 1i )^( 1 - varsigma ) ...
        * ( thetaH * ( H_LEAD( 1i ) / N( 1i ) ) ^ ( 1 + nu ) / ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H_LEAD( 1i ) / N( 1i ) ) ^ ( 1 + nu ) ) ...
            - thetaN * log( N( 1i ) / N ) / ( 1 / 2 * Omega ^ 2 - 1 / 2 * log( N( 1i ) / N ) ^ 2 ) ...
            + psi1 * SN_LEAD( 1i ) / ( N( 1i ) - SN_LEAD( 1i ) ) - ( thetaC + thetaF + thetaL + psi3 ) ) );


    Xi_LEAD = beta * ( GN / GFTrend ) * ( GUTrend ) ^ ( 1 - varsigma );
% loop over j > 1
    E( 1 ) / N_LAG( 1 ) / U( 1 ) ^ ( 1 - varsigma ) = E( 1j ) / N_LAG( 1j ) / U( 1j ) ^ ( 1 - varsigma );

% loop over i > 1
    N( 1i ) = GN * N_LAG( 1i ) - SN( 1i ) + sum( Weight( 1j ) * SN( 1j , 1i ) );

% loop over j
    1 = sum( Weight( 1j ) * N( 1j ) );

% loop over i
    P( 1i ) = ( 1 + lambda ) * ( sum( Weight( 1j ) * J( 1j ) * ( SP( 1j ) * exp( tau * Distance( 1i , 1j ) ) ) ^ ( - 1 / lambda ) ) )^ ( - lambda );
    Z( 1i ) = ( ( K_LAG( 1i )^ alpha * ( A( 1i ) * H( 1i ) ) ^ ( 1 - alpha ) ) ^ ( 1 - kappa ) * ( kappa * SP( 1i ) / P( 1i ) ) ^ kappa ) ^ ( 1 / ( 1 - kappa ) );
    M( 1i ) = kappa * SP( 1i ) * Z( 1i ) / P( 1i );
    SRK( 1i ) = ( 1 - kappa ) * alpha * SP( 1i ) * Z( 1i ) / K_LAG( 1i );
    W( 1i ) = ( 1 - kappa ) * ( 1 - alpha ) * SP( 1i ) * Z( 1i ) / H( 1i );
    1 = Xi_LEAD * ( GPTrend * SRK( 1i ) + GPTrend * Q( 1i ) * ( 1 - deltaK ) ) / Q( 1i );
    P( 1i ) = Q( 1i ) * ( 1 - Phi2 / 2 * ( GYTrend - 1 ) ^ 2 - Phi2 * ( GYTrend - 1 ) * GYTrend ) + Xi_LEAD * Q( 1i ) * GPTrend * Phi2 * ( GYTrend - 1 ) * GYTrend^ 2;
    Y( 1i ) = C( 1i ) + I( 1i ) + M( 1i );
    thetaC * E( 1i ) = thetaF * P( 1i ) * C( 1i );
    thetaL * E( 1i ) = thetaF * SRL( 1i ) * ( 1 - L( 1i ) );
    thetaH * ( H( 1i ) / N_LAG( 1i ) )^nu = thetaF * N( 1i ) / E( 1i ) * W( 1i ) * ( 1 / ( 1 + nu ) * Gamma^( 1 + nu ) - 1 / ( 1 + nu ) * ( H( 1i ) / N_LAG( 1i) )^( 1 + nu ) );
%  thetaF * N@{Index1} / E@{Index1} * W@{Index1} * ( 1 / ( 1 + nu ) * Gamma ^ ( 1 + nu ) - 1 / ( 1 + nu ) * ( H@{Index1} / N@{Index1}_LAG ) ^ ( 1 + nu ) ) - thetaH * ( H@{Index1} / N@{Index1}_LAG ) ^ nu;
%   thetaF .* N_x_LAG_ ./ E_x_ .* W_x_ .* ( 1 ./ ( 1 + nu ) .* Gamma .^ ( 1 + nu ) - 1 ./ ( 1 + nu ) .* ( H_x ./ N_x_LAG_ ).^ ( 1 + nu ) ) - thetaH .* ( H_x ./ N_x_LAG_ ).^ nu
    
    
% loop over i
    YBar( 1i ) = sum( Weight( 1j ) * Y( 1j ) * P( 1j ) ^ ( ( 1 + lambda ) / lambda ) * exp( - tau / lambda * Distance( 1i , 1j ) ) );
Pi( 1i ) = lambda / ( 1 + lambda ) * ( 1 + lambda ) ^ ( - 1 / lambda ) * SP( 1i ) ^ ( - 1 / lambda ) * YBar( 1i );
    phi * SP( 1i ) = Pi( 1i ) + ( 1 - deltaJ ) * Xi_LEAD * phi_LEAD * SP_LEAD( 1i );
    Z( 1i ) = ZF( 1i ) + J( 1i ) *  ( phi * ( 1  - ( 1 - deltaJ ) / GZTrend ) +  ( ( 1 + lambda ) * SP( 1i ) ) ^ ( - ( 1 + lambda ) / lambda ) *  YBar( 1i ) );

    1 = R * Xi_LEAD;
    0 = sum( Weight( 1i ) * ( E( 1i ) - F( 1i ) ) );
    C = sum( Weight( 1i ) * ( C( 1i ) ) );
    K = sum( Weight( 1i ) * ( K( 1i ) ) );
    I = sum( Weight( 1i ) * ( I( 1i ) ) );
    F = sum( Weight( 1i ) * ( F( 1i ) ) );
    J = sum( Weight( 1i ) * ( J( 1i ) ) );
    L = sum( Weight( 1i ) * ( L( 1i ) ) );
    H = sum( Weight( 1i ) * ( H( 1i ) ) );
    SN = sum( Weight( 1i ) * ( SN( 1i ) ) );
    SD = sum( Weight( 1i ) * ( SD( 1i ) ) );
    U = sum( Weight( 1i ) * ( U( 1i ) ) );





















