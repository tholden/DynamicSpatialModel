function [c,ceq] = constraints_sym( E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_, shape, SpatialPointsPerDimension, d, GN_, nu, gamma, Gamma, Omega, lambda, phi_, deltaJ, GJTrend_, thetaC, thetaF, thetaH, thetaL, thetaN, kappa, alpha, GYTrend_, GSRKTrend_, Xi_LEAD_, deltaK, Phi2, GSPTrend_, GmuNTrend_, GUTrend_, psi1, psi2, psi3, tau_, dBar, beta_, varsigma )

SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;

E_by_F_x    = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 1:SpatialPoints );
N_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( SpatialPoints+1:2*SpatialPoints );
F_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 2*SpatialPoints+1:3*SpatialPoints );
K_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 3*SpatialPoints+1:4*SpatialPoints );
H_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 4*SpatialPoints+1:5*SpatialPoints );
Q_x         = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 5*SpatialPoints+1:6*SpatialPoints );
SN_xx_       = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( 6*SpatialPoints+1:end );
SN_xx_       = reshape(SN_xx_,[SpatialPoints,SpatialPoints]);

num_ceq_it = 1;

ceq( num_ceq_it ) = sum( N_x ) - SpatialPoints;
num_ceq_it = num_ceq_it+1;

for ii=1:length(E_by_F_x)
    ceq( num_ceq_it ) = E_by_F_x( ii ) - 1;
    num_ceq_it = num_ceq_it+1;
end
  
for ii=1:length(N_x)
    ceq( num_ceq_it ) = N_x( ii ) - 1;
    num_ceq_it = num_ceq_it+1;
end

for ii=1:length(K_x)-1
    ceq( num_ceq_it ) = K_x( 1 ) - K_x( ii+1 );
    num_ceq_it = num_ceq_it+1;
end

for ii=1:length(F_x)-1
    ceq( num_ceq_it ) = F_x( 1 ) - F_x( ii+1 );
    num_ceq_it = num_ceq_it+1;
end

for ii=1:length(H_x)-1
    ceq( num_ceq_it ) = H_x( 1 ) - H_x( ii+1 );
    num_ceq_it = num_ceq_it+1;
end

for ii=1:length(Q_x)-1
    ceq( num_ceq_it ) = Q_x( 1 ) - Q_x( ii+1 );
    num_ceq_it = num_ceq_it+1;
end

c = [];
 
end        