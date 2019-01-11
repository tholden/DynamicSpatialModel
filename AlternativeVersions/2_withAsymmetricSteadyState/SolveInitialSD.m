function SD_0_ = SolveInitialSD( psi1, psi2, psi3, N_0_LAG, N_LAG , SN_0_, dBar, d )
    SD_0_ = exp( fzero( @( log_SD_0_ ) GetResidual( exp( log_SD_0_ ), psi1, psi2, psi3, N_0_LAG, N_LAG , SN_0_, dBar, d ), -8.2, optimset( 'Display', 'iter' ) ) );
end

function Residual = GetResidual( SD_0_, psi1, psi2, psi3, N_0_LAG, N_LAG , SN_0_, dBar, d )
    SpatialNumPoints = length( d );
    Residual = 0;
    for i = 1 : SpatialNumPoints
        SN_0_i_ =  ( N_0_LAG / N_LAG ) * psi3 / ( psi1 / ( N_0_LAG - SN_0_ ) + psi2 * ( d( i ) * SN_0_ - SD_0_ ) / ( dBar * SN_0_ * SN_0_ - SN_0_ * SD_0_ ) );
        Residual = Residual + d( i ) * SN_0_i_;
    end
    Residual = Residual / SpatialNumPoints;
    Residual = Residual - SD_0_;
end
