function SD_1_ = GetSD_1_( psi1, psi2, psi3, N_1_LAG_, SN_1_, dBar, varargin )
    SD_1_ = exp( fzero( @( log_SD_1_ ) GetResidual( exp( log_SD_1_ ), psi1, psi2, psi3, N_1_LAG_, SN_1_, dBar, varargin{:} ), -9, optimset( 'Display', 'iter' ) ) );
end

function Residual = GetResidual( SD_1_, psi1, psi2, psi3, N_1_LAG_, SN_1_, dBar, varargin )
    SpatialNumPoints = length( varargin );
    Residual = 0;
    for i = 1 : SpatialNumPoints
        SN_1_i_ = psi3 / ( psi1 / ( N_1_LAG_ - SN_1_ ) + psi2 * ( varargin{i} * SN_1_ - SD_1_ ) / ( dBar * SN_1_ * SN_1_ - SN_1_ * SD_1_ ) );
        Residual = Residual + varargin{i} * SN_1_i_;
    end
    Residual = Residual / SpatialNumPoints;
    Residual = Residual - SD_1_;
end
