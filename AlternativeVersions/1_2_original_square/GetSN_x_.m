function SN_x_ = GetSN_x_( Index , SpatialPointsPerDimension )
    % Only for shape-type Plane!
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
%     i = ij; j = ij;
%     Index = j * SpatialPointsPerDimension - ( SpatialPointsPerDimension - i);
    SNxx_ = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( Index , 7:end );
    Weight = 1 / (SpatialPointsPerDimension^2);
    integral_SN_ = zeros( 1,SpatialPoints );
    for tildex=1:SpatialPoints
        integral_SN_( tildex ) = Weight * ( SNxx_( Index , tildex ) );
    end
    SN_x_ = sum( integral_SN_ );
end
