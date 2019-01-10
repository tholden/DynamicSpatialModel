function SD_x_ = GetSD_x_( Index , SpatialPointsPerDimension )
    % Only for shape-type Plane!
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
%     i = ij; j = ij;
%     Index = j * SpatialPointsPerDimension - ( SpatialPointsPerDimension - i);
    SNxx_ = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( Index , 7:end );
    Weight = 1 / (SpatialPointsPerDimension^2);
    d = getDistanceMatrix( SpatialPointsPerDimension );
    integral_SD = zeros( 1,SpatialPoints );
    for tildex=1:SpatialPoints
        integral_SD( tildex ) = Weight * d( x , tildex ) * SNxx_( x , tildex  );
    end
    SD_x_ = sum( integral_SD );
end
