function H_x_ = GetH_x_( Index )
    % Only for shape-type Plane!
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
%     i = ij; j = ij;
%     Index = j * SpatialPointsPerDimension - ( SpatialPointsPerDimension - i);
    H_x_ = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( Index , 5 );
end
