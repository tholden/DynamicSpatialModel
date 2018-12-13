function F_x_ = GetF_x_( Index )
    % Only for shape-type Plane!
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
%     i = ij; j = ij;
%     Index = j * SpatialPointsPerDimension - ( SpatialPointsPerDimension - i);
    F_x_ = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( Index , 3 );
end
