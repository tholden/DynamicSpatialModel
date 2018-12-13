function SNxx_ = GetSNxx_( Index1 , Index2 )
    % Only for shape-type Plane!
    global E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_
    SNxx_ = E_x_by_F_x_N_x_F_x_K_x_H_x_Q_x_SN_xx_( Index , 7:end );
    SNxx_ = SNxx_( Index1 , Index2 );
end
