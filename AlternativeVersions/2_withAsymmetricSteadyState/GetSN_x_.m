function SN_x = GetSN_x_( Index , SpatialPointsPerDimension )
    global SNxx_
    % change shape if needed:
    Weight = getWeightMatrix( 'Torus' , SpatialPointsPerDimension );
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    integral_SN_ = zeros( 1,SpatialPoints );
    for tildex=1:SpatialPoints
        integral_SN_( tildex ) = Weight( tildex ) * ( SNxx_( Index , tildex ) );
    end
    SN_x = sum( integral_SN_ );
end
