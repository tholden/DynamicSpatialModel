function SD_x_ = GetSD_x_( Index , SpatialPointsPerDimension )
    global SNxx_
    % change shape if needed:
    shape = 'Torus';
    Weight = getWeightMatrix( shape , SpatialPointsPerDimension );
    d = getDistanceMatrix( shape, SpatialPointsPerDimension ); 
    SpatialPoints = SpatialPointsPerDimension * SpatialPointsPerDimension;
    integral_SD = zeros( 1,SpatialPoints );
    for tildex=1:SpatialPoints
        integral_SD( tildex ) = Weight( tildex ) * d( Index , tildex ) * SNxx_( Index , tildex  );
    end
    SD_x_ = sum( integral_SD );
end
