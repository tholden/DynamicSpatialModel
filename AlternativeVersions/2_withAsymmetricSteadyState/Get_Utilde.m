function Utilde_x = Get_Utilde( SpatialPointsPerDimension, index, homotopy )

% figure;
% SpatialPointsPerDimension=3;

sigma = 0.01;
%scale = 0.26328;
%scale = 0.3;
%scale = 0.2637;
scale = 0.2;
center = [round(SpatialPointsPerDimension/2) round(SpatialPointsPerDimension/2)];

gsize = [SpatialPointsPerDimension SpatialPointsPerDimension];
[x,y] = ndgrid(1:gsize(1), 1:gsize(2));
Utilde = (exp(-((((x-center(1))/SpatialPointsPerDimension).^2 + ((y-center(2))/SpatialPointsPerDimension).^2)./(2*sigma))));
Utilde = scale * Utilde;
Utilde = 1+Utilde-mean(mean(Utilde));

% surf(Utilde)
% ratio = max(max(Utilde))/min(min(Utilde))
% mean(mean(Utilde))

Utilde_x = 1 + homotopy*(Utilde( index )-1);

end
