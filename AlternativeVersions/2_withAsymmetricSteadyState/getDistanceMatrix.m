function d = getDistanceMatrix( shape, P )
% Gets distance between all elements in a uniform 2-dimensional plane
% Takes the number of points per dimension, P,  and returns a  P x P matrix
% with the distance between points represented by the index.

X = P * P;
d = zeros( X , X );

for i=1:X
    for j=1:X

        [ i1 , i2 ] = ind2sub([P,P],i);
        [ j1 , j2 ] = ind2sub([P,P],j);

        if strcmp( shape , 'P')
          d( i , j ) = sqrt( ( (j1-1)/(P-1) - (i1-1)/(P-1) )^2 + ( (j2-1)/(P-1) - (i2-1)/(P-1) )^2 );
        else
          d( i , j ) = sqrt( min( abs((j1-1)/(P) - (i1-1)/(P)) , 1 - abs((j1-1)/(P) - (i1-1)/(P)) )^2  + min( abs((j2-1)/(P) - (i2-1)/(P)) , 1 - abs( (j2-1)/(P) - (i2-1)/(P) ) )^2 );
        end
    end
end


end
