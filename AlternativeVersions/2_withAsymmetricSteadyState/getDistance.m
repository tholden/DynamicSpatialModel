function d = getDistance( i , j , P, shape )

D = getDistanceMatrix( shape, P );

d = D( i , j );

end
