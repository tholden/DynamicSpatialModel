function W = getWeightMatrix( shape, P )
% Gets integration weight of each element 
% Takes the number of points per dimension, P,  and returns a  P x P matrix


W = ones(P , P);
    
if strcmp('P',shape)
    for i=1:P
        for j=1:P
            if i==1 || i==P
                if j==1 || j== P
                    W( i , j ) = (.5 * (1 / ( P - 1 )) ) ^2;
                else
                    W( i , j ) = (.5 * (1 / ( P - 1 )) ) * (1 / ( P - 1 ));
                end
            else
                if j==1 || j== P
                    W( i , j ) = (.5 * (1 / ( P - 1 )) ) * (1 / ( P - 1 ));
                else
                    W( i , j ) = (1 / ( P - 1 ))^2;
                end
            end
        end
    end
else
    W = ( 1 / ( P * P ) ) * W;
end


end
