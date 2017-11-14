function [ NewVariableOrder, NewMLVOrder ] = VariableRotationFunction( StatePreShockTotal )
    global M_
    EndoNames = M_.EndoNames;
    NIndices = find( all( EndoNames( :, 1:6 ) == 'log_N_', 2 ) );
    Suffixes = EndoNames( NIndices, 5:end );
    Dimension = length( find( Suffixes( 1, : ) == '_' ) );
    Indexes = sscanf( Suffixes, repmat( '_%d', 1, Dimension ), [ Dimension, size( Suffixes, 1 ) ] )';
    [ ~, NMaxIndex ] = max( StatePreShockTotal( NIndices ) );
    NewIndexes = mod( bsxfun( @minus, Indexes, Indexes( NMaxIndex, : ) ), max( Indexes(:) ) ) + 1;
    NewSuffixes = fprintf( [ repmat( '_%d', 1, Dimension ) '\n' ], NewIndexes );
    NewVariableOrder = 1 : size( EndoNames, 1 );
    
end
