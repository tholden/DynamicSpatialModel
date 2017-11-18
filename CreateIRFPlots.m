ShockNames = { 'epsilon_AT_4_4', 'epsilon_GA' 'epsilon_GN' 'epsilon_tau' 'epsilon_phi' 'epsilon_beta' };
VariableNames = { 'C', 'K', 'I', 'E', 'F', 'Q', 'J', 'L', 'H', 'N', 'SN', 'SD', 'muN', 'U' };
SpecificIndices = [ 1, 1; 4, 4 ];

IRFLength = 400;
ShockScale = 10;

SpatialDimensions = 2;
SpatialPointsPerDimension = 7;

SpatialNumPoints = SpatialPointsPerDimension ^ SpatialDimensions;

SpatialIndices = cell( 1, SpatialDimensions );
[ SpatialIndices{:} ] = ndgrid( 1:SpatialPointsPerDimension );
SpatialIndices = cellfun( @( c ) c(:), SpatialIndices, 'UniformOutput', false );
SpatialIndices = cell2mat( SpatialIndices );

GridSizeCell = repmat( { SpatialPointsPerDimension }, 1, SpatialDimensions );

ZeroIRF = zeros( IRFLength, 1 );

FigureHandle = figure( 1 );
figure( FigureHandle );
robot = java.awt.Robot;
robot.keyPress(java.awt.event.KeyEvent.VK_ALT);      %// send ALT
robot.keyPress(java.awt.event.KeyEvent.VK_SPACE);    %// send SPACE
robot.keyRelease(java.awt.event.KeyEvent.VK_SPACE);  %// release SPACE
robot.keyRelease(java.awt.event.KeyEvent.VK_ALT);    %// release ALT
robot.keyPress(java.awt.event.KeyEvent.VK_X);        %// send X
robot.keyRelease(java.awt.event.KeyEvent.VK_X);      %// release X
pause( 0.1 );
FigureHandle.PaperPositionMode = 'auto';

for ShockIdx = 1 : length( ShockNames )
    ShockName = ShockNames{ ShockIdx };
    for VariableIdx = 1 : length( VariableNames )
        VariableName = VariableNames{ VariableIdx };
        PercentIRF = zeros( GridSizeCell{:}, IRFLength );
        for Point = 1 : SpatialNumPoints
            CurrentIndices = SpatialIndices( Point, : );
            CurrentIndicesCell = num2cell( CurrentIndices );
            CurrentIndicesString = sprintf( repmat( '_%d', 1, SpatialDimensions ), CurrentIndices );
            PercentIRF( CurrentIndicesCell{:}, : ) = 100 * ShockScale * oo_.irfs.( [ 'log_' VariableName CurrentIndicesString '_' ShockName ] )( 1:IRFLength );
        end
        FolderName = [ VariableName '_' ShockName ];
        mkdirStatus = mkdir( FolderName );
        assert( mkdirStatus == 1 );
        cd( FolderName );
        for SpecificIndexIdx = 1 : size( SpecificIndices, 1 )
            CurrentIndices = SpecificIndices( SpecificIndexIdx, : );
            CurrentIndicesCell = num2cell( CurrentIndices );
            CurrentIndicesString = sprintf( repmat( '_%d', 1, SpatialDimensions ), CurrentIndices );
            CurrentPercentIRF = squeeze( PercentIRF( CurrentIndicesCell{:}, : ) );
            plot( 1:IRFLength, CurrentPercentIRF, 'k-' );
            hold on;
            plot( 1:IRFLength, ZeroIRF, 'r--' );
            hold off;
            drawnow;
            savefig( FigureHandle, CurrentIndicesString, 'compact' );
            saveas( FigureHandle, CurrentIndicesString, 'meta' );
        end
        AggregatedVariableName = [ 'log_' VariableName '_' ShockName ];
        if isfield( oo_.irfs, AggregatedVariableName )
            CurrentPercentIRF = 100 * ShockScale * oo_.irfs.( AggregatedVariableName );
            plot( 1:IRFLength, CurrentPercentIRF, 'k-' );
            hold on;
            plot( 1:IRFLength, ZeroIRF, 'r--' );
            hold off;
            drawnow;
            savefig( FigureHandle, 'Aggregate', 'compact' );
            saveas( FigureHandle, 'Aggregate', 'meta' );
        end
        cd( '..' );
    end
end
