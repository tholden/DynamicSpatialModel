ShockNames = { 'epsilon_AT_5_5', 'epsilon_GA' 'epsilon_GN' 'epsilon_tau' 'epsilon_phi' 'epsilon_beta' };
VariableNames = { 'C', 'K', 'I', 'E', 'F', 'Q', 'J', 'L', 'H', 'N', 'SN', 'SD', 'muN', 'U' };
SpecificIndices = [ 1, 1; 4, 4 ];

IRFLength = 400;
ShockScale = 10;

SpatialDimensions = 2;
SpatialPointsPerDimension = 8;

SpatialNumPoints = SpatialPointsPerDimension ^ SpatialDimensions;

SpatialIndices = cell( 1, SpatialDimensions );
[ SpatialIndices{:} ] = ndgrid( 1:SpatialPointsPerDimension );
SpatialIndices = cellfun( @( c ) c(:), SpatialIndices, 'UniformOutput', false );
SpatialIndices = cell2mat( SpatialIndices );

GridSizeCell = repmat( { SpatialPointsPerDimension }, 1, SpatialDimensions );

XIRF = ( 1:400 ) / 4;
ZeroIRF = zeros( IRFLength, 1 );

if SpatialDimensions == 2
    [ SurfaceX, SurfaceY ] = meshgrid( ( 0:SpatialPointsPerDimension ) / SpatialPointsPerDimension );
end

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

        FolderName = [ VariableName '_' ShockName ];
        [ mkdirStatus, ~, ~ ] = mkdir( FolderName );
        assert( mkdirStatus == 1 );
        cd( FolderName );

        AggregatedVariableName = [ 'log_' VariableName '_' ShockName ];
        if isfield( oo_.irfs, AggregatedVariableName )
            CurrentPercentIRF = 100 * ShockScale * oo_.irfs.( AggregatedVariableName );
            plot( XIRF, CurrentPercentIRF, 'k-' );
            hold on;
            plot( XIRF, ZeroIRF, 'r--' );
            hold off;
            drawnow;
            savefig( FigureHandle, 'Aggregate', 'compact' );
            saveas( FigureHandle, 'Aggregate', 'meta' );
        end
        
        PercentIRF = zeros( GridSizeCell{:}, IRFLength );
        for Point = 1 : SpatialNumPoints
            CurrentIndices = SpatialIndices( Point, : );
            CurrentIndicesCell = num2cell( CurrentIndices );
            CurrentIndicesString = sprintf( repmat( '_%d', 1, SpatialDimensions ), CurrentIndices );
            PercentIRF( CurrentIndicesCell{:}, : ) = 100 * ShockScale * oo_.irfs.( [ 'log_' VariableName CurrentIndicesString '_' ShockName ] )( 1:IRFLength );
        end
        
        for SpecificIndexIdx = 1 : size( SpecificIndices, 1 )
            CurrentIndices = SpecificIndices( SpecificIndexIdx, : );
            CurrentIndicesCell = num2cell( CurrentIndices );
            CurrentIndicesString = sprintf( repmat( '_%d', 1, SpatialDimensions ), CurrentIndices );
            CurrentPercentIRF = squeeze( PercentIRF( CurrentIndicesCell{:}, : ) );
            plot( XIRF, CurrentPercentIRF, 'k-' );
            hold on;
            plot( XIRF, ZeroIRF, 'r--' );
            hold off;
            drawnow;
            savefig( FigureHandle, CurrentIndicesString, 'compact' );
            saveas( FigureHandle, CurrentIndicesString, 'meta' );
        end
        
        if SpatialDimensions == 2
            
            Video = VideoWriter( 'Video', 'MPEG-4' ); %#ok<TNMLP>
            Video.FrameRate = 4;
            Video.Quality = 90;
            open( Video );

            [ mkdirStatus, ~, ~ ] = mkdir( 'Frames' );
            assert( mkdirStatus == 1 );
            cd( 'Frames' );

            for Period = 1 : IRFLength
                CurrentSurface = [ PercentIRF( :, :, Period ), PercentIRF( :, 1, Period ); PercentIRF( 1, :, Period ), PercentIRF( 1, 1, Period ) ];
                pcolor( SurfaceX, SurfaceY, CurrentSurface );
                shading interp;
                drawnow;
                FileName = num2str( Period );
                savefig( FigureHandle, FileName, 'compact' );
                saveas( FigureHandle, FileName, 'meta' );
                
                % Axis = gca;
                % Axis.Units = 'pixels';
                % A xisPosition = Axis.Position;
                % AxisTightInset = Axis.TightInset;
                % Rectangle = [ -AxisTightInset(1), -AxisTightInset(2), AxisPosition(3)+AxisTightInset(1)+AxisTightInset(3), AxisPosition(4)+AxisTightInset(2)+AxisTightInset(4) ];
                % Frame = getframe( Axis, Rectangle );
                
                Frame = getframe( FigureHandle );
                writeVideo( Video, Frame );
            end
            
            cd( '..' );
            
            close( Video );
            
        end
        
        cd( '..' );
    end
end
