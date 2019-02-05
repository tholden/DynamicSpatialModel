disp( 'This code requires ffmpeg.exe from http://ffmpeg.zeranoe.com/builds/ to be placed in the same folder.' );

VariableNames = { 'C', 'K', 'I', 'E', 'F', 'N' };

SimulationLength = 4000;

InterpolationAmount = 3;

SpatialDimensions = 2;
SpatialPointsPerDimension = 7;

SpatialNumPoints = SpatialPointsPerDimension ^ SpatialDimensions;

SpatialIndices = cell( 1, SpatialDimensions );
[ SpatialIndices{:} ] = ndgrid( 1:SpatialPointsPerDimension );
SpatialIndices = cellfun( @( c ) c(:), SpatialIndices, 'UniformOutput', false );
SpatialIndices = cell2mat( SpatialIndices );

GridSizeCell = repmat( { SpatialPointsPerDimension }, 1, SpatialDimensions );

XSimulation = ( 1:SimulationLength ) / 4;
ZeroSimulation = zeros( SimulationLength, 1 );

if SpatialDimensions == 2
    InterpolationMultiplier = 2 ^ InterpolationAmount;
    [ SurfaceX, SurfaceY ] = meshgrid( ( 0 : ( SpatialPointsPerDimension * InterpolationMultiplier ) ) / ( SpatialPointsPerDimension * InterpolationMultiplier ) );
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
FigureHandle.Color = [ 1 1 1 ];
colormap( parula( 2 ^ 16 ) );

EndoNames = cellstr( M_.endo_names );
PercentSimulations = 100 * bsxfun( @minus, oo_.endo_simul( :, 1:SimulationLength ), oo_.steady_state );

FolderName = 'Results';
[ mkdirStatus, ~, ~ ] = mkdir( FolderName );
assert( mkdirStatus == 1 );
cd( FolderName );

for VariableIdx = 1 : length( VariableNames )
    VariableName = VariableNames{ VariableIdx };

    FolderName = [ VariableName '_Simulation' ];
    [ mkdirStatus, ~, ~ ] = mkdir( FolderName );
    assert( mkdirStatus == 1 );
    cd( FolderName );

    EndoVariableIndex = find( strcmp( [ 'log_' VariableName ], EndoNames ), 1 );
    if ~isempty( EndoVariableIndex )
        CurrentPercentSimulation = PercentSimulations( EndoVariableIndex, : );
        plot( XSimulation, CurrentPercentSimulation, 'k-' );
        hold on;
        plot( XSimulation, ZeroSimulation, 'r--' );
        hold off;
        drawnow;
        savefig( FigureHandle, 'Aggregate', 'compact' );
        saveas( FigureHandle, 'Aggregate', 'meta' );
    end

    PercentSimulation = zeros( GridSizeCell{:}, SimulationLength );
    for Point = 1 : SpatialNumPoints
        CurrentIndices = SpatialIndices( Point, : );
        CurrentIndicesCell = num2cell( CurrentIndices );
        CurrentIndicesString = sprintf( repmat( '_%d', 1, SpatialDimensions ), CurrentIndices );
        EndoVariableIndex = find( strcmp( [ 'log_' VariableName CurrentIndicesString ], EndoNames ), 1 );
        assert( ~isempty( EndoVariableIndex ) );
        PercentSimulation( CurrentIndicesCell{:}, : ) = PercentSimulations( EndoVariableIndex, : );
    end

    if SpatialDimensions == 2

        Video = VideoWriter( 'Video', 'Archival' ); %#ok<TNMLP>
        Video.FrameRate = 40;
        % Video.Quality = 100;
        open( Video );

        [ mkdirStatus, ~, ~ ] = mkdir( 'Frames' );
        assert( mkdirStatus == 1 );
        cd( 'Frames' );

        PercentSimulation = [ PercentSimulation( :, :, : ), PercentSimulation( :, 1, : ); PercentSimulation( 1, :, : ), PercentSimulation( 1, 1, : ) ];
        minPercentSimulation = min( PercentSimulation(:) );
        maxPercentSimulation = max( PercentSimulation(:) );

        for Period = 1 : size( PercentSimulation, 3 )
            CurrentSurface = interp2( PercentSimulation( :, :, Period ), InterpolationAmount );
            pcolor( SurfaceX, SurfaceY, CurrentSurface );
            shading interp;
            axis square;
            caxis( [ minPercentSimulation, maxPercentSimulation ] );
            colorbar;
            drawnow;

            if ( mod( Period - 1, 40 ) == 0 ) || ( Period == size( PercentSimulation, 3 ) )
                FileName = num2str( Period );
                savefig( FigureHandle, FileName, 'compact' );
                saveas( FigureHandle, FileName, 'png' );
            end

            % Axis = gca;
            % Axis.Units = 'pixels';
            % AxisPosition = ceil( Axis.Position );
            % Margin = max( ceil( Axis.TightInset ) );
            % Rectangle = [ -Margin, -Margin, AxisPosition(3) + 2 * Margin, AxisPosition(4) + 2 * Margin ];
            % Frame = getframe( Axis, Rectangle );

            Frame = getframe( FigureHandle );

            writeVideo( Video, Frame );
        end

        cd( '..' );

        close( Video );

        dos( '..\ffmpeg.exe -i Video.mj2 -filter "setsar=sar=1" -crf 15 -preset slow Video.mp4', '-echo' );
        delete Video.mj2

    end

    cd( '..' );
end

cd( '..' );
