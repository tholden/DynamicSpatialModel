disp( 'Note, plotting code requires subtightplot from: https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot' );

SimulationLength = 1000;
SpatialPointsPerDimension = 100;

Gap = [ 0.05, 0.05 ];

XSimulation = ( 1:SimulationLength ) / 4;
[ SurfaceX, SurfaceY ] = meshgrid( XSimulation, ( 0 : SpatialPointsPerDimension ) / SpatialPointsPerDimension );

EndoSimulation = oo_.endo_simul( :, ( end-SimulationLength+1):end );

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

Axis1 = subtightplot( 3, 1, 1, Gap );
Indices = find( all( bsxfun( @eq, M_.endo_names( :, 1:6 ), 'log_A_' ), 2 ) );
CurrentSurface = [ EndoSimulation( Indices, : ); EndoSimulation( Indices( 1 ), : ) ];
yyaxis left;
pcolor( SurfaceX, SurfaceY, CurrentSurface );
shading interp;
ColourBar1 = colorbar;
yyaxis right;
plot( XSimulation, EndoSimulation( find( all( bsxfun( @eq, M_.endo_names( :, 1:6 ), 'log_A ' ), 2 ), 1 ), : ), 'Color', 'k', 'LineWidth', 1 );
title( 'Log productivity (left) and aggregate log productivity (right)' );

Axis3 = subtightplot( 3, 1, 2, Gap );
Indices = find( all( bsxfun( @eq, M_.endo_names( :, 1:6 ), 'log_C_' ), 2 ) );
CurrentSurface = [ EndoSimulation( Indices, : ); EndoSimulation( Indices( 1 ), : ) ];
yyaxis left;
pcolor( SurfaceX, SurfaceY, CurrentSurface );
shading interp;
ColourBar3 = colorbar;
yyaxis right;
plot( XSimulation, EndoSimulation( find( all( bsxfun( @eq, M_.endo_names( :, 1:6 ), 'log_C ' ), 2 ), 1 ), : ), 'Color', 'k', 'LineWidth', 1 );
title( 'Log consumption (left) and log aggregate consumption (right)' );

Axis2 = subtightplot( 3, 1, 3, Gap );
Indices = find( all( bsxfun( @eq, M_.endo_names( :, 1:8 ), 'level_B_' ), 2 ) );
CurrentSurface = [ EndoSimulation( Indices, : ); EndoSimulation( Indices( 1 ), : ) ];
pcolor( SurfaceX, SurfaceY, CurrentSurface );
shading interp;
ColourBar2 = colorbar;
title( 'Bond holdings' );

MinWidth = min( [ Axis1.Position( 3 ) Axis2.Position( 3 ) Axis3.Position( 3 ) ] );

Axis1.Position( 3 ) = MinWidth;
Axis2.Position( 3 ) = MinWidth;
Axis3.Position( 3 ) = MinWidth;

MaxLeft = max( [ ColourBar1.Position( 1 ) ColourBar2.Position( 1 ) ColourBar3.Position( 1 ) ] );

ColourBar1.Position( 1 ) = MaxLeft;
ColourBar2.Position( 1 ) = MaxLeft;
ColourBar3.Position( 1 ) = MaxLeft;
