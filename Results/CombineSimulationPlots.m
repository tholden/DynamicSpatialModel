function CombineSimulationPlots

    disp( 'Note, this code requires subtightplot from: https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot' );
    
    addpath( cd );

    NewFigure = figure;
    figure( NewFigure );
    robot = java.awt.Robot;
    robot.keyPress(java.awt.event.KeyEvent.VK_ALT);      %// send ALT
    robot.keyPress(java.awt.event.KeyEvent.VK_SPACE);    %// send SPACE
    robot.keyRelease(java.awt.event.KeyEvent.VK_SPACE);  %// release SPACE
    robot.keyRelease(java.awt.event.KeyEvent.VK_ALT);    %// release ALT
    robot.keyPress(java.awt.event.KeyEvent.VK_X);        %// send X
    robot.keyRelease(java.awt.event.KeyEvent.VK_X);      %// release X
    pause( 0.1 );
    
    NewFigure.PaperPositionMode = 'auto';
    NewFigure.Color = [ 1 1 1 ];
    colormap( parula( 2 ^ 16 ) );

    cd( 'N_Simulation' );
    cd( 'Frames' );

    for PlotNumber = 1 : 50
        FileNumber = ( PlotNumber - 1 ) * 80 + 1;
        FileName = [ num2str( FileNumber ) '.fig' ];
        
        ProcessFigure( FileName, NewFigure, PlotNumber, [ 't = ' num2str( ( FileNumber - 1 ) / 4 ) ] );
    end
    
    cd( '..' );
    cd( '..' );
    
end

function ProcessFigure( FileName, NewFigure, PlotNumber, PlotTitle )
    Gap = [ ];
    OldFigure = openfig( FileName, 'new' );
    OldAxes = gca;
    figure( NewFigure );
    NewAxes = subtightplot( 10, 5, PlotNumber, Gap );
    Children = OldAxes.Children;
    for i = 1 : numel( Children )
        copyobj( Children( i ), NewAxes );
        hold on;
    end
    caxis( NewAxes, caxis( OldAxes ) );
    if PlotNumber < 46
        NewAxes.XTickLabel = [];
    end
    if mod( PlotNumber, 5 ) ~= 1
        NewAxes.YTickLabel = [];
    end
    axis( NewAxes, 'square' );
    title( NewAxes, PlotTitle );
    close( OldFigure );
end
