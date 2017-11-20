function CombineIRFPlots

    disp( 'Note, this code requires subtightplot from: https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot' );
    
    addpath( cd );

    VariableNames = { 'C', 'K', 'I', 'E', 'F', 'Q', 'J', 'L', 'H', 'N', 'SN', 'SD', 'muN', 'U' };
    VariableLabels = { 'C_{x,t}', 'K_{x,t}', 'I_{x,t}', 'E_{x,t}', 'F_{x,t}', 'Q_{x,t}', 'J_{x,t}', 'L_{x,t}', 'H_{x,t}', 'N_{x,t}', '\mathcal{N}_{x,t}', '\mathcal{D}_{x,t}', '\mu_{N,x,t}', 'U_{x,t}' };

    K = length( VariableNames );

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

    for VariableIdx = 1 : length( VariableNames )
        VariableName = VariableNames{ VariableIdx };
        VariableLabel = VariableLabels{ VariableIdx };
        
        FolderName = [ VariableName '_epsilon_AT_5_5' ];
        cd( FolderName );
        
        cd( 'Frames' );
        ProcessFigure( '1.fig', NewFigure, K, VariableIdx, 1, VariableLabel );
        ProcessFigure( '33.fig', NewFigure, K, VariableIdx, 2, VariableLabel );
        ProcessFigure( '321.fig', NewFigure, K, VariableIdx, 3, VariableLabel );
        ProcessFigure( '3193.fig', NewFigure, K, VariableIdx, 4, VariableLabel );
        cd( '..' );
        ProcessFigure( '_5_5.fig', NewFigure, K, VariableIdx, 5, VariableLabel );
        ProcessFigure( '_1_1.fig', NewFigure, K, VariableIdx, 6, VariableLabel );
        ProcessFigure( 'Aggregate.fig', NewFigure, K, VariableIdx, 7, VariableLabel );
        
        cd( '..' );
    end

end

function ProcessFigure( FileName, NewFigure, NumVariables, VariableIdx, Column, VariableLabel )
    Gap = [ ];
    OldFigure = openfig( FileName, 'new' );
    OldAxes = gca;
    figure( NewFigure );
    NewAxes = subtightplot( NumVariables, 7, (VariableIdx-1)*7 + Column, Gap );
    Children = OldAxes.Children;
    for i = 1 : numel( Children )
        copyobj( Children( i ), NewAxes );
        hold on;
    end
    if VariableIdx < NumVariables
        NewAxes.XTickLabel = [];
    end
    if ( Column >= 2 ) && ( Column <= 4 )
        NewAxes.YTickLabel = [];
    end
    axis( NewAxes, 'square' );
    if Column == 1
        ylabel( NewAxes, [ '{\boldmath\fontsize{18}{24}\selectfont{$$' VariableLabel '$$}}' ], 'Interpreter', 'latex' );
    end
    if VariableIdx == 1
        Titles = { 'Initial impact', '1 year impact', '10 year impact', '100 year impact', '(0.5,0.5)-IRF', '(0,0)-IRF', 'Aggregate IRF' };
        title( NewAxes, Titles{ Column } );
    end
    close( OldFigure );
end
