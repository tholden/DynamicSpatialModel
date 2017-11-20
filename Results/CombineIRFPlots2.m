function CombineIRFPlots2

    disp( 'Note, this code requires subtightplot from: https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot' );
    
    addpath( cd );

    VariableNames = { 'C', 'K', 'I', 'F', 'J', 'L', 'H', 'SD', 'U' };
    VariableLabels = { 'C_{x,t}', 'K_{x,t}', 'I_{x,t}', 'F_{x,t}', 'J_{x,t}', 'L_{x,t}', 'H_{x,t}', '\mathcal{D}_{x,t}', 'U_{x,t}' };
    
    ShockNames = { 'epsilon_GA' 'epsilon_GN' 'epsilon_tau' 'epsilon_phi' 'epsilon_beta' };

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
    
    NewFigure.PaperPositionMode = 'auto';
    NewFigure.Color = [ 1 1 1 ];
    colormap( parula( 2 ^ 16 ) );

    for VariableIdx = 1 : K
        VariableName = VariableNames{ VariableIdx };
        VariableLabel = VariableLabels{ VariableIdx };
        
        for ShockIdx = 1 : length( ShockNames )
            ShockName = ShockNames{ ShockIdx };
            FolderName = [ VariableName '_' ShockName ];
            cd( FolderName );
            ProcessFigure( 'Aggregate.fig', NewFigure, K, VariableIdx, ShockIdx, VariableLabel );
            cd( '..' );
        end
    end
    
end

function ProcessFigure( FileName, NewFigure, NumVariables, VariableIdx, Column, VariableLabel )
    Gap = [ ];
    OldFigure = openfig( FileName, 'new' );
    OldAxes = gca;
    figure( NewFigure );
    NewAxes = subtightplot( NumVariables, 5, (VariableIdx-1)*5 + Column, Gap );
    Children = OldAxes.Children;
    for i = 1 : numel( Children )
        copyobj( Children( i ), NewAxes );
        hold on;
    end
    axis( NewAxes, 'square' );
    if VariableIdx < NumVariables
        NewAxes.XTickLabel = [];
    end
    if Column == 1
        ylabel( NewAxes, [ '{\boldmath\fontsize{18}{24}\selectfont{$$' VariableLabel '$$}}' ], 'Interpreter', 'latex' );
    end
    if VariableIdx == 1
        Titles = { 'G_A', 'G_N', '\tau', '\phi', '\beta' };
        title( NewAxes, [ '{\boldmath\fontsize{18}{24}\selectfont{$$\epsilon_{{' Titles{ Column } '},t}$$}}' ], 'Interpreter', 'latex' );
    end
    close( OldFigure );
end
