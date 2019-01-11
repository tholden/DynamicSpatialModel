function PatchDynareCholToSqrtm
    FilesToPatch = { 'conditional_variance_decomposition.m', 'disp_moments.m', 'extended_path_initialization.m', 'setup_stochastic_perfect_foresight_model_solver.m', 'solve_stochastic_perfect_foresight_model.m', 'forcst2.m', 'imcforecast.m', 'irf.m', 'PosteriorIRF_core1.m', 'simult.m', 'simul_backward_model.m', 'simul_backward_nonlinear_model.m', 'stoch_simul.m', 'th_autocovariances.m' };
    for FileIdx = 1 : length( FilesToPatch )
        FileName = which( FilesToPatch{ FileIdx } );
        disp( [ 'Patching: ' FileName ] );
        FileText = fileread( FileName );
        FileText = regexprep( FileText, '\<chol\>', 'sqrtm' );
        FileHandle = fopen( FileName, 'w' );
        fprintf( FileHandle, '%s', FileText );
        fclose( FileHandle );
    end
end
