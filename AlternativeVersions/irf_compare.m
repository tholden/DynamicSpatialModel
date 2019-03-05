clear; %close all;

opts.results_files = {'highlambda_highomega_withIRF' ; 'highlambda_lowomega_withIRF' ; 'lowlambda_highomega_withIRF' ; 'lowlambda_lowomega_withIRF'};
opts.model_names = {'high lambda high omega' ; 'high lambda low omega' ; 'low lambda high omega' ; 'low lambda low omega'};
opts.results_files = {'model_2_thetaN4_PhiL2_Phi24_Omega2.3191_zeta8_lambda0.1' ; 'model_2_thetaN1.598_PhiL2_Phi24_Omega1.7_zeta8_lambda0.1' ; 'model_2_thetaN2.5369_PhiL2_Phi24_Omega2.3191_zeta8_lambda0.05' ; 'model_2_thetaN0.7994_PhiL2_Phi24_Omega1.7_zeta8_lambda0.05'};
opts.results_files = {'model_2_irf_thetaN0.7994_PhiL2_Phi24_Omega1.7_zeta8_lambda0.05_AT0.0125' ; 'model_2_irf_thetaN0.7994_PhiL0_Phi24_Omega1.7_zeta8_lambda0.05_AT0.0125' ; 'model_2_irf_thetaN0.7994_PhiL4_Phi24_Omega1.7_zeta8_lambda0.05_AT0.0125'};
opts.model_names = {'Baseline' ; 'Land adj costs = 0' ; 'High Land adj costs = 4'};
opts.results_files = {'model_2_irf_thetaN0.7994_PhiL4_Phi24_Omega1.7_zeta8_lambda0.05_AT0.047'};
opts.model_names = {'High phiL'};
opts.shock = { 'epsilon_AT_4_4' };

% Variables to plot (file 1 and 2)
opts.vars(:,1) = { 'log_Lx' ; 'log_Lx' ; 'log_Cx' ; 'log_Cx' ; 'log_I'  ; 'log_Ix' ; 'log_Nx' ; 'log_Nx' ; 'log_Fx' ; 'log_Fx' ; 'log_Ex' ; 'log_Ex' ; 'log_Wx' ; 'log_Wx' ; 'relWx' ; 'relWx' };
opts.vars(:,2) = { 'log_Lx' ; 'log_Lx' ; 'log_Cx' ; 'log_Cx' ; 'log_I'  ; 'log_Ix' ; 'log_Nx' ; 'log_Nx' ; 'log_Fx' ; 'log_Fx' ; 'log_Ex' ; 'log_Ex' ; 'log_Wx' ; 'log_Wx' ; 'relWx' ; 'relWx'  };
opts.vars(:,3) = { 'log_Lx' ; 'log_Lx' ; 'log_Cx' ; 'log_Cx' ; 'log_I'  ; 'log_Ix' ; 'log_Nx' ; 'log_Nx' ; 'log_Fx' ; 'log_Fx' ; 'log_Ex' ; 'log_Ex' ; 'log_Wx' ; 'log_Wx' ; 'relWx' ; 'relWx'  };
opts.vars(:,4) = { 'log_Lx' ; 'log_Lx' ; 'log_Cx' ; 'log_Cx' ; 'log_I'  ; 'log_Ix' ; 'log_Nx' ; 'log_Nx' ; 'log_Fx' ; 'log_Fx' ; 'log_Ex' ; 'log_Ex' ; 'log_Wx' ; 'log_Wx' ; 'relWx' ; 'relWx'  };

opts.location = {'1_1' ; '4_4' ; '1_1' ; '4_4' ; '1_1' ; '4_4' ; '1_1' ; '4_4' ; '1_1' ; '4_4' ; '1_1' ; '4_4' ; '1_1' ; '4_4' ; '1_1' ; '4_4'}; 
opts.index = [1 ; 25 ; 1 ; 25 ; 1 ; 25 ; 1 ; 25 ; 1 ; 25 ; 1 ; 25 ; 1 ; 25 ; 1 ; 25 ]; 
opts.per_capita = logical([ 0 ; 0 ; 1 ; 1 ; 1 ; 1 ; 0 ; 0 ; 0 ; 0 ; 1 ; 1 ; 0 ; 0 ; 0 ; 0 ]); 
    
opts.vars_names = { '$log L_{1,1}$' ; '$log L_{4,4}$' ; '$log C_{1,1}$' ; '$log C_{4,4}$' ; '$log I_{1,1}$' ; '$log I_{4,4}$' ; '$log N_{1,1}$' ; '$log N_{4,4}$' ; '$log F_{1,1}$' ; '$log F_{4,4}$' ; '$log E_{1,1}$' ; '$log E_{4,4}$' ; '$log W_{1,1}$' ; '$log W_{4,4}$' ; '$log W_{1,1} - log W$' ; '$log W_{4,4} - log W$' };

opts.SpatialPointsPerDimension = 7;
opts.SpatialPoints = opts.SpatialPointsPerDimension*opts.SpatialPointsPerDimension;
opts.simulation_length = 400;

opts.plot_rows = 4;
opts.plot_cols = 4;

opts.num_models = length(opts.results_files);
opts.num_vars = length(opts.vars);

%% Prepare data
for models=1:opts.num_models 
    load(['Results/',char(opts.results_files(models)),'.mat']);
    
    numpar = M_.param_nbr;
   
    for iter = 1:numpar
        eval(['results.parameters{',num2str(models),'}.',deblank(M_.param_names(iter,:)),' = M_.params(', num2str(iter) ,');']);
    end

    numvar = M_.endo_nbr;
    for iter = 1:numvar
        eval(['results.model{',num2str(models),'}.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_',char(opts.shock),';']);
        eval(['results.SS{',num2str(models),'}.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
    end
    
    mlv_simulation_irf_firstorder;
end

clearvars -except opts results
%% Plot results

opts.plots_per_figure = opts.plot_rows * opts.plot_cols;
opts.num_figures = ceil(opts.num_vars / opts.plots_per_figure);
counts.plot_num = 0;

for figures=1:opts.num_figures
    counts.remaining_plots = opts.num_vars - counts.plot_num;
    if counts.remaining_plots < opts.plots_per_figure
        counts.plots_to_plot = counts.remaining_plots;
    else
        counts.plots_to_plot = opts.plots_per_figure;
    end
    figure;
    for subpl=1:counts.plots_to_plot
        subplot(opts.plot_rows,opts.plot_cols,subpl);
        for models=1:opts.num_models
            if opts.per_capita(subpl)
                eval(['x = results.model{',num2str(models),'}.',char(opts.vars(subpl,models)),'(opts.index(subpl),:) - results.SS{',num2str(models),'}.',char(opts.vars(subpl,models)),'(:,opts.index(subpl))'';']);
                n = results.model{models}.log_Nx(opts.index(subpl),:)-results.SS{models}.log_Nx(:,opts.index(subpl))';
                plot( x-n ); hold on;
                title([opts.vars_names(subpl),'$-log N$'],'Interpreter','Latex')
            else
                eval(['plot(results.model{',num2str(models),'}.',char(opts.vars(subpl,models)),'(opts.index(subpl),:)-results.SS{',num2str(models),'}.',char(opts.vars(subpl,models)),'(:,opts.index(subpl))''); hold on;']);
                title(opts.vars_names(subpl),'Interpreter','Latex')
            end
        end
    end
    legend(opts.model_names)
end
