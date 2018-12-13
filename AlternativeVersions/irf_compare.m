clear; close all;

opts.results_files = {'model_5b' ; 'model_7b' ; 'model_7b'};
opts.model_names = {'Base' ; 'With renters (savers)' ; 'with renters (renters)'};
opts.shock = { 'epsilon_AT_2_2' };

% Variables to plot (file 1 and 2)
opts.vars(:,1) = { 'log_L_1_1' ; 'log_L_2_2' ; 'log_C_1_1' ; 'log_C_2_2' ; 'log_I_1_1'  ; 'log_I_2_2' ; 'log_N_1_1' ; 'log_N_2_2' ; 'log_F_1_1' ; 'log_F_2_2' ; 'log_E_1_1' ; 'log_E_2_2' };
opts.vars(:,2) = { 'log_L_1_1' ; 'log_L_2_2' ; 'log_Cs_1_1' ; 'log_Cs_2_2' ; 'log_I_1_1'  ; 'log_I_2_2' ; 'log_Ns_1_1' ; 'log_Ns_2_2' ; 'log_F_1_1' ; 'log_F_2_2' ; 'log_Es_1_1' ; 'log_Es_2_2' };
opts.vars(:,3) = { 'log_L_1_1' ; 'log_L_2_2' ; 'log_Cr_1_1' ; 'log_Cr_2_2' ; 'log_I_1_1'  ; 'log_I_2_2' ; 'log_Nr_1_1' ; 'log_Nr_2_2' ; 'log_F_1_1' ; 'log_F_2_2' ; 'log_Er_1_1' ; 'log_Er_2_2' };

opts.vars_names = { '$log L_{1,1}$' ; '$log L_{2,2}$' ; '$log C_{1,1}$' ; '$log C_{2,2}$' ; '$log I_{1,1}$' ; '$log I_{2,2}$' ; '$log N_{1,1}$' ; '$log N_{2,2}$' ; '$log F_{1,1}$' ; '$log F_{2,2}$' ; '$log E_{1,1}$' ; '$log E_{2,2}$' };

opts.plot_rows = 3;
opts.plot_cols = 4;

opts.num_models = length(opts.results_files);
opts.num_vars = length(opts.vars);

%% Prepare data
for models=1:opts.num_models 
    eval(['load(''Results/',char(opts.results_files(models)),'.mat'')']);
    
    numpar = M_.param_nbr;
    for iter = 1:numpar
        eval(['results.parameters{',num2str(models),'}.',deblank(M_.param_names(iter,:)),' = M_.params(', num2str(iter) ,');']);
    end

    numvar = M_.endo_nbr;
    for iter = 1:numvar
        eval(['results.model{',num2str(models),'}.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_',char(opts.shock),';']);
        eval(['results.SS{',num2str(models),'}.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
    end

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
            eval(['plot(results.model{',num2str(models),'}.',char(opts.vars(subpl,models)),'); hold on;']);
            title(opts.vars_names(subpl),'Interpreter','Latex')
        end
    end
    legend(opts.model_names)
end
