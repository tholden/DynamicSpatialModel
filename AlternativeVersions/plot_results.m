clear;% close all;

opt.result_files = { 'model_5_uneven' , 'model_5' };% 'model_5_MP_even' };
opt.result_names = { 'Uneven starting pt' , 'Symmetric starting pt' };
opt.shock = 'epsilon_AT_2_2';
opt.irf_horizon = 100;

opt.vars_to_plot = { 'log_C','log_C_2_2','log_C_1_1','C_by_N',...
    'C_by_N_2_2','C_by_N_1_1','log_I','log_I_2_2','log_I_1_1',...
    'log_H','log_H_2_2','log_H_1_1','log_F','log_E_2_2','log_E_1_1',...
    'E_by_N','E_by_N_2_2','E_by_N_1_1','log_L','log_L_2_2',...
    'log_L_1_1','log_SN','log_N_2_2','log_N_1_1','res_by_N',...
    'res_by_N_2_2','res_by_N_1_1','log_F','log_F_2_2',...
    'log_F_1_1','log_U','log_U_2_2','log_U_1_1','log_ZF_2_2',...
    'log_ZF_1_1','log_J','log_J_2_2','log_J_1_1' };
opt.vars_names = { 'log C','log C_2_2','log C_1_1','log C by N',...
    'log C by N_2_2','log C by N_1_1','log I','log I_2_2','log I_1_1',...
    'log H','log H_2_2','log H_1_1','log E','log E_2_2','log E_1_1',...
    'log E by N','log E by N_2_2','log E by N_1_1','log L','log L_2_2',...
    'log L_1_1','log SN','log N_2_2','log N_1_1','log res by N',...
    'log res by N_2_2','log res by N_1_1','log F','log F_2_2',...
    'log F_1_1','log U','log U_2_2','log U_1_1','log ZF_2_2',...
    'log ZF_1_1','log J','log J_2_2','log J_1_1' };
    
opt.fig_rows = 4;
opt.fig_cols = 4;

opt.num_models = length(opt.result_files);
results.model = cell(1,opt.num_models);
results.SS = cell(1,opt.num_models);
results.init = cell(1,opt.num_models);
opt.num_vars = length(opt.vars_names);


for iter=1:opt.num_models
eval(['load(''Results/',char(opt.result_files(iter)),'.mat'')'])

if exist('IRF_spatial','var')
    [results.model{iter},results.SS{iter},results.init{iter}] = save_results(M_,oo_,opt.shock,IRF_spatial,starting_point);
else
    [results.model{iter},results.SS{iter},results.init{iter}] = save_results(M_,oo_,opt.shock);
end

clearvars -except results opt
end

subplot_num = 1;
subplots_per_fig = opt.fig_rows*opt.fig_cols;
figure;
for var_iter=1:opt.num_vars
    if subplot_num > subplots_per_fig
        subplot_num = 1;
        figure;
    end
    subplot(opt.fig_rows,opt.fig_cols,subplot_num);
    for model_iter=1:opt.num_models
        eval(['plot(results.model{1,model_iter}.',char(opt.vars_to_plot(var_iter)),'(1:',num2str(opt.irf_horizon),')); hold on;']);
    end
    title(opt.vars_names(var_iter))
    if subplot_num==subplots_per_fig || var_iter==opt.num_vars
        legend(opt.result_names)
    end
    subplot_num=subplot_num+1;
end
