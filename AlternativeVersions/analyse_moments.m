clear; %close all;

opts.results_files = 'model_2_sim_thetaN0.7994_PhiL2_Phi24_Omega1.7_zeta8_lambda0.05_AT0.0125';
% opts.results_files = 'model_2_sim_thetaN4_PhiL2_Phi24_Omega2.3191_zeta8_lambda0.1';
% opts.results_files = 'model_2_sim_thetaN1.598_PhiL2_Phi24_Omega1.7_zeta8_lambda0.1';
% opts.results_files = 'model_2_sim_thetaN2.5369_PhiL2_Phi24_Omega2.3191_zeta8_lambda0.05';
opts.SpatialPointsPerDimension = 7;
opts.SpatialPoints = opts.SpatialPointsPerDimension*opts.SpatialPointsPerDimension;
opts.simulation_length = 10000;

%% 
load(['Results/',opts.results_files,'.mat'],'oo_','M_')

mlv_simulation;
simulation_plots;


%% Persistence

% Plot auto-covariance functions
% figure;
% subplot(3,1,1); autocorr(Nx(1,:)); title('ACF N')
% subplot(3,1,2); autocorr(Kx(1,:)); title('ACF K')
% subplot(3,1,3); autocorr(Ix(1,:)); title('ACF I')
% 
% acf_N = autocorr(Nx(1,:));
% acf_K = autocorr(Kx(1,:);
% acf_I = autocorr(Ix(1,:));

%% Moments
% variance
average_wage = mean( exp(log_Hx) .* Wx ./ repmat(sum(exp(log_Hx)),opts.SpatialPoints,1) );
wage_growth = Wx(:,2:end) ./ Wx(:,1:end-1);
average_wage_growth = average_wage(2:end) ./ average_wage(1:end-1);

variance.Nx = var(exp(log_Nx),1,2);
variance.Yx = var(log(Yx),1,2);
variance.Yx = var(log(Yx),1,2);
variance.wage_growth = var(wage_growth,1,2);

variance.average_wage_growth = var(average_wage_growth);


spatial_variance.Nx = var(exp(log_Nx)./mean(exp(log_Nx)));
spatial_variance.Yx = var((Yx)./mean(Yx));

variance.Yx = reshape(variance.Yx,[opts.SpatialPointsPerDimension,opts.SpatialPointsPerDimension]);


plot_moments;

% Aggrgegate variables
variance.Y = var(log(Y));%./mean(Y));

metric.Y = variance.Y;
metric.Yx_1_1 = variance.Yx(1);
metric.Yx_4_4 = variance.Yx(25);
metric.Nx_1_1 = variance.Nx(1);
metric.Nx_4_4 = variance.Nx(25);
metric.wage_growth_1_1 = variance.wage_growth(1);
metric.wage_growth_4_4 = variance.wage_growth(25);
variable_names = { 'log_Y' , 'log_Y_1_1' , 'log_Y_4_4' , 'log_N_1_1' , 'log_N_4_4' , 'wagegrowth_1_1' , 'wagegrowth_4_4' };

disp('**---- Dynamic distribution ---**')
disp( table( metric.Y , metric.Yx_1_1  , metric.Yx_4_4 , metric.Nx_1_1  , metric.Nx_4_4 , metric.wage_growth_1_1 , metric.wage_growth_4_4 , 'VariableNames' , variable_names , 'RowNames' , {'var'} ));


spatial_metric.Nx = mean(spatial_variance.Nx);
spatial_metric.Yx = mean(spatial_variance.Yx);
variable_names = { 'Nx' , 'Yx' };

disp('**---- Spatial distribution ---**')
disp( table( spatial_metric.Nx , spatial_metric.Yx , 'VariableNames' , variable_names , 'RowNames' , {'var'} ));

figure;
hist(400*(average_wage_growth-1)); title('Growth rate in average wage');
    xlabel('% oty')

figure;
for t=1:9
    subplot(3,3,t); hist(400*(wage_growth(:,t+5000)-1));
    title(['Period t = ',num2str(t+5000)])
    xlabel('% oty')
end
sgtitle('Spatial wage growth distribution')

   