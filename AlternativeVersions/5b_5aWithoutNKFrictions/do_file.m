clear; %close all;
dbstop if error

run_indet = 0;
run_simult = 0;

param_thetaN = 14;
param_PhiL = 4;
save('param_vals.mat')

%% Solve model
if ~logical(run_indet) && ~logical(run_simult)

dynare DynamicSpatialModel
save('../Results/model_5b.mat')

delete('param_vals.mat')
end

%% Indeterminacy
if logical(run_indet)
    
% Choose 2 parameters
param1_name = 'param_PhiL';
param1_string='\Phi_{L}';
param1_min = 0;
param1_step = 1;
param1_max = 4;
param2_name = 'param_thetaN';
param2_string='\theta_N';
param2_min = 5;
param2_step = 1;
param2_max = 10;

%
param1_vec = [param1_min : param1_step : param1_max];
param2_vec = [param2_min : param2_step : param2_max];
ind = zeros(length(param1_vec),length(param2_vec));
iter.current = 1; iter.total = length(param1_vec) * length(param2_vec);
save('indet.mat','ind','iter')

for jj=1:length(param1_vec)
    for ii=1:length(param2_vec)
    disp(['Iteration ',num2str(iter.current),' of ',num2str(iter.total)]);
    
    eval([param1_name,'= param1_vec(jj);']);
    eval([param2_name,'= param2_vec(ii);']);
    save('param_vals.mat','-append',param1_name,param2_name)    
    %try
        dynare DynamicSpatialModel 
    %catch
     %   disp('Failure...')
    %end
        
    %the following corrected for 4.3.3
    qz_criterium = 1.000001;
    %no of eigenvalue(s) larger than 1 in modulus - 4.3.3
    n_explod = nnz(abs(oo_.dr.eigval) > qz_criterium);
    %new definition for nfwrd for 4.3.3 - M_.endo_nbr-M_.nstatic-M_.npred
    nfwrd = M_.endo_nbr-M_.nstatic-M_.npred;
    
    load('indet.mat')
    if nfwrd == n_explod     %%IF # OF NON-PREDETERMINED VARIABLES == RANK -->UNIQUE EQUILIBRIUM, SADDLE PATH STABLE
       ind(jj,ii)=0;
    elseif nfwrd < n_explod   %%IF # OF NON-PREDETERMINED VARIABLES  < RANK -->INSTABILITY  
       ind(jj,ii)=3;        
    elseif nfwrd > n_explod   %%IF # OF NON-PREDETERMINED VARIABLES  > RANK -->INDETERMINACY
       ind(jj,ii)=4;
    end
    iter.current = iter.current+1;
    save('indet.mat','ind','iter')
    
    end
end
   
figure;
markersize=1;
spy(ind(:,:)==0,'r.');
hold on 
spy(ind(:,:)==4,'k.');
hold on 
spy(ind(:,:)==3,'g.');
axis xy;
set(gca,'XTick',1:1:length(param2_vec));
set(gca,'XTicklabel',{param2_min:param2_step:param2_max});
set(gca,'YTick',1:1:length(param1_vec));
set(gca,'YTicklabel',{param1_min:param1_step:param1_max});
title('Indeterminacy','fontsize',9); 
xlabel(param2_string,'fontsize',9);
ylabel(param1_string,'fontsize',9);    
legend('Determinacy','Indeterminacy','Instability');


end

%% Simulation model given initial values -- MP IRF with capital and housing stock at 1,1 20% above SS
if logical(run_simult)

dynare DynamicSpatialModel

%define options
irf_periods=100;
drop_periods=0; %burnin
irf_replication=1000;

impulse_vec=zeros(1,M_.exo_nbr); 
impulse_vec(strmatch('epsilon_MP',M_.exo_names,'exact'))=.01; 

starting_point=oo_.dr.ys; %define starting point at steady state
starting_point(strmatch('log_K_1_1',M_.endo_names,'exact'))=log(1.1) + starting_point(strmatch('log_K_1_1',M_.endo_names,'exact')); %set capital to 10% above steady state
starting_point(strmatch('log_L_1_1',M_.endo_names,'exact'))=log(0.9) + starting_point(strmatch('log_L_1_1',M_.endo_names,'exact')); %set land to 10% below steady state

% initialize shock matrices to 0
shocks_baseline = zeros(irf_periods+drop_periods,M_.exo_nbr); %baseline
shocks_impulse = shocks_baseline;

%  eliminate shocks with 0 variance
i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 )); %finds shocks with 0 variance
nxs = length(i_exo_var); %number of those shocks
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));%get Cholesky of covariance matrix to generate random numbers

IRF_mat=NaN(M_.endo_nbr,irf_periods,irf_replication);
for irf_iter = 1: irf_replication
    shocks_baseline(:,i_exo_var) = randn(irf_periods+drop_periods,nxs)*chol_S; %generate baseline shocks
    shocks_impulse = shocks_baseline; %use same shocks in impulse simulation
    shocks_impulse(drop_periods+1,:) = shocks_impulse(drop_periods+1,:)+impulse_vec; %add deterministic impulse
    y_baseline = simult_(starting_point,oo_.dr,shocks_baseline,options_.order); %baseline simulation
    y_shock = simult_(starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
    IRF_mat(:,:,irf_iter) = (y_shock(:,M_.maximum_lag+drop_periods+1:end)-y_baseline(:,M_.maximum_lag+drop_periods+1:end)); %add up difference between series
end
IRF_spatial=mean(IRF_mat,3); %take average

save('../Results/model_5_MP_uneven.mat')
%define options
irf_periods=400;
drop_periods=0; %burnin
irf_replication=1000;

impulse_vec=zeros(1,M_.exo_nbr); 
impulse_vec(strmatch('epsilon_MP',M_.exo_names,'exact'))=.01; 

starting_point=oo_.dr.ys; %define starting point at steady state

% initialize shock matrices to 0
shocks_baseline = zeros(irf_periods+drop_periods,M_.exo_nbr); %baseline
shocks_impulse = shocks_baseline;

%  eliminate shocks with 0 variance
i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 )); %finds shocks with 0 variance
nxs = length(i_exo_var); %number of those shocks
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));%get Cholesky of covariance matrix to generate random numbers

IRF_mat=NaN(M_.endo_nbr,irf_periods,irf_replication);
for irf_iter = 1: irf_replication
    shocks_baseline(:,i_exo_var) = randn(irf_periods+drop_periods,nxs)*chol_S; %generate baseline shocks
    shocks_impulse = shocks_baseline; %use same shocks in impulse simulation
    shocks_impulse(drop_periods+1,:) = shocks_impulse(drop_periods+1,:)+impulse_vec; %add deterministic impulse
    y_baseline = simult_(starting_point,oo_.dr,shocks_baseline,options_.order); %baseline simulation
    y_shock = simult_(starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
    IRF_mat(:,:,irf_iter) = (y_shock(:,M_.maximum_lag+drop_periods+1:end)-y_baseline(:,M_.maximum_lag+drop_periods+1:end)); %add up difference between series
end
IRF_spatial=mean(IRF_mat,3); %take average

save('../Results/model_5_MP_even.mat')

delete('param_vals.mat')
    
    
end

%%
clear;
