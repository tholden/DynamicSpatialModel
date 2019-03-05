clear global; clear;  
%close all;
%dbstop if error
global first_run
first_run=1;

% options
run_indet = 0; % indeterminacy/instability region testing         
run_singlesimulation = 1; % single time-series simulation
run_steadysolver = 1; % stand-alone steady state solver

% Parameters
load('../Results/lowlambda_lowomega_SSinit.mat')
param_thetaN = 0.7994; 
param_Omega = 1.7; 
param_lambda = 0.05; 

% load('../Results/highlambda_highomega_SSinit.mat')
% param_thetaN = 4;  
% param_Omega = 2.31908;
% param_lambda = 0.1; 

% load('../Results/highlambda_lowomega_SSinit.mat')
% param_thetaN = 1.598;  
% param_Omega = 1.7;  
% param_lambda = 0.1; 

% load('../Results/lowlambda_highomega_SSinit.mat')
% param_thetaN = 2.5369;  
% param_Omega = 2.31908;  
% param_lambda = 0.05;



param_Phi2 = 4; %4
param_PhiL = 2; %2
param_zeta = 8; %8
save('param_vals.mat')

%% Solve steady state 
if logical(run_steadysolver)% && ~logical(run_singlesimulation)
global E_by_F_x N_x F_x K_x H_x Q_x par  
SpatialPointsPerDimension = 7; % need to also change in mod file if simulating model
param.thetaN = param_thetaN;
param.Omega=param_Omega;
param.lambda=param_lambda;
param.Phi2=param_Phi2;
param.PhiL=param_PhiL;
param.zeta=param_zeta;
E_by_F_1_ = solve_ss(SpatialPointsPerDimension,param);
E_x = E_by_F_x .* F_x;
Nx = reshape(N_x,[SpatialPointsPerDimension,SpatialPointsPerDimension]);
Ex = reshape(E_x,[SpatialPointsPerDimension,SpatialPointsPerDimension]);
Fx = reshape(F_x,[SpatialPointsPerDimension,SpatialPointsPerDimension]);
Kx = reshape(K_x,[SpatialPointsPerDimension,SpatialPointsPerDimension]);
Hx = reshape(H_x,[SpatialPointsPerDimension,SpatialPointsPerDimension]);
Qx = reshape(Q_x,[SpatialPointsPerDimension,SpatialPointsPerDimension]);

Px_by_Qx = ( 1 - par.Phi2 / 2 * ( par.GYTrend_ - 1 ) ^ 2 - par.Phi2 * ( par.GYTrend_ - 1 ) * par.GYTrend_ ) + par.Xi_LEAD_ * par.GSRKTrend_ * par.Phi2 * ( par.GYTrend_ - 1 ) * par.GYTrend_ ^ 2; 
Px = Px_by_Qx .* Qx;
Cx = par.thetaC .* Ex ./ ( par.thetaF .* Px );

cent_pt = round(SpatialPointsPerDimension/2);
disp(['City vs rural statistics:'])
disp(['Max population ratio: ',num2str(Nx(cent_pt,cent_pt)/Nx(1,1))]);
disp(['Max capital ratio: ',num2str(Kx(cent_pt,cent_pt)/Kx(1,1))]);
disp(['Max hours per head ratio: ',num2str(Hx(cent_pt,cent_pt) * Nx(1,1)/( Hx(1,1) * Nx(cent_pt,cent_pt)))]);
disp(['Max eating per head ratio: ',num2str(Ex(cent_pt,cent_pt) * Nx(1,1)/( Ex(1,1) * Nx(cent_pt,cent_pt)))]);
disp(['Max consumption per head ratio: ',num2str(Cx(cent_pt,cent_pt) * Nx(1,1)/( Cx(1,1) * Nx(cent_pt,cent_pt)))]);
disp(['Max food production ratio: ',num2str(Fx(1,1)/Fx(cent_pt,cent_pt))]);

rel_Kx = Kx / mean(mean(Kx));
rel_Hx = Hx / mean(mean(Hx));

figure;
scatter(Nx(:),rel_Kx(:),30,'filled','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]); hold on;
scatter(rel_Hx(:),rel_Kx(:),30,'filled','MarkerEdgeColor',[165/255, 165/255, 165/255],'MarkerFaceColor',[165/255, 165/255, 165/255]); 
lin_K1 = polyfit(Nx(:),rel_Kx(:),1);
lin_K1 = polyval(lin_K1,Nx(:));
lin_K2 = polyfit(rel_Hx(:),rel_Kx(:),1);
lin_K2 = polyval(lin_K2,rel_Hx(:));
plot([min(Nx(:)) max(Nx(:))],[min(lin_K1) max(lin_K1)],':','color',[0, 0.4470, 0.7410]);
plot([min(rel_Hx(:)) max(rel_Hx(:))],[min(lin_K2) max(lin_K2)],':','color',[165/255, 165/255, 165/255]);
xlabel('Relative population'); ylabel('Relative capital')
set(gca, 'XTick', 0:5)
set(gca, 'YTick', 0:6)
legend({'Total pop','Labour force'},'location','southoutside','orientation','horizontal')
axis tight
grid on

disp(['d capital / d labour force: ',num2str((min(lin_K2)-max(lin_K2))/(min(rel_Hx(:))-max(rel_Hx(:))))])
disp(['d capital / d population: ',num2str((min(lin_K1)-max(lin_K1))/(min(rel_Hx(:))-max(rel_Hx(:))))])

% update metrics in excel file
global scale
paramvals = [param.thetaN param.Omega param.lambda param.Phi2 param.PhiL param.zeta scale];
varvals = [ Nx(cent_pt,cent_pt)/Nx(1,1) Kx(cent_pt,cent_pt)/Kx(1,1) Hx(cent_pt,cent_pt)*Nx(1,1)/(Hx(1,1)*Nx(cent_pt,cent_pt)) Ex(cent_pt,cent_pt)*Nx(1,1)/(Ex(1,1)*Nx(cent_pt,cent_pt)) Cx(cent_pt,cent_pt)*Nx(1,1)/(Cx(1,1)*Nx(cent_pt,cent_pt)) Fx(1,1)/Fx(cent_pt,cent_pt) (min(lin_K1)-max(lin_K1))/(min(rel_Hx(:))-max(rel_Hx(:)))];
existing = xlsread('steady_state_numbers.xlsx');
new_data = [existing(:,1:14) ; paramvals varvals];
%dlmwrite('steady_state_numbers.xlsx',[paramvals varvals],'delimiter',',','-append')
xlswrite('steady_state_numbers.xlsx',new_data,1,'A2')
                
end

%% Solve model
if logical(run_singlesimulation)

dynare DynamicSpatialModel

% % change for alternative grid sizes
% disp('*-- City vs rural statistics:  --------*')
% disp(['Average population ratio: ',num2str(mean(exp(log_N_4_4))/mean(exp(log_N_1_1)))]);
% disp(['Average capital ratio: ',num2str(mean(exp(log_K_4_4))/mean(exp(log_K_1_1)))]);
% disp(['Average hours per head ratio: ',num2str(mean(exp(log_H_4_4)) * mean(exp(log_N_1_1))/( mean(exp(log_H_1_1)) * mean(exp(log_N_4_4))))]);
% disp(['Average eating per head ratio: ',num2str(mean(exp(log_E_4_4)) * mean(exp(log_N_1_1))/( mean(exp(log_E_1_1)) * mean(exp(log_N_4_4))))]);
% disp(['Average consumption per head ratio: ',num2str(mean(exp(log_C_4_4)) * mean(exp(log_N_1_1))/( mean(exp(log_C_1_1)) * mean(exp(log_N_4_4))))]);
% disp(['Average food production ratio: ',num2str(mean(exp(log_F_1_1))/mean(exp(log_F_4_4)))]);

%save(['../Results/model_2_sim_thetaN',num2str(param_thetaN),'_PhiL',num2str(param_PhiL),'_Phi2',num2str(param_Phi2),'_Omega',num2str(param_Omega),'_zeta',num2str(param_zeta),'_lambda',num2str(param_lambda),'_AT0.0125.mat'])
save(['../Results/model_2_irf_thetaN',num2str(param_thetaN),'_PhiL',num2str(param_PhiL),'_Phi2',num2str(param_Phi2),'_Omega',num2str(param_Omega),'_zeta',num2str(param_zeta),'_lambda',num2str(param_lambda),'_AT0.047.mat'])
delete('param_vals.mat')
end

%% Indeterminacy
if logical(run_indet)
    
% Choose 2 parameters
parameter1_string=['\theta_y'];
parameter2_string=['\theta_{\pi}'];
param1_min = 0.01;
param1_step = 0.5;
param1_max = 1.5;
param2_min = 0.01;
param2_step = 0.5;
param2_max = 2;
param1_name = 'param_theta_pi';
param2_name = 'param_theta_y';

%
ind = zeros(length(param1_vec),length(param2_vec));
save('indet.mat','ind')
param1_vec = [param1_min : param1_step : param1_max];
param2_vec = [param2_min : param2_step : param2_max];

for jj=1:length(param1_vec)
    for ii=1:length(param2_vec)

    eval([param1_name,'= param1_vec(jj);']);
    eval([param2_name,'= param2_vec(ii);']);
    save('param_vals.mat','-append',param1_name,param2_name)

    dynare DynamicSpatialModel 
 
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
    
    save('indet.mat','ind')
    
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
xlabel(parameter2_string,'fontsize',9);
ylabel(parameter1_string,'fontsize',9);    
legend('Determinacy','Indeterminacy','Instability',-1);


end
