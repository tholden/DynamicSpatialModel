clear; close all;
%dbstop if error

run_indet = 0;
run_single = 1;
global first_run
first_run=1;

%% Solve model
if logical(run_single)
param_thetaN = 10; %10
param_Omega = 3; %3
param_lambda = 0.1; %0.1

param_Phi2 = 4; %4
param_PhiL = 2; %2
param_zeta = 8; %8
save('param_vals.mat')

dynare DynamicSpatialModel
save(['../Results/model_2_thetaN',num2str(param_thetaN),'_PhiL',num2str(param_PhiL),'_Phi2',num2str(param_Phi2),'_Omega',num2str(param_Omega),'_zeta',num2str(param_zeta),'_lambda',num2str(param_lambda),'.mat'])


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
