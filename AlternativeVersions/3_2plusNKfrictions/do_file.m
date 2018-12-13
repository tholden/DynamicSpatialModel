clear; %close all;
%dbstop if error

run_indet = 1;

param_thetaN = 10;
param_PhiL = 4;
param_xi = .8;
param_rho_R = 8;
param_theta_y = 0.2;
param_theta_pi = 1.9;
save('param_vals.mat')

%% Solve model
if ~logical(run_indet)

dynare DynamicSpatialModel
save('../Results/model_3_pic.mat')

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
legend('Determinacy','Indeterminacy','Instability',-1);


end

%%
clear;
