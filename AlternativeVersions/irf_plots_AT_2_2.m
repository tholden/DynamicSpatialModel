clear; close all;

%% File 1
load('Results/model_1.mat')

numpar = M_.param_nbr;
for iter = 1:numpar
    eval([deblank(M_.param_names(iter,:)),' = M_.params(', num2str(iter) ,');']);
end

numvar = M_.endo_nbr;
for iter = 1:numvar
    eval(['model1.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_epsilon_AT_2_2;']);
    eval(['SS1.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
end

model1.E_by_N = model1.log_F;
model1.C_by_N =  model1.log_C;
model1.res_by_N = log( 1 - exp(model1.log_L) );
model1.lagN_2_2 = [ SS1.log_N_2_2 model1.log_N_2_2(1:end-1)  ]  ;
model1.lagN_1_1 = [ SS1.log_N_1_1 model1.log_N_1_1(1:end-1)  ]  ;
model1.E_by_N_2_2 = model1.log_E_2_2 - model1.lagN_2_2;
model1.C_by_N_2_2 =  model1.log_C_2_2 - model1.lagN_2_2;
SS1.res_by_N_2_2 = log( 1 - exp(SS1.log_L_2_2 ) ) - SS1.log_N_2_2;
model1.res_by_N_2_2 = log( 1 - exp(SS1.log_L_2_2 + model1.log_L_2_2) ) - (SS1.log_N_2_2 + model1.lagN_2_2) - SS1.res_by_N_2_2;
model1.E_by_N_1_1 = model1.log_E_1_1 - model1.lagN_1_1;
model1.C_by_N_1_1 =  model1.log_C_1_1 - model1.lagN_1_1;
SS1.res_by_N_1_1 = log( 1 - exp(SS1.log_L_1_1 ) ) - SS1.log_N_1_1;
model1.res_by_N_1_1 =  log( 1 - exp( SS1.log_L_1_1 + model1.log_L_1_1 ) ) - ( SS1.log_N_1_1 + model1.lagN_1_1 ) - SS1.res_by_N_1_1;

F_2_2 = exp(model1.log_F_2_2 + SS1.log_F_2_2);
L_2_2 = exp(model1.log_L_2_2 + SS1.log_L_2_2);
F_1_1 = exp(model1.log_F_1_1 + SS1.log_F_1_1);
L_1_1 = exp(model1.log_L_1_1 + SS1.log_L_1_1);
ZF_2_2 = ( F_2_2 ./ L_2_2.^gamma ).^ ( 1 / ( 1 - gamma ) );
ZF_1_1 = ( F_1_1 ./ L_1_1.^gamma ).^ ( 1 / ( 1 - gamma ) );

SS1.ZF_2_2 = ( exp(SS1.log_F_2_2) ./ exp(SS1.log_L_2_2).^gamma ).^ ( 1 / ( 1 - gamma ) );
SS1.ZF_1_1 = ( exp(SS1.log_F_1_1) ./ exp(SS1.log_L_1_1).^gamma ).^ ( 1 / ( 1 - gamma ) );
model1.log_ZF_2_2 = log(ZF_2_2) - log(SS1.ZF_2_2);
model1.log_ZF_1_1 = log(ZF_1_1) - log(SS1.ZF_1_1);

%% File 2
load('Results/model_3_pic.mat')

numpar = M_.param_nbr;
for iter = 1:numpar
    eval([deblank(M_.param_names(iter,:)),' = M_.params(', num2str(iter) ,');']);
end

numvar = M_.endo_nbr;
for iter = 1:numvar
    eval(['model2.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_epsilon_AT_2_2;']);
    eval(['SS2.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
end

model2.E_by_N = model2.log_F;
model2.C_by_N =  model2.log_C;
model2.res_by_N = log( 1 - exp(model2.log_L) );
model2.lagN_2_2 = [ SS2.log_N_2_2 model2.log_N_2_2(1:end-1)  ]  ;
model2.lagN_1_1 = [ SS2.log_N_1_1 model2.log_N_1_1(1:end-1)  ]  ;
model2.E_by_N_2_2 = model2.log_E_2_2 - model2.lagN_2_2;
model2.C_by_N_2_2 =  model2.log_C_2_2 - model2.lagN_2_2;
SS2.res_by_N_2_2 = log( 1 - exp(SS2.log_L_2_2 ) ) - SS2.log_N_2_2;
model2.res_by_N_2_2 = log( 1 - exp( SS2.log_L_2_2 + model2.log_L_2_2) ) - ( SS2.log_N_2_2 + model2.lagN_2_2 )  - SS2.res_by_N_2_2;
model2.E_by_N_1_1 = model2.log_E_1_1 - model2.lagN_1_1;
model2.C_by_N_1_1 =  model2.log_C_1_1 - model2.lagN_1_1;
SS2.res_by_N_1_1 = log( 1 - exp(SS2.log_L_1_1 ) ) - SS2.log_N_1_1;
model2.res_by_N_1_1 =  log( 1 - exp( SS2.log_L_1_1 + model2.log_L_1_1 ) ) - ( SS2.log_N_1_1 + model2.lagN_1_1 ) - SS2.res_by_N_1_1;

F_2_2 = exp(model2.log_F_2_2 + SS2.log_F_2_2);
L_2_2 = exp(model2.log_L_2_2 + SS2.log_L_2_2);
F_1_1 = exp(model2.log_F_1_1 + SS2.log_F_1_1);
L_1_1 = exp(model2.log_L_1_1 + SS2.log_L_1_1);

model2.ZF_2_2 = ( F_2_2 ./ L_2_2.^gamma ).^ ( 1 / ( 1 - gamma ) );
model2.ZF_1_1 = ( F_1_1 ./ L_1_1.^gamma ).^ ( 1 / ( 1 - gamma ) );
SS2.ZF_2_2 = ( exp(SS2.log_F_2_2) ./ exp(SS2.log_L_2_2).^gamma ).^ ( 1 / ( 1 - gamma ) );
SS2.ZF_1_1 = ( exp(SS2.log_F_1_1) ./ exp(SS2.log_L_1_1).^gamma ).^ ( 1 / ( 1 - gamma ) );
model2.log_ZF_2_2 = log(model2.ZF_2_2) - log(SS2.ZF_2_2);
model2.log_ZF_1_1 = log(model2.ZF_1_1) - log(SS2.ZF_1_1);
%% File 3
load('Results/model_3.mat')

numpar = M_.param_nbr;
for iter = 1:numpar
    eval([deblank(M_.param_names(iter,:)),' = M_.params(', num2str(iter) ,');']);
end

numvar = M_.endo_nbr;
for iter = 1:numvar
    eval(['model3.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_epsilon_AT_2_2;']);
    eval(['SS2.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
end

model3.E_by_N = model3.log_F;
model3.C_by_N =  model3.log_C;
model3.res_by_N = log( 1 - exp(model3.log_L) );
model3.lagN_2_2 = [ SS2.log_N_2_2 model3.log_N_2_2(1:end-1)  ]  ;
model3.lagN_1_1 = [ SS2.log_N_1_1 model3.log_N_1_1(1:end-1)  ]  ;
model3.E_by_N_2_2 = model3.log_E_2_2 - model3.lagN_2_2;
model3.C_by_N_2_2 =  model3.log_C_2_2 - model3.lagN_2_2;
SS2.res_by_N_2_2 = log( 1 - exp(SS2.log_L_2_2 ) ) - SS2.log_N_2_2;
model3.res_by_N_2_2 = log( 1 - exp( SS2.log_L_2_2 + model3.log_L_2_2) ) - ( SS2.log_N_2_2 + model3.lagN_2_2 )  - SS2.res_by_N_2_2;
model3.E_by_N_1_1 = model3.log_E_1_1 - model3.lagN_1_1;
model3.C_by_N_1_1 =  model3.log_C_1_1 - model3.lagN_1_1;
SS2.res_by_N_1_1 = log( 1 - exp(SS2.log_L_1_1 ) ) - SS2.log_N_1_1;
model3.res_by_N_1_1 =  log( 1 - exp( SS2.log_L_1_1 + model3.log_L_1_1 ) ) - ( SS2.log_N_1_1 + model3.lagN_1_1 ) - SS2.res_by_N_1_1;

F_2_2 = exp(model3.log_F_2_2 + SS2.log_F_2_2);
L_2_2 = exp(model3.log_L_2_2 + SS2.log_L_2_2);
F_1_1 = exp(model3.log_F_1_1 + SS2.log_F_1_1);
L_1_1 = exp(model3.log_L_1_1 + SS2.log_L_1_1);

model3.ZF_2_2 = ( F_2_2 ./ L_2_2.^gamma ).^ ( 1 / ( 1 - gamma ) );
model3.ZF_1_1 = ( F_1_1 ./ L_1_1.^gamma ).^ ( 1 / ( 1 - gamma ) );
SS2.ZF_2_2 = ( exp(SS2.log_F_2_2) ./ exp(SS2.log_L_2_2).^gamma ).^ ( 1 / ( 1 - gamma ) );
SS2.ZF_1_1 = ( exp(SS2.log_F_1_1) ./ exp(SS2.log_L_1_1).^gamma ).^ ( 1 / ( 1 - gamma ) );
model3.log_ZF_2_2 = log(model3.ZF_2_2) - log(SS2.ZF_2_2);
model3.log_ZF_1_1 = log(model3.ZF_1_1) - log(SS2.ZF_1_1);
%% Plot results

figure;
subplot(3,3,1); 
    plot(model1.log_C); hold on;
    plot(model2.log_C);
    plot(model3.log_C);
    legend('Original','Land adjustment costs','NK')
    title('log C')
subplot(3,3,2); 
    plot(model1.log_C_2_2); hold on;
    plot(model2.log_C_2_2);
    plot(model3.log_C_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log C_2_2')
subplot(3,3,3); 
    plot(model1.log_C_1_1); hold on;
    plot(model2.log_C_1_1);
    plot(model3.log_C_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log C_1_1')
subplot(3,3,4); 
    plot(model1.C_by_N); hold on;
    plot(model2.C_by_N);
    plot(model3.C_by_N);
    legend('Original','Land adjustment costs','NK')
    title('log C_by_N')
subplot(3,3,5); 
    plot(model1.C_by_N_2_2); hold on;
    plot(model2.C_by_N_2_2);
    plot(model3.C_by_N_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log C_by_N_2_2')
subplot(3,3,6); 
    plot(model1.C_by_N_1_1); hold on;
    plot(model2.C_by_N_1_1);
    plot(model3.C_by_N_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log C_by_N_1_1')
subplot(3,3,7); 
    plot(model1.log_I); hold on;
    plot(model2.log_I);
    plot(model3.log_I);
    legend('Original','Land adjustment costs','NK')
    title('log I')
subplot(3,3,8); 
    plot(model1.log_I_2_2); hold on;
    plot(model2.log_I_2_2);
    plot(model3.log_I_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log I_2_2')
subplot(3,3,9); 
    plot(model1.log_I_1_1); hold on;
    plot(model2.log_I_1_1);
    plot(model3.log_I_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log I_1_1')
    
figure;
subplot(3,3,1); 
    plot(model1.log_H); hold on;
    plot(model2.log_H);
    plot(model3.log_H);
    legend('Original','Land adjustment costs','NK')
    title('log H')
subplot(3,3,2); 
    plot(model1.log_H_2_2); hold on;
    plot(model2.log_H_2_2);
    plot(model3.log_H_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log H_2_2')
subplot(3,3,3); 
    plot(model1.log_H_1_1); hold on;
    plot(model2.log_H_1_1);
    plot(model3.log_H_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log H_1_1')
subplot(3,3,4); 
    plot(model1.log_F); hold on;
    plot(model2.log_F);
    plot(model3.log_F);
    legend('Original','Land adjustment costs','NK')
    title('log E')
subplot(3,3,5); 
    plot(model1.log_E_2_2); hold on;
    plot(model2.log_E_2_2);
    plot(model3.log_E_2_2);
    legend('Original','Land adjustment costs')
    title('log E_2_2')
subplot(3,3,6); 
    plot(model1.log_E_1_1); hold on;
    plot(model2.log_E_1_1);
    plot(model3.log_E_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log E_1_1')
subplot(3,3,7); 
    plot(model1.E_by_N); hold on;
    plot(model2.E_by_N);
    plot(model3.E_by_N);
    legend('Original','Land adjustment costs','NK')
    title('log E_by_N')
subplot(3,3,8);
    plot(model1.E_by_N_2_2); hold on;
    plot(model2.E_by_N_2_2);
    plot(model3.E_by_N_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log E_by_N_2_2')
subplot(3,3,9); 
    plot(model1.E_by_N_1_1); hold on;
    plot(model2.E_by_N_1_1);
    plot(model3.E_by_N_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log E_by_N_1_1')
    
    
figure;
subplot(3,3,1); 
    plot(model1.log_L); hold on;
    plot(model2.log_L);
    plot(model3.log_L);
    legend('Original','Land adjustment costs','NK')
    title('log L')
subplot(3,3,2); 
    plot(model1.log_L_2_2); hold on;
    plot(model2.log_L_2_2);
    plot(model3.log_L_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log L_2_2')
subplot(3,3,3); 
    plot(model1.log_L_1_1); hold on;
    plot(model2.log_L_1_1);
    plot(model3.log_L_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log L_1_1')
subplot(3,3,4); 
    plot(model1.log_SN); hold on;
    plot(model2.log_SN);
    plot(model3.log_SN);
    legend('Original','Land adjustment costs','NK')
    title('log SN')
subplot(3,3,5); 
    plot(model1.log_N_2_2); hold on;
    plot(model2.log_N_2_2);
    plot(model3.log_N_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log N_2_2')
subplot(3,3,6); 
    plot(model1.log_N_1_1); hold on;
    plot(model2.log_N_1_1);
    plot(model3.log_N_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log N_1_1')
subplot(3,3,7); 
    plot(model1.res_by_N); hold on;
    plot(model2.res_by_N);
    plot(model3.res_by_N);
    legend('Original','Land adjustment costs','NK')
    title('log res_by_N')
subplot(3,3,8); 
    plot(model1.res_by_N_2_2); hold on;
    plot(model2.res_by_N_2_2);
    plot(model3.res_by_N_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log res_by_N_2_2')
subplot(3,3,9); 
    plot(model1.res_by_N_1_1); hold on;
    plot(model2.res_by_N_1_1);
    plot(model3.res_by_N_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log res_by_N_1_1')

    
figure;
subplot(3,3,1); 
    plot(model1.log_F); hold on;
    plot(model2.log_F);
    plot(model3.log_F);
    legend('Original','Land adjustment costs','NK')
    title('log F')
subplot(3,3,2); 
    plot(model1.log_F_2_2); hold on;
    plot(model2.log_F_2_2);
    plot(model3.log_F_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log F_2_2')
subplot(3,3,3); 
    plot(model1.log_F_1_1); hold on;
    plot(model2.log_F_1_1);
    plot(model3.log_F_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log F_1_1')
subplot(3,3,4); 
    plot(model1.log_U); hold on;
    plot(model2.log_U);
    plot(model3.log_U);
    legend('Original','Land adjustment costs','NK')
    title('log U')
subplot(3,3,5); 
    plot(model1.log_U_2_2); hold on;
    plot(model2.log_U_2_2);
    plot(model3.log_U_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log U_2_2')
subplot(3,3,6); 
    plot(model1.log_U_1_1); hold on;
    plot(model2.log_U_1_1);
    plot(model3.log_U_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log U_1_1')
subplot(3,3,7); 
    plot(model1.log_ZF_2_2); hold on;
    plot(model2.log_ZF_2_2);
    plot(model3.log_ZF_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log ZF_2_2')
subplot(3,3,8); 
    plot(model1.log_ZF_1_1); hold on;
    plot(model2.log_ZF_1_1);
    plot(model3.log_ZF_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log ZF_1_1')
    
 
figure;
subplot(3,3,1); 
    plot(model1.log_J); hold on;
    plot(model2.log_J);
    plot(model3.log_J);
    legend('Original','Land adjustment costs','NK')
    title('log J')
subplot(3,3,2); 
    plot(model1.log_J_2_2); hold on;
    plot(model2.log_J_2_2);
    plot(model3.log_J_2_2);
    legend('Original','Land adjustment costs','NK')
    title('log J_2_2')
subplot(3,3,3); 
    plot(model1.log_J_1_1); hold on;
    plot(model2.log_J_1_1);
    plot(model3.log_J_1_1);
    legend('Original','Land adjustment costs','NK')
    title('log J_1_1')   
    