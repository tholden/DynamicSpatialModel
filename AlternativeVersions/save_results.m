function [model_results,SS_results,init_results] = save_results(M_,oo_,shock,varargin)
%SAVE_RESULTS Stores impulse response functions in desired format for plotting  
%   save_results(M_,oo_,IRF_matrix,init_values)
%   M_ and oo_ are the standard structs outputted by dynare
%   shock is the choice of IRF shock, e.g. epsilon_MP or epsilon_AT_2_2
%   IRF_matrix is an optional NxT matrix where N is number of variables
%   and T is IRF periods, in the same order as M_.endo_names
%   init_values is optional starting values which can be different from
%   steay state

numpar = M_.param_nbr;
for iter = 1:numpar
    eval(['p.',deblank(M_.param_names(iter,:)),' = M_.params(', num2str(iter) ,');']);
end

numvar = M_.endo_nbr;
if nargin ==3   
    for iter = 1:numvar
        eval(['model_results.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_',shock,';']);
        eval(['SS_results.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
    end
    init_results = SS_results;
elseif nargin == 4
    IRF_matrix = varargin{1,1};
    for iter = 1:numvar
        eval(['model_results.',deblank(M_.endo_names(iter,:)),' = IRF_matrix(',num2str(iter),',:);']);
        eval(['SS_results.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
    end
    init_results = SS_results;
elseif nargin == 5
    IRF_matrix = varargin{1,1};
    init_values = varargin{1,2};
    for iter = 1:numvar
        eval(['model_results.',deblank(M_.endo_names(iter,:)),' = IRF_matrix(',num2str(iter),',:);']);
        eval(['SS_results.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
        eval(['init_results.',deblank(M_.endo_names(iter,:)),' = init_values(',num2str(iter),',:);']);
    end
else 
    disp('Incorrect inputs')
    return
end

model_results.E_by_N = model_results.log_F;
model_results.C_by_N =  model_results.log_C;
model_results.res_by_N = log( 1 - exp(model_results.log_L + SS_results.log_L) ) - log( 1 - exp(SS_results.log_L) );
model_results.lagN_2_2 = [ init_results.log_N_2_2 model_results.log_N_2_2(1:end-1)  ]  ;
model_results.lagN_1_1 = [ init_results.log_N_1_1 model_results.log_N_1_1(1:end-1)  ]  ;
model_results.E_by_N_2_2 = model_results.log_E_2_2 - model_results.lagN_2_2;
model_results.C_by_N_2_2 =  model_results.log_C_2_2 - model_results.lagN_2_2;
SS_results.res_by_N_2_2 = log( 1 - exp(SS_results.log_L_2_2 ) ) - SS_results.log_N_2_2;
model_results.res_by_N_2_2 = log( 1 - exp(SS_results.log_L_2_2 + model_results.log_L_2_2) ) - (SS_results.log_N_2_2 + model_results.lagN_2_2) - SS_results.res_by_N_2_2;
model_results.E_by_N_1_1 = model_results.log_E_1_1 - model_results.lagN_1_1;
model_results.C_by_N_1_1 =  model_results.log_C_1_1 - model_results.lagN_1_1;
SS_results.res_by_N_1_1 = log( 1 - exp(SS_results.log_L_1_1 ) ) - SS_results.log_N_1_1;
model_results.res_by_N_1_1 =  log( 1 - exp( SS_results.log_L_1_1 + model_results.log_L_1_1 ) ) - ( SS_results.log_N_1_1 + model_results.lagN_1_1 ) - SS_results.res_by_N_1_1;

F_2_2 = exp(model_results.log_F_2_2 + SS_results.log_F_2_2);
L_2_2 = exp(model_results.log_L_2_2 + SS_results.log_L_2_2);
F_1_1 = exp(model_results.log_F_1_1 + SS_results.log_F_1_1);
L_1_1 = exp(model_results.log_L_1_1 + SS_results.log_L_1_1);
ZF_2_2 = ( F_2_2 ./ L_2_2.^p.gamma ).^ ( 1 / ( 1 - p.gamma ) );
ZF_1_1 = ( F_1_1 ./ L_1_1.^p.gamma ).^ ( 1 / ( 1 - p.gamma ) );

SS_results.ZF_2_2 = ( exp(SS_results.log_F_2_2) ./ exp(SS_results.log_L_2_2).^p.gamma ).^ ( 1 / ( 1 - p.gamma ) );
SS_results.ZF_1_1 = ( exp(SS_results.log_F_1_1) ./ exp(SS_results.log_L_1_1).^p.gamma ).^ ( 1 / ( 1 - p.gamma ) );
model_results.log_ZF_2_2 = log(ZF_2_2) - log(SS_results.ZF_2_2);
model_results.log_ZF_1_1 = log(ZF_1_1) - log(SS_results.ZF_1_1);

end