%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'steady_state_solver';
M_.dynare_version = '4.5.1';
oo_.dynare_version = '4.5.1';
options_.dynare_version = '4.5.1';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('steady_state_solver.log');
M_.exo_names = 'epsilon_GA';
M_.exo_names_tex = 'epsilon\_GA';
M_.exo_names_long = 'epsilon_GA';
M_.exo_names = char(M_.exo_names, 'epsilon_GN');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_GN');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_GN');
M_.exo_names = char(M_.exo_names, 'epsilon_tau');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_tau');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_tau');
M_.exo_names = char(M_.exo_names, 'epsilon_phi');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_phi');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_phi');
M_.exo_names = char(M_.exo_names, 'epsilon_beta');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_beta');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_beta');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_1_1');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_1\_1');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_1_1');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_1_2');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_1\_2');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_1_2');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_1_3');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_1\_3');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_1_3');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_2_1');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_2\_1');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_2_1');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_2_2');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_2\_2');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_2_2');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_2_3');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_2\_3');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_2_3');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_3_1');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_3\_1');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_3_1');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_3_2');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_3\_2');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_3_2');
M_.exo_names = char(M_.exo_names, 'epsilon_AT_3_3');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsilon\_AT\_3\_3');
M_.exo_names_long = char(M_.exo_names_long, 'epsilon_AT_3_3');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names_long = 'alpha';
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names_long = char(M_.param_names_long, 'kappa');
M_.param_names = char(M_.param_names, 'nu');
M_.param_names_tex = char(M_.param_names_tex, 'nu');
M_.param_names_long = char(M_.param_names_long, 'nu');
M_.param_names = char(M_.param_names, 'varsigma');
M_.param_names_tex = char(M_.param_names_tex, 'varsigma');
M_.param_names_long = char(M_.param_names_long, 'varsigma');
M_.param_names = char(M_.param_names, 'zeta');
M_.param_names_tex = char(M_.param_names_tex, 'zeta');
M_.param_names_long = char(M_.param_names_long, 'zeta');
M_.param_names = char(M_.param_names, 'lambda');
M_.param_names_tex = char(M_.param_names_tex, 'lambda');
M_.param_names_long = char(M_.param_names_long, 'lambda');
M_.param_names = char(M_.param_names, 'deltaJ');
M_.param_names_tex = char(M_.param_names_tex, 'deltaJ');
M_.param_names_long = char(M_.param_names_long, 'deltaJ');
M_.param_names = char(M_.param_names, 'deltaK');
M_.param_names_tex = char(M_.param_names_tex, 'deltaK');
M_.param_names_long = char(M_.param_names_long, 'deltaK');
M_.param_names = char(M_.param_names, 'Phi2');
M_.param_names_tex = char(M_.param_names_tex, 'Phi2');
M_.param_names_long = char(M_.param_names_long, 'Phi2');
M_.param_names = char(M_.param_names, 'thetaC');
M_.param_names_tex = char(M_.param_names_tex, 'thetaC');
M_.param_names_long = char(M_.param_names_long, 'thetaC');
M_.param_names = char(M_.param_names, 'thetaF');
M_.param_names_tex = char(M_.param_names_tex, 'thetaF');
M_.param_names_long = char(M_.param_names_long, 'thetaF');
M_.param_names = char(M_.param_names, 'thetaL');
M_.param_names_tex = char(M_.param_names_tex, 'thetaL');
M_.param_names_long = char(M_.param_names_long, 'thetaL');
M_.param_names = char(M_.param_names, 'thetaH');
M_.param_names_tex = char(M_.param_names_tex, 'thetaH');
M_.param_names_long = char(M_.param_names_long, 'thetaH');
M_.param_names = char(M_.param_names, 'thetaN');
M_.param_names_tex = char(M_.param_names_tex, 'thetaN');
M_.param_names_long = char(M_.param_names_long, 'thetaN');
M_.param_names = char(M_.param_names, 'psi1');
M_.param_names_tex = char(M_.param_names_tex, 'psi1');
M_.param_names_long = char(M_.param_names_long, 'psi1');
M_.param_names = char(M_.param_names, 'psi2');
M_.param_names_tex = char(M_.param_names_tex, 'psi2');
M_.param_names_long = char(M_.param_names_long, 'psi2');
M_.param_names = char(M_.param_names, 'psi3');
M_.param_names_tex = char(M_.param_names_tex, 'psi3');
M_.param_names_long = char(M_.param_names_long, 'psi3');
M_.param_names = char(M_.param_names, 'Gamma');
M_.param_names_tex = char(M_.param_names_tex, 'Gamma');
M_.param_names_long = char(M_.param_names_long, 'Gamma');
M_.param_names = char(M_.param_names, 'Omega');
M_.param_names_tex = char(M_.param_names_tex, 'Omega');
M_.param_names_long = char(M_.param_names_long, 'Omega');
M_.param_names = char(M_.param_names, 'dBar');
M_.param_names_tex = char(M_.param_names_tex, 'dBar');
M_.param_names_long = char(M_.param_names_long, 'dBar');
M_.param_names = char(M_.param_names, 'UtilityParamSum');
M_.param_names_tex = char(M_.param_names_tex, 'UtilityParamSum');
M_.param_names_long = char(M_.param_names_long, 'UtilityParamSum');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 14;
M_.endo_nbr = 0;
M_.param_nbr = 22;
M_.orig_endo_nbr = 0;
M_.aux_vars = [];
M_.Sigma_e = zeros(14, 14);
M_.Correlation_matrix = eye(14, 14);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('steady_state_solver_static');
erase_compiled_function('steady_state_solver_dynamic');
M_.orig_eq_nbr = 0;
M_.eq_nbr = 0;
M_.ramsey_eq_nbr = 0;
clear;
load('param_vals.mat');
M_.params( 1 ) = 0.3;
alpha = M_.params( 1 );
M_.params( 2 ) = 0.5;
gamma = M_.params( 2 );
M_.params( 3 ) = 0.5;
kappa = M_.params( 3 );
M_.params( 4 ) = 2;
nu = M_.params( 4 );
M_.params( 5 ) = 1.5;
varsigma = M_.params( 5 );
M_.params( 6 ) = 8;
zeta = M_.params( 6 );
M_.params( 7 ) = 0.1;
lambda = M_.params( 7 );
M_.params( 8 ) = 0.01;
deltaJ = M_.params( 8 );
M_.params( 9 ) = 0.03;
deltaK = M_.params( 9 );
M_.params( 10 ) = 4;
Phi2 = M_.params( 10 );
M_.params( 11 ) = 4;
thetaC = M_.params( 11 );
M_.params( 12 ) = 1;
thetaF = M_.params( 12 );
M_.params( 13 ) = 0.3333333333333333*M_.params(12)*M_.params(2);
thetaL = M_.params( 13 );
M_.params( 14 ) = 4;
thetaH = M_.params( 14 );
M_.params( 15 ) = param_thetaN;
thetaN = M_.params( 15 );
M_.params( 16 ) = 0.5;
psi1 = M_.params( 16 );
M_.params( 17 ) = 0.5;
psi2 = M_.params( 17 );
M_.params( 18 ) = M_.params(16)*0.02/0.98;
psi3 = M_.params( 18 );
M_.params( 22 ) = M_.params(16)+M_.params(12)+M_.params(11)+M_.params(13)+M_.params(14)+M_.params(15)+M_.params(17)+M_.params(18);
UtilityParamSum = M_.params( 22 );
M_.params( 11 ) = M_.params(11)/M_.params(22);
thetaC = M_.params( 11 );
M_.params( 12 ) = M_.params(12)/M_.params(22);
thetaF = M_.params( 12 );
M_.params( 13 ) = M_.params(13)/M_.params(22);
thetaL = M_.params( 13 );
M_.params( 14 ) = M_.params(14)/M_.params(22);
thetaH = M_.params( 14 );
M_.params( 15 ) = M_.params(15)/M_.params(22);
thetaN = M_.params( 15 );
M_.params( 16 ) = M_.params(16)/M_.params(22);
psi1 = M_.params( 16 );
M_.params( 17 ) = M_.params(17)/M_.params(22);
psi2 = M_.params( 17 );
M_.params( 18 ) = M_.params(18)/M_.params(22);
psi3 = M_.params( 18 );
M_.params( 19 ) = 1;
Gamma = M_.params( 19 );
M_.params( 20 ) = 3;
Omega = M_.params( 20 );
M_.params( 21 ) = 1.414213562373095;
dBar = M_.params( 21 );
save('steady_state_solver_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('steady_state_solver_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('steady_state_solver_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('steady_state_solver_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('steady_state_solver_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('steady_state_solver_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('steady_state_solver_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
disp('Note: 158 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
