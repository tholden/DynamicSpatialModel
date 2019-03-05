
Weight = getWeightMatrix( 'T' , opts.SpatialPointsPerDimension );
Weight = reshape(Weight,[opts.SpatialPoints,1]);
distance = getDistanceMatrix( 'T', opts.SpatialPointsPerDimension ); 
% 
Weight_sim = repmat(Weight,1,opts.simulation_length);
% 
empty_location_variable = zeros( opts.SpatialPoints , opts.simulation_length);
% 
Nx = empty_location_variable;
Kx = empty_location_variable;
Ix = empty_location_variable;
index = 1;
for ii=1:opts.SpatialPointsPerDimension
    for jj=1:opts.SpatialPointsPerDimension
        
        
        eval(['results.model{',num2str(models),'}.',deblank(M_.endo_names(iter,:)),' = oo_.irfs.',deblank(M_.endo_names(iter,:)),'_',char(opts.shock),';']);
        eval(['results.SS{',num2str(models),'}.',deblank(M_.endo_names(iter,:)),' = ',num2str(oo_.steady_state(iter,:)),';']);
        
        results.SS{models}.log_Nx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_N_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Kx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_K_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Ix(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_I_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Cx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_I_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Hx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_H_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Ex(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_E_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Fx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_F_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Qx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_Q_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Jx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_J_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Lx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_L_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_SNx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_SN_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_SDx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_SD_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_muNx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_muN_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Ux(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_U_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_Utildex(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_Utilde_',num2str(ii),'_',num2str(jj)]),:);
        results.SS{models}.log_ATx(index) = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_AT_',num2str(ii),'_',num2str(jj)]),:);
        
        eval(['results.model{models}.log_Nx(index,:) = log_N_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Nx(index);']);
        eval(['results.model{models}.log_Kx(index,:) = log_K_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Kx(index);']);
        eval(['results.model{models}.log_Ix(index,:) = log_I_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Ix(index);']);
        eval(['results.model{models}.log_Cx(index,:) = log_C_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Cx(index);']);
        eval(['results.model{models}.log_Hx(index,:) = log_H_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Hx(index);']);
        eval(['results.model{models}.log_Ex(index,:) = log_E_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Ex(index);']);
        eval(['results.model{models}.log_Fx(index,:) = log_F_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Fx(index);']);
        eval(['results.model{models}.log_Qx(index,:) = log_Q_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Qx(index);']);
        eval(['results.model{models}.log_Jx(index,:) = log_J_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Jx(index);']);
        eval(['results.model{models}.log_Lx(index,:) = log_L_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Lx(index);']);
        eval(['results.model{models}.log_SNx(index,:) = log_SN_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_SNx(index);']);
        eval(['results.model{models}.log_SDx(index,:) = log_SD_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_SDx(index);']);
        eval(['results.model{models}.log_muNx(index,:) = log_muN_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_muNx(index);']);
        eval(['results.model{models}.log_Ux(index,:) = log_U_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Ux(index);']);
        eval(['results.model{models}.log_Utildex(index,:) = log_Utilde_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_Utildex(index);']);
        eval(['results.model{models}.log_ATx(index,:) = log_AT_',num2str(ii),'_',num2str(jj),'_',char(opts.shock),' + results.SS{models}.log_ATx(index);']);
        
        index=index+1;
    end
end
% 
SS.log_tau = oo_.steady_state(strcmp(cellstr(M_.endo_names),'log_tau'),:);
eval(['log_tau = log_tau_',char(opts.shock),' + results.SS{models}.log_tau;']);

for ii=1:length(M_.param_names)
    eval([deblank(M_.param_names(ii,:)),' = ',num2str(M_.params(ii)),';']);
end

Kx_LAG = [ exp(results.SS{models}.log_Kx') exp(results.model{models}.log_Kx(:,1:end-1)) ];

AP = 1;
Ax = AP .* exp(results.model{models}.log_ATx);
ZFx = ( exp(results.model{models}.log_Fx) ./ exp(results.model{models}.log_Lx) .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) );       
SPx = ( 1 - gamma ) .* exp(results.model{models}.log_Fx) ./ ZFx;
Px = empty_location_variable;
for ii=1:opts.SpatialPoints
    Px(ii,:) = ( 1 + lambda ) .* ( sum ( Weight_sim .* exp(results.model{models}.log_Jx) .* ( SPx .* exp( exp(repmat(results.model{models}.log_tau,opts.SpatialPoints,1)) .* repmat( distance(ii,:)',1,opts.simulation_length) ) ) .^ ( - 1 / lambda ) ) ) .^ ( - lambda );
end
Zx = ( ( Kx_LAG .^ alpha .* ( Ax .* exp(results.model{models}.log_Hx) ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* SPx ./ Px ) .^ kappa ) .^ ( 1 ./ ( 1 - kappa ) );
Mx = kappa .* SPx .* Zx ./ Px;
Yx = exp(results.model{models}.log_Cx) + exp(results.model{models}.log_Ix) + Mx;
results.model{models}.Wx = ( 1 - kappa ) .* ( 1 - alpha ) .* SPx .* Zx ./ exp(results.model{models}.log_Hx);
results.model{models}.log_Wx = log(results.model{models}.Wx);
results.model{models}.log_W  = mean(results.model{models}.log_Wx);
results.model{models}.relWx = results.model{models}.log_Wx - results.model{models}.log_W;

results.SS{models}.ZFx = ( exp(results.SS{models}.log_Fx) ./ exp(results.SS{models}.log_Lx) .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) );       
results.SS{models}.SPx = ( 1 - gamma ) .* exp(results.SS{models}.log_Fx) ./ results.SS{models}.ZFx;
for ii=1:opts.SpatialPoints
    results.SS{models}.Px(:,ii) = ( 1 + lambda ) .* ( sum ( Weight' .* exp(results.SS{models}.log_Jx) .* ( results.SS{models}.SPx .* exp( exp(repmat(results.SS{models}.log_tau,opts.SpatialPoints,1))' .* distance(ii,:)) ) .^ ( - 1 / lambda ) ) ) .^ ( - lambda );
end
results.SS{models}.Zx = ( ( exp(results.SS{models}.log_Kx) .^ alpha .* ( exp(results.SS{models}.log_Hx) ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* results.SS{models}.SPx ./ results.SS{models}.Px ) .^ kappa ) .^ ( 1 ./ ( 1 - kappa ) );
results.SS{models}.Wx = ( 1 - kappa ) .* ( 1 - alpha ) .* results.SS{models}.SPx .* results.SS{models}.Zx ./ exp(results.SS{models}.log_Hx);
results.SS{models}.log_Wx = log(results.SS{models}.Wx);
results.SS{models}.log_W  = mean(results.SS{models}.log_Wx);
results.SS{models}.relWx = results.SS{models}.log_Wx - results.SS{models}.log_W;

results.model{models}.Wgrowthx = results.model{models}.Wx ./ [ results.SS{models}.Wx' results.model{models}.Wx(:,1:end-1) ];
results.SS{models}.Wgrowthx = ones(size(results.SS{models}.Wx));

Y = sum( Weight_sim .* Yx  );
