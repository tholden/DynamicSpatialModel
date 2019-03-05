
Weight = getWeightMatrix( 'T' , opts.SpatialPointsPerDimension );
Weight = reshape(Weight,[opts.SpatialPoints,1]);
distance = getDistanceMatrix( 'T', opts.SpatialPointsPerDimension ); 

Weight_sim = repmat(Weight,1,opts.simulation_length-1);

empty_location_variable = zeros( opts.SpatialPoints , opts.simulation_length);
empty_location_variable_trunc = zeros( opts.SpatialPoints , opts.simulation_length-1);

Nx = empty_location_variable;
Kx = empty_location_variable;
Ix = empty_location_variable;
index = 1;
for ii=1:opts.SpatialPointsPerDimension
    for jj=1:opts.SpatialPointsPerDimension
        SS.log_Nx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_N_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Kx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_K_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Ix = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_I_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Cx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_I_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Hx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_H_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Ex = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_E_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Fx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_F_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Qx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_Q_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Jx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_J_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Lx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_L_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_SNx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_SN_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_SDx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_SD_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_muNx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_muN_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Ux = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_U_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_Utildex = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_Utilde_',num2str(ii),'_',num2str(jj)]),:);
        SS.log_ATx = oo_.steady_state(strcmp(cellstr(M_.endo_names),['log_AT_',num2str(ii),'_',num2str(jj)]),:);
        
        log_Nx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_N_',num2str(ii),'_',num2str(jj)]),:);
        log_Kx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_K_',num2str(ii),'_',num2str(jj)]),:);
        log_Ix(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_I_',num2str(ii),'_',num2str(jj)]),:);
        log_Cx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_C_',num2str(ii),'_',num2str(jj)]),:);
        log_Hx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_H_',num2str(ii),'_',num2str(jj)]),:);
        log_Ex(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_E_',num2str(ii),'_',num2str(jj)]),:);
        log_Fx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_F_',num2str(ii),'_',num2str(jj)]),:);
        log_Qx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_Q_',num2str(ii),'_',num2str(jj)]),:);
        log_Jx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_J_',num2str(ii),'_',num2str(jj)]),:);
        log_Lx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_L_',num2str(ii),'_',num2str(jj)]),:);
        log_SNx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_SN_',num2str(ii),'_',num2str(jj)]),:);
        log_SDx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_SD_',num2str(ii),'_',num2str(jj)]),:);
        log_muNx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_muN_',num2str(ii),'_',num2str(jj)]),:);
        log_Ux(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_U_',num2str(ii),'_',num2str(jj)]),:);
        log_Utildex(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_Utilde_',num2str(ii),'_',num2str(jj)]),:);
        log_ATx(index,:) = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_AT_',num2str(ii),'_',num2str(jj)]),:);
        
        index=index+1;
    end
end

SS.log_I = oo_.steady_state(strcmp(cellstr(M_.endo_names),'log_I'),:);
SS.log_SN = oo_.steady_state(strcmp(cellstr(M_.endo_names),'log_SN'),:);
SS.log_tau = oo_.steady_state(strcmp(cellstr(M_.endo_names),'log_tau'),:);

log_I = oo_.endo_simul(strcmp(cellstr(M_.endo_names),'log_I'),:);
log_SN = oo_.endo_simul(strcmp(cellstr(M_.endo_names),'log_SN'),:);
log_tau = oo_.endo_simul(strcmp(cellstr(M_.endo_names),'log_tau'),:);


for ii=1:length(M_.param_names)
    eval([deblank(M_.param_names(ii,:)),' = ',num2str(M_.params(ii)),';']);
end

% Lags and truncating
Kx_LAG = exp( log_Kx(:,1:end-1) );
log_Kx = log_Kx(:,2:end);
log_Nx = log_Nx(:,2:end);
log_Ix = log_Ix(:,2:end);
log_Cx = log_Cx(:,2:end);
log_Hx = log_Hx(:,2:end);
log_Ex = log_Ex(:,2:end);
log_Fx = log_Fx(:,2:end);
log_Qx = log_Qx(:,2:end);
log_Jx = log_Jx(:,2:end);
log_Lx = log_Lx(:,2:end);
log_SNx = log_SNx(:,2:end);
log_SDx = log_SDx(:,2:end);
log_muNx = log_muNx(:,2:end);
log_Ux = log_Ux(:,2:end);
log_ATx = log_ATx(:,2:end);
log_Utildex = log_Utildex(:,2:end);

log_I = log_I(2:end);
log_SN = log_SN(2:end);
log_tau = log_tau(2:end);


AP = 1;
Ax = AP .* exp(log_ATx);
ZFx = ( exp(log_Fx) ./ exp(log_Lx) .^ gamma ) .^ ( 1 ./ ( 1 - gamma ) );       

SPx = ( 1 - gamma ) .* exp(log_Fx) ./ ZFx;

Px = empty_location_variable_trunc;
for ii=1:opts.SpatialPoints
    Px(ii,:) = ( 1 + lambda ) .* ( sum ( Weight_sim .* exp(log_Jx) .* ( SPx .* exp( exp(repmat(log_tau,opts.SpatialPoints,1)) .* repmat( distance(ii,:)',1,opts.simulation_length-1) ) ) .^ ( - 1 / lambda ) ) ) .^ ( - lambda );
end
Zx = ( ( Kx_LAG .^ alpha .* ( Ax .* exp(log_Hx) ) .^ ( 1 - alpha ) ) .^ ( 1 - kappa ) .* ( kappa .* SPx ./ Px ) .^ kappa ) .^ ( 1 ./ ( 1 - kappa ) );
Mx = kappa .* SPx .* Zx ./ Px;
Yx = exp(log_Cx) + exp(log_Ix) + Mx;

Wx = ( 1 - kappa ) .* ( 1 - alpha ) .* SPx .* Zx ./ exp(log_Hx);


Y = sum( Weight_sim .* Yx  );


%% Variables in time - 3-D matrices, in form (i,j,t)

% Nt = [];
% for t=1:10000
%     for ii=1:opts.SpatialPointsPerDimension
%         Nt_row = [];
%         for jj=1:opts.SpatialPointsPerDimension
%             Nt_row = [ Nt_row N{ii,jj}(t) ]; 
%         end
%         Nt = [Nt ; Nt_row];
%     end
% end
% 
% for t=1:10000
%     Nt(:,:,t) = [ ...
%          N{1,1}(t) N{1,2}(t) N{1,3}(t) N{1,4}(t) N{1,5}(t) N{1,6}(t)  N{1,7}(t) ; ...
%          N{2,1}(t) N{2,2}(t) N{2,3}(t) N{2,4}(t) N{2,5}(t) N{2,6}(t)  N{2,7}(t) ; ...
%          N{3,1}(t) N{3,2}(t) N{3,3}(t) N{3,4}(t) N{3,5}(t) N{3,6}(t)  N{3,7}(t) ; ...
%          N{4,1}(t) N{4,2}(t) N{4,3}(t) N{4,4}(t) N{4,5}(t) N{4,6}(t)  N{4,7}(t) ; ...
%          N{5,1}(t) N{5,2}(t) N{5,3}(t) N{5,4}(t) N{5,5}(t) N{5,6}(t)  N{5,7}(t) ; ...
%          N{6,1}(t) N{6,2}(t) N{6,3}(t) N{6,4}(t) N{6,5}(t) N{6,6}(t)  N{6,7}(t) ; ...
%          N{7,1}(t) N{7,2}(t) N{7,3}(t) N{7,4}(t) N{7,5}(t) N{7,6}(t)  N{7,7}(t) ];
% end
% 
% figure;
% subplot(3,3,1); contourf(exp(Nt(:,:,1000))); colorbar; title('N (t=1000)');
% subplot(3,3,2); contourf(exp(Nt(:,:,2000))); colorbar; title('N (t=2000)');
% subplot(3,3,3); contourf(exp(Nt(:,:,3000))); colorbar; title('N (t=3000)');
% subplot(3,3,4); contourf(exp(Nt(:,:,4000))); colorbar; title('N (t=4000)');
% subplot(3,3,5); contourf(exp(Nt(:,:,5000))); colorbar; title('N (t=5000)');
% subplot(3,3,6); contourf(exp(Nt(:,:,6000))); colorbar; title('N (t=6000)');
% subplot(3,3,7); contourf(exp(Nt(:,:,7000))); colorbar; title('N (t=7000)');
% subplot(3,3,8); contourf(exp(Nt(:,:,8000))); colorbar; title('N (t=8000)');
% subplot(3,3,9); contourf(exp(Nt(:,:,9000))); colorbar; title('N (t=9000)');
