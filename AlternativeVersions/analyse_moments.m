clear; close all;

opts.results_files = 'model_2_sim_thetaN0.7994_PhiL2_Phi24_Omega1.7_zeta8_lambda0.05';
opts.SpatialPointsPerDimension = 7;
%%
load(['Results/',opts.results_files,'.mat'],'oo_','M_')

for ii=1:opts.SpatialPointsPerDimension
    for jj=1:opts.SpatialPointsPerDimension
        N{ii,jj} = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_N_',num2str(ii),'_',num2str(jj)]),:);
        K{ii,jj} = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_K_',num2str(ii),'_',num2str(jj)]),:);
        Y{ii,jj} = oo_.endo_simul(strcmp(cellstr(M_.endo_names),['log_Y_',num2str(ii),'_',num2str(jj)]),:);
    end
end

for t=1:10000
    Nt(:,:,t) = [ ...
         N{1,1}(t) N{1,2}(t) N{1,3}(t) N{1,4}(t) N{1,5}(t) N{1,6}(t)  N{1,7}(t) ; ...
         N{2,1}(t) N{2,2}(t) N{2,3}(t) N{2,4}(t) N{2,5}(t) N{2,6}(t)  N{2,7}(t) ; ...
         N{3,1}(t) N{3,2}(t) N{3,3}(t) N{3,4}(t) N{3,5}(t) N{3,6}(t)  N{3,7}(t) ; ...
         N{4,1}(t) N{4,2}(t) N{4,3}(t) N{4,4}(t) N{4,5}(t) N{4,6}(t)  N{4,7}(t) ; ...
         N{5,1}(t) N{5,2}(t) N{5,3}(t) N{5,4}(t) N{5,5}(t) N{5,6}(t)  N{5,7}(t) ; ...
         N{6,1}(t) N{6,2}(t) N{6,3}(t) N{6,4}(t) N{6,5}(t) N{6,6}(t)  N{6,7}(t) ; ...
         N{7,1}(t) N{7,2}(t) N{7,3}(t) N{7,4}(t) N{7,5}(t) N{7,6}(t)  N{7,7}(t) ];
end

figure;
subplot(3,3,1); contourf(exp(Nt(:,:,1000))); colorbar; title('N (t=1000)');
subplot(3,3,2); contourf(exp(Nt(:,:,2000))); colorbar; title('N (t=2000)');
subplot(3,3,3); contourf(exp(Nt(:,:,3000))); colorbar; title('N (t=3000)');
subplot(3,3,4); contourf(exp(Nt(:,:,4000))); colorbar; title('N (t=4000)');
subplot(3,3,5); contourf(exp(Nt(:,:,5000))); colorbar; title('N (t=5000)');
subplot(3,3,6); contourf(exp(Nt(:,:,6000))); colorbar; title('N (t=6000)');
subplot(3,3,7); contourf(exp(Nt(:,:,7000))); colorbar; title('N (t=7000)');
subplot(3,3,8); contourf(exp(Nt(:,:,8000))); colorbar; title('N (t=8000)');
subplot(3,3,9); contourf(exp(Nt(:,:,9000))); colorbar; title('N (t=9000)');

figure;
subplot(3,3,1); plot(N{1,1}(4000:4400)); title('log N{1,1}')
subplot(3,3,2); plot(N{1,4}(4000:4400)); title('log N{1,4}')
subplot(3,3,3); plot(N{1,7}(4000:4400)); title('log N{1,7}')
subplot(3,3,4); plot(N{4,1}(4000:4400)); title('log N{4,1}')
subplot(3,3,5); plot(N{4,4}(4000:4400)); title('log N{4,4}')
subplot(3,3,6); plot(N{4,7}(4000:4400)); title('log N{4,7}')
subplot(3,3,7); plot(N{7,1}(4000:4400)); title('log N{7,1}')
subplot(3,3,8); plot(N{7,4}(4000:4400)); title('log N{7,4}')
subplot(3,3,9); plot(N{7,7}(4000:4400)); title('log N{7,7}')

figure;
subplot(3,3,1); plot(K{1,1}(4000:4400)); title('log K{1,1}')
subplot(3,3,2); plot(K{1,4}(4000:4400)); title('log K{1,4}')
subplot(3,3,3); plot(K{1,7}(4000:4400)); title('log K{1,7}')
subplot(3,3,4); plot(K{4,1}(4000:4400)); title('log K{4,1}')
subplot(3,3,5); plot(K{4,4}(4000:4400)); title('log K{4,4}')
subplot(3,3,6); plot(K{4,7}(4000:4400)); title('log K{4,7}')
subplot(3,3,7); plot(K{7,1}(4000:4400)); title('log K{7,1}')
subplot(3,3,8); plot(K{7,4}(4000:4400)); title('log K{7,4}')
subplot(3,3,9); plot(K{7,7}(4000:4400)); title('log K{7,7}')

figure;
subplot(3,3,1); plot(Y{1,1}(4000:4400)); title('log Y{1,1}')
subplot(3,3,2); plot(Y{1,4}(4000:4400)); title('log Y{1,4}')
subplot(3,3,3); plot(Y{1,7}(4000:4400)); title('log Y{1,7}')
subplot(3,3,4); plot(Y{4,1}(4000:4400)); title('log Y{4,1}')
subplot(3,3,5); plot(Y{4,4}(4000:4400)); title('log Y{4,4}')
subplot(3,3,6); plot(Y{4,7}(4000:4400)); title('log Y{4,7}')
subplot(3,3,7); plot(Y{7,1}(4000:4400)); title('log Y{7,1}')
subplot(3,3,8); plot(Y{7,4}(4000:4400)); title('log Y{7,4}')
subplot(3,3,9); plot(Y{7,7}(4000:4400)); title('log Y{7,7}')