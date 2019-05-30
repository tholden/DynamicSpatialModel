%clear; 

% Data file
%load('../Results/model_2_thetaN5_PhiL0.mat')

% Options
SpatialPointsPerDimension = 4;


%%
capital = [];
population = [];
for ii=1:SpatialPointsPerDimension
    for jj=1:SpatialPointsPerDimension
        eval( [ 'capital = [ capital ; log_K_',num2str(ii),'_',num2str(jj),'-log_K ];']);
        eval( [ 'population = [ population ; log_N_',num2str(ii),'_',num2str(jj),' ];']);
    end
end

figure;
subplot(1,2,1)
hist(capital); title('Capital')
subplot(1,2,2)
hist(population); title('Population')
