%% Options

% choose plots:
time_series_for_selection_of_locations = 0;
time_series_for_all_locations = 0;

timestart = 1;
timeend = 10000;

%% time series for selection of locations:
if logical(time_series_for_selection_of_locations)

    % Nx
    figure;
    subplot(3,3,1); plot(Nx(1,timestart:timeend)); title('log Nx{1,1}')
    subplot(3,3,2); plot(Nx(4,timestart:timeend)); title('log Nx{1,4}')
    subplot(3,3,3); plot(Nx(7,timestart:timeend)); title('log Nx{1,7}')
    subplot(3,3,4); plot(Nx(22,timestart:timeend)); title('log Nx{4,1}')
    subplot(3,3,5); plot(Nx(25,timestart:timeend)); title('log Nx{4,4}')
    subplot(3,3,6); plot(Nx(28,timestart:timeend)); title('log Nx{4,7}')
    subplot(3,3,7); plot(Nx(43,timestart:timeend)); title('log Nx{7,1}')
    subplot(3,3,8); plot(Nx(46,timestart:timeend)); title('log Nx{7,4}')
    subplot(3,3,9); plot(Nx(49,timestart:timeend)); title('log Nx{7,7}')


    % Kx
    figure;
    subplot(3,3,1); plot(Kx(1,timestart:timeend)); title('log Kx{1,1}')
    subplot(3,3,2); plot(Kx(4,timestart:timeend)); title('log Kx{1,4}')
    subplot(3,3,3); plot(Kx(7,timestart:timeend)); title('log Kx{1,7}')
    subplot(3,3,4); plot(Kx(22,timestart:timeend)); title('log Kx{4,1}')
    subplot(3,3,5); plot(Kx(25,timestart:timeend)); title('log Kx{4,4}')
    subplot(3,3,6); plot(Kx(28,timestart:timeend)); title('log Kx{4,7}')
    subplot(3,3,7); plot(Kx(43,timestart:timeend)); title('log Kx{7,1}')
    subplot(3,3,8); plot(Kx(46,timestart:timeend)); title('log Kx{7,4}')
    subplot(3,3,9); plot(Kx(49,timestart:timeend)); title('log Kx{7,7}')

    % Ix
    figure;
    subplot(3,3,1); plot(Ix(1,timestart:timeend)); title('log Ix{1,1}')
    subplot(3,3,2); plot(Ix(4,timestart:timeend)); title('log Ix{1,4}')
    subplot(3,3,3); plot(Ix(7,timestart:timeend)); title('log Ix{1,7}')
    subplot(3,3,4); plot(Ix(22,timestart:timeend)); title('log Ix{4,1}')
    subplot(3,3,5); plot(Ix(25,timestart:timeend)); title('log Ix{4,4}')
    subplot(3,3,6); plot(Ix(28,timestart:timeend)); title('log Ix{4,7}')
    subplot(3,3,7); plot(Ix(43,timestart:timeend)); title('log Ix{7,1}')
    subplot(3,3,8); plot(Ix(46,timestart:timeend)); title('log Ix{7,4}')
    subplot(3,3,9); plot(Ix(49,timestart:timeend)); title('log Ix{7,7}')

end
%% time series for all locations
if logical(time_series_for_all_locations)
    
    % Nx 
    figure;
    for ii=1:opts.SpatialPointsPerDimension*opts.SpatialPointsPerDimension
        subplot(opts.SpatialPointsPerDimension,opts.SpatialPointsPerDimension,ii); plot(Ix(ii,timestart:timeend)); title(['log Nx{',num2str(ii),',',num2str(jj),'}'])
    end
    
end