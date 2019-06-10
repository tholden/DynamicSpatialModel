omegar = 100;
rbar = 0.04;
r = [-0.02:0.00001:0.08];

P_true = zeros(size(r)); P_true(r<=0)=max(P);
P = exp(-omegar .* r) ./ omegar;
P_2nd = exp( -omegar * rbar ) / omegar + ( r - rbar ) .* (-exp(-omegar * rbar)) + 0.5*( r - rbar ).^2 .* omegar * exp(-omegar * rbar);
P_3rd = P_2nd - (1/6)*( r - rbar ).^3 .* omegar^2 * exp(-omegar * rbar);


figure;
plot(r,P_true); hold on; 
plot(r,P); 
plot(r,P_2nd);
plot(r,P_3rd);


legend({'True','global penalty function','2nd order approx','3rd order approx'});
xlabel('r'); ylabel('P');
title('omega = 100, steady state r = 0.04')
ylim([-0.002,P(r==0)]);
