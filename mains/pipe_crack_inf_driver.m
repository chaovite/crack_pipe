% pipe_crack_inf_driver
% close all
source_dir = '../source';
addpath(genpath(source_dir));
%%
dt = 0.01;
t   = [0:dt:1000];
t0 = 1;
T  = 0.1;
g  = exp(-(t-t0).^2/T^2);
%%
L    =  2000;
a    =  3;
d    =  1e-1;
rho =  1e3;
cp0  = 2e3;
cc0  = 2.5e2;
z      = 0;
mu   = 1;
R  = [10:50:1000];
figure
for i = 1: length(R)
    r = R(i);
    [p, v, t,~, omega] = pipe_crack_inf(L, a, d, rho, cp0, cc0, g, dt, z, r, mu);
    p = real(p);
%     plot(t,p+(i-1)*0.05);
    plot(t,v+(i-1)*0.05/(rho*cc0));
    hold on;
    shg
end
% xlim([1, L/cp0 + R(end)/cc0 + t0]);

