%plot acoustics2D_analytical
cd('/Users/Chaovite/Documents/code/magmagpu/projects/chao/test');
addpath('../helper/')
T = 10;
dt = 0.001;
T0 = 2;
t = 0:dt:T;
Tc = 0.05;
g = ricker(t, 2, T0);
c = 1;
L = 0.8*(T-T0)*c;
r = [0.001:0.001:1]*L;
N = length(r);
p = cell(N,1);
p = zeros(N, length(t)-1);
tic
for i = 1 : N
    [p(i,:), t] = acoustics2D_pointsource(r(i), c, g, dt);
end
toc
%%
% plot pressure distribution.
for i = 1:100: length(t); plot(r, real(p(:,i))); title(num2str(t(i)));ylim([-3 3]); pause(0.2);drawnow;end
%% plot pressure time series at different r.
for i = 1:100: N; plot(t, real(p(i,:)) + r(i),'k-'); title(num2str(r(i)/c)); hold on;end
hold off;
daspect([1 1 1]);
shg









