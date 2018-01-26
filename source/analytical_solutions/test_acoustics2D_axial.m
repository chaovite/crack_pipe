% acoustics2D_axial_sym_driver
addpath('../helper');
%% pressure boundary condition.
a     = 1;   % position to apply boundary condition.
c0   = 1e3;   % acoustic wave speed.
rho = 1e3;   % density
r     = a;   % query point position
d    = 1e0;  % crack full width.
mu = 0;       % fluid viscosity.
dt = 0.01;
t   = [0:dt:1000];
t0 = 1;
T  = 0.1;
g  = exp(-(t-t0).^2/T^2);
g_ps = exp(-(t-t0).^2/T^2)*rho/d;

% plot(t,g);shg
[ghat, f] = fft_dim(g,dt);
% plot(f,abs(ghat));
 [p, v,  t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu, 'p');
 assert(max(abs(imag(p)))<1e-10);
p = real(p);
v = real(v);
 [p_ps,  t_ps] = acoustics2D_pointsource(r, c0, g_ps, dt);
 plot(t,p);
%% velocity boundary condition:
g  = exp(-(t-t0).^2/T^2)/(2*pi*d*a);
g_ps = exp(-(t-t0).^2/T^2)*rho/d;
R = [a:a:10*a];
P = cell(length(R), 1);
V = cell(length(R), 1);
P_ps = cell(length(R), 1);
for i = 1:length(R)
    r = R(i);
    [P{i}, V{i},  t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu,'v');
    [P_ps{i},  t_ps] = acoustics2D_pointsource(r, c0, g_ps, dt);
    assert(max(abs(imag(P{i})))<1e-10);
    assert(max(abs(imag(P_ps{i})))<1e-10);
    P_ps{i} = real(P_ps{i});
    P{i} = real(P{i});
    V{i} = real(V{i});
end
%%
figure;
for i = 1:length(R)
    r = R(i);
    p = P{i};
    p_ps = P_ps{i};
    plot(t, p + (i-1)*30,'k-');
    hold on;
    plot(t_ps, p_ps + (i-1)*30,'r-.');
    ylabel('p');
    hold on;
end
hold off;
xlim([0 6]);
xlabel('time (sec)');
ylabel('p (Pa)');
title('Solution compare, black: velocity bc, red: point source');



