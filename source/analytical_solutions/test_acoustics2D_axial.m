% acoustics2D_axial_sym_driver

%% pressure boundary condition.
a     = 1e-2;   % position to apply boundary condition.
c0   = 1e3;   % acoustic wave speed.
rho = 2e3;   % density
r     = 1e3;   % query point position
d    = 1;  % crack full width.
mu = 0;       % fluid viscosity.
dt = 0.01;
t   = [0:dt:100];
t0 = 1;
T  = 0.1;
g  = exp(-(t-t0).^2/T^2)/(2*pi*d*a);
g_ps = exp(-(t-t0).^2/T^2)*rho/d;

% plot(t,g);shg
[ghat, f] = fft_dim(g,dt);
% plot(f,abs(ghat));
 [p, v,  t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu,'v');
 assert(max(abs(imag(p)))<1e-10);
p = real(p);
v = real(v);
 [p_ps,  t_ps] = acoustics2D_pointsource(r, c0, g_ps, dt);
%% velocity boundary condition:
R = [1e2:2e2:4e3];
P = cell(length(R), 1);
V = cell(length(R), 1);
P_ps = cell(length(R), 1);
for i = 1:length(R)
    r = R(i);
    [P{i}, V{i},  t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu,'v');
    [P_ps{i},  t_ps] = acoustics2D_pointsource(r, c0, g_ps, dt);
    assert(max(abs(imag(P{i})))<1e-10);
    P{i} = real(P{i});
    V{i} = real(V{i});
end
%%
figure;
for i = 1:length(R)
    r = R(i);
    p = P{i};
    p_ps = P_ps{i};
    plot(t, p + i*200,'k-');
    hold on;
    plot(t_ps, p_ps + i*200,'r-.');
     ylabel('p');
     hold on;
end
hold off;
xlim([0 10]);



