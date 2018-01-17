% acoustics2D_axial_sym_driver

%% pressure boundary condition.
a     = 1e-1;   % position to apply boundary condition.
c0   = 1e3;   % acoustic wave speed.
rho = 2e3;   % density
r     = 3e2;   % query point position
d    = 1;  % crack full width.
mu = 0;       % fluid viscosity.
dt = 0.01;
t   = [0:dt:100];
t0 = 1;
T  = 0.1;
g  = exp(-(t-t0).^2/T^2)/(2*pi*d*a);
g_ps = exp(-(t-t0).^2/T^2)/d;

plot(t,g);shg
[ghat, f] = fft_dim(g,dt);
figure;
plot(f,abs(ghat));
 [p, v,  t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu,'v');
 [p_ps,  t_ps] = acoustics2D_pointsource(r, c0, g_ps, dt);
subplot(2,1,1);
 plot(t, p./max(p), t_ps, p_ps./max(p_ps));
 subplot(2,1,2);
plot(t,v);

%% velocity boundary condition:
dt = 0.01;
t   = [0:dt:100];
t0 = 1;
T  = 0.1;
g  = exp(-(t-t0).^2/T^2)/rho/c0/(2*pi*a*d);

R = [1:1:100];
P = cell(length(R), 1);
V = cell(length(R), 1);
for i = 1:length(R)
    r = R(i);
    [P{i}, V{i},  t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu,'v');
    assert(max(abs(imag(P{i})))<1e-10);
    P{i} = real(P{i});
    V{i} = real(V{i});
end
%
figure;
for i = 1:length(R)
    r = R(i);
    p = P{i};
    v = V{i};
    subplot(1,2,1)
    plot(t, p + i*0.01,'k-');
%     xlim([0 7]);
    hold on;
     ylabel('p');
    subplot(1,2,2);
    plot(t, v + i*0.01/(rho*c0),'k-');
    ylabel('v');
    hold on;
%     xlim([0 7]);
end
hold off;
% xlim([0 7]);
%% compare solution with acoustics2D with point source.
g_ps = g*2*pi*d*a/(pi*a^2*d)*rho*c0^2;
P_ps = cell(length(R), 1);
for i = 1:length(R)
    r = R(i);
    [P_ps{i},  t_ps] = acoustics2D_pointsource(r, c0, g_ps, dt);
    assert(max(abs(imag(P_ps{i})))<1e-10);
    P_ps{i} = real(P_ps{i});
end
%%
figure;
for i = 1:length(R)
    r = R(i);
    p = P{i};
    p_ps = P_ps{i};
%     plot(t, p + i*0.01,'k-');
    plot(t_ps, p_ps + i*0.01,'r-');
     ylabel('p');
     hold on;
end
hold off;
xlim([0 7]);



