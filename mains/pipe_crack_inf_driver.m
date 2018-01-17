% pipe_crack_inf_driver
source_dir = '../source';
addpath(genpath(source_dir));
%%
dt = 0.001;
t   = [0:dt:100];
t0 = 1;
T  = 0.1;
g  = exp(-(t-t0).^2/T^2);
figure
plot(t,g);
[ghat, f] = fft_dim(g, dt);
figure
plot(f,abs(ghat));
%%
omega = 2*pi*f;
L    = 1000;
a    = 0.01;
d    = 0;
rho = 1e3;
c0  = 1e3;
alpha = 8*1e-3/a^2;
c    = c0*sqrt(1./(1+1i*alpha./(rho*omega))); 
z    = -L/2;
% r    = 2000;
mu = 1e-2;
[p, v, t, R, T, omega] = pipe_inf_crack_inf(L, a, d, rho, c0, g, dt, z, r, mu);
% assert(max(abs(imag(p)./real(p)))<1e-6);
p = real(p);
v = real(v);
subplot(2,1,1);
plot(t,p);
subplot(2,1,2);
plot(t,v);
