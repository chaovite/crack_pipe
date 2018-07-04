% test_conduit_interal_g inviscid incompressible
% compare the inviscid period and that from analytical solution.
% It's not absolutely incompressible but can be approximated very well.
%

tic
sourcedir = '../source';
addpath(genpath(sourcedir));

Mc.R = 5; 
Mc.L = 300; 
Mc.nz = 20000; 
Mc.nr = 12; 
Mc.order = 6;
Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1); 
Mc.g = 10;  
Mc.with_exsolution=false;
Mc.interface_split=false;

% inviscid.
Mc.mu = 0*ones(Mc.nz+1, 1);
z = Mc.L/Mc.nz*[0:Mc.nz]';

% density profile
rho0 = 500;
rho1 = 2000;
alpha = (log(rho1)-log(rho0))/Mc.L;
Mc.rho = rho0*exp(alpha*(Mc.L-z));
Mc.c    = 20000*ones(Mc.nz+1,1);
Mc.K    = Mc.rho.*Mc.c.^2;
Mc.Mg = alpha - Mc.rho*Mc.g./Mc.K;
%
Mc.epsilon = 0.3;
cond_g = conduit_internal_g(Mc);

% resonant period from analytical solution.
rho_m = mean(Mc.rho);
drho    =   rho1-(1-Mc.epsilon)*rho0;
g_prime = drho/rho_m*10;
omega0 = sqrt(g_prime/Mc.L);
%
% solve the eigenvalue problem.
omega = eigs(cond_g.A, 1, omega0*1i);
fprintf('Analytical solution T=%8.4f, numerical solution T=%8.4f \n',2*pi/omega0, 2*pi/abs(imag(omega)));
toc