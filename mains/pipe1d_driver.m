% pipe1d_driver, multi-section pipe model.
source_dir = '../source';
addpath(genpath(source_dir));
%%
%-----------------------------------pipe parameter-----------------------------------
clear
% three section pipe. [0 200], [200, 400], [400, 1000]
Mp.name  = 'pipe';
Mp.L         = [500, 500];        % length

% radius does not vary in depth for now.
Mp.R         = [20];                                         % radius

Mp.rho      = [2800, 1000];  % density
Mp.c         = [1800, 1000];     % wave speed 
Mp.K         = Mp.rho.*Mp.c.^2;                       % fluid bulk modulus.
Mp.mu      = [20];                                           %  viscosity.
Mp.S         = pi*Mp.R.^2;                               % pipe surface area.
Mp.nz       = [32, 32];               % number grid points in each section.
Mp.g         = 10;                            % gravitational acceleration
Mp.order  = 4;                             % order of spatial discretization.

% interface_type: 
% 'c': crack (v+ - v- = -Q/A, p+ = p- = pc)
% 'm': material properties (v+ = v-, p+ = p-)
% if the interface is a crack, then must specify the name of the crack that it links to.

% Mp.interfaces = [];

Mp.interfaces{1}.type    =  'm';

% boundary conditions for x=0 and x=L 
% 'c': crack.
% 'h': moving surface, with extra pressure perturbation.
% 'p': pressure boundary condition, p = 0
% 'v' : velocity boundary condition, v = 0

% Note I only implement moving surface on top bc.
% bt bc can be p=0, v=0 or coupled to crack.

Mp.bc.bt.type = 'p';
Mp.bc.bt.link  = 'bt';
% Mp.bc.bt.type = 'p';
% Mp.bc.bt.f       = @(t) 0;
% Mp.bc.bt.type  = 'v';
% Mp.bc.bt.f       = @(t) 0;

Mp.pT.A = 5e3;  % pressure perturbation amplitude Pa.
Mp.pT.T = 0.25;      % pressure perturbation duration
Mp.pT.t  = 2;      % pressure perturbation center time

Mp.bc.tp.type = 'h'; 
Mp.bc.tp.f       = @(t) Mp.pT.A*exp(-0.5*((t-Mp.pT.t)/Mp.pT.T)^2); % external force from the conduit.

model = pipe1d(Mp);
%% ------------------------------------Time integration-------------------------------------
CFL = 0.5;
skip = 10;
T = 200;
plot_simu = true;

% time stepping
[cmax, hmin] = model.getCFL();
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

[L,U,p,q,B] = imex_ark4_get_lu(model.Ai, dt);
fun = @(u, t) model.fun_integrate(u, t);
tic
Us = zeros(3, nt);
itr = 0;
model.init();
u = model.u;
for i=1:nt
    t = (i-1)*dt;
    u = imex_ark4_lu(u,t,dt,fun, model.Ai, L, U, p, q);
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i, skip) == 0
            vz           = u(model.indu.indv);
            pz          = u(model.indu.indp);
            eta         = u(model.indu.indetas);
            z = model.grd.z;
            figure(1)
            subplot(2,1,1)
            plot(z, pz);
            ylim([-Mp.pT.A, Mp.pT.A]*2);
            xlim([0 sum(Mp.L)])
            title(sprintf('t = %8.2f', t));
            subplot(2,1,2)
            plot(z, vz);
            ylim([-0.01 0.01]);
            xlim([0 sum(Mp.L)])
            title(sprintf('t = %8.2f', t));
           drawnow();
        end
    end
    
end






