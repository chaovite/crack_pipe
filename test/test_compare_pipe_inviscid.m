% tests to check if pipe1d and conduit produces same solution at inviscid
% limit with gravity.
%
source_dir = '../source';
addpath(genpath(source_dir));
%% pipe1d solution with one material interface.
clear
Mp.name  = 'pipe';
Mp.L         = [500, 500];        % length

% radius does not vary in depth for now.
Mp.R         = [10, 10];                                         % radius

Mp.rho      = [2800, 1800];  % density
Mp.c         =  [2000, 1200];     % wave speed 
Mp.K         = Mp.rho.*Mp.c.^2;                       % fluid bulk modulus.
Mp.mu      = 0;                                           %  viscosity.
Mp.S         = pi*Mp.R.^2;                               % pipe surface area.
Mp.nz        = [32,32];               % number grid points in each section.
Mp.g         = 10;                            % gravitational acceleration
Mp.order   = 4;                             % order of spatial discretization.

Mp.interfaces{1}.type = 'm';

Mp.bc.bt.type = 'p';
Mp.bc.bt.link  = 'bt';

Mp.pT.A = 5e3;  % pressure perturbation amplitude Pa.
Mp.pT.T = 0.25;      % pressure perturbation duration
Mp.pT.t  = 2;      % pressure perturbation center time

Mp.bc.tp.type = 'h'; 
Mp.bc.tp.f       = @(t) Mp.pT.A*exp(-0.5*((t-Mp.pT.t)/Mp.pT.T)^2); % external force from the conduit.

model = pipe1d(Mp);
%
% ------------------------------------Time integration-------------------------------------
CFL  = 0.5;
T       = 5;
plot_simu = true;

% time stepping
[cmax, hmin] = model.getCFL();
dt      = CFL*hmin/cmax;
nt      = ceil(T/dt);
skip   = floor(nt/20);

nkeep = floor(nt/skip);
pzs       = cell(nkeep, 1);
vzs       = cell(nkeep, 1);
tkeep   = zeros(nkeep, 1);

fun = @(u, t) model.fun_integrate(u, t);
tic
model=model.init();
u = model.u;
cnt = 0;

for i=1:nt
    t = (i-1)*dt;
    u = lsrk4(fun,u,t,dt);
    
    if plot_simu
        if mod(i, skip) == 0
            cnt = cnt+1;
            vz           = u(model.indu.indv);
            pz          = u(model.indu.indp);
            pzs{cnt, 1} = pz;
            vzs{cnt, 1} = vz;
            tkeep(cnt) = t+dt;
            z = model.grd.z;
            figure(1)
            subplot(2,1,1)
            plot(z, pz);
            ylim(2*[-Mp.pT.A, Mp.pT.A]);
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

d1.name = 'pipe1d';
d1.bctype = Mp.bc.bt.type;
d1.pz = pzs;
d1.vz = vzs;
d1.z     = z;
d1.t     = tkeep;
%% conduit2d solution with one material interface.
% conduit parameters
Mc.R = 10;
Mc.L = 1000;
Mc.nz = 64;
Mc.nr = 8;
Mc.order = 4;
Mc.S = pi*Mc.R^2;
Mc.g = 10;
Mc.mu = 0;
dz = Mc.L/Mc.nz;
z = dz*[0: Mc.nz]';

% conduit wave speed, density for the upper and lower section.
Mc.rho_upper = 1800;
Mc.rho_lower  = 2800;
Mc.c_upper     = 1200;
Mc.c_lower     = 2000;
Mc.L_upper    = 0.5*Mc.L;
Mc.L_lower    = Mc.L - Mc.L_upper;

Mc.interface_split = true;
Mc.with_exsolution = false;

Mc.split_index = round(Mc.L_lower/dz) + 1;
Mc.rho = [Mc.rho_lower*ones(Mc.split_index, 1); Mc.rho_upper*ones(Mc.nz+2 - Mc.split_index, 1)];
Mc.c    = [Mc.c_lower*ones(Mc.split_index, 1); Mc.c_upper*ones(Mc.nz+2 - Mc.split_index, 1)];
Mc.K   = Mc.rho.*(Mc.c).^2;

Mc.pT.A = 5e3; % pressure perturbation amplitude
Mc.pT.T = 0.25; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time

Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
%
model = conduit(Mc);

%
CFL = 0.5;
T     = 5;
plot_simu = true;

% time stepping
cmax = max(Mc.c);
hmin  = Mc.L/Mc.nz; 
dt = CFL*hmin/cmax;
nt      = ceil(T/dt);
skip   = floor(nt/20);
nkeep = floor(nt/skip);
pzs       = cell(nkeep, 1);
vzs       = cell(nkeep, 1);
tkeep   = zeros(nkeep, 1);

fun = @(u, t) model.A*u + model.Fp*Mc.G(t);
cnt = 0;
u = model.u;
for i=1:nt
    t = (i-1)*dt;
    % invicid, integrate explicitly
    u = lsrk4(fun,u,t,dt);
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end

    if plot_simu
        if mod(i, skip) == 0
            cnt = cnt + 1;
           indv = field_indices(model.dimensions(),1);
           indp = field_indices(model.dimensions(),2);
            vz          = u(indv);
            uz          = model.op.W2*vz;
            pz          = u(indp);
            
            pzs{cnt, 1} = pz;
            vzs{cnt, 1} = uz;
            tkeep(cnt) = t+dt;
            
            z = model.geom.z;
            figure(1)
            subplot(2,1,1)
            plot(z, pz,'k-');
            ylim(2*[-Mc.pT.A, Mc.pT.A]);
            xlim([0, Mc.L]);
            ylabel('pressure');
            title(sprintf('t = %8.2f', t));
            
            subplot(2,1,2)
            plot(z, uz,'k-');
            ylim([-0.01, 0.01]);
            xlim([0, Mc.L]);
            ylabel('velocity');
            drawnow();
        end
    end
end

d2.name = 'conduit2s';
d2.bctype = 'p';
d2.pz = pzs;
d2.vz = vzs;
d2.z     = z;
d2.t     = tkeep;
save('test_compare_pipe_inviscid_mat','d1','d2','nkeep');
%%
load('test_compare_pipe_inviscid_mat','d1','d2','nkeep');
for i = 1: nkeep
    subplot(2,1,1)
    plot(d1.z, d1.pz{i}, 'k-', d2.z,d2.pz{i},'r-');
    ylim([-5e3, 5e3]*2);
    legend('pipe1d','conduit2d');
    ylabel('pressure');
    subplot(2,1,2)
    plot(d1.z, d1.vz{i}, 'k-', d2.z, d2.vz{i},'r-');
    legend('pipe1d','conduit2d');
    ylim([-0.01, 0.01]);
    ylabel('velocity');
    pause();
end





