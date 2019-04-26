% test_compare_coupled_model solution
%
% test if the frac3d+conduit and frac2dxy+pipe1d produces the same
% solution in time domain under the limit of zero viscosity.
%

source_dir = '../source';
addpath(genpath(source_dir));
%% frac2dxy+pipe1d pipCrack solution.
clear
% three section pipe. [0 200], [200, 400], [400, 1000]
Mp.name  = 'pipe';
Mp.L         = 1000;        % length

% radius does not vary in depth for now.
Mp.R         = 10;                                         % radius

Mp.rho      = 1800;  % density
Mp.c         =  1200;     % wave speed 
Mp.K         = Mp.rho.*Mp.c.^2;                       % fluid bulk modulus.
Mp.mu      = 0;                                           %  viscosity.
Mp.S         = pi*Mp.R.^2;                               % pipe surface area.
Mp.nz        = 64;               % number grid points in each section.
Mp.g         = 0;                            % gravitational acceleration
Mp.order   = 2;                             % order of spatial discretization.

Mp.interfaces = [];

Mp.bc.bt.type = 'c';
Mp.bc.bt.link  = 'bt';

Mp.pT.A = 5e3;  % pressure perturbation amplitude Pa.
Mp.pT.T = 0.25;      % pressure perturbation duration
Mp.pT.t  = 2;      % pressure perturbation center time

Mp.bc.tp.type = 'h'; 
Mp.bc.tp.f       = @(t) Mp.pT.A*exp(-0.5*((t-Mp.pT.t)/Mp.pT.T)^2); % external force from the conduit.

pipe = pipe1d(Mp);
% -----------------------------------crack parameters-----------------------------------
 % define the properties of crack.
% fluid-filled fracture parameters.

% the name of the interface of the pipe that this crack is coupled to. 
Mc.name  = 'bt';
Mc.w0      = 2;
Mc.Lx       = 2e3;
Mc.Ly       = 2e3;
Mc.nx       = 32;
Mc.ny       = 32;
Mc.order   = 2;
Mc.interp_order = 2;

% Injecting source location at the crack.
Mc.xs = 0.75*Mc.Lx;
Mc.ys = 0.5*Mc.Ly;
Mc.G  = @(t) 0;

% coupling location at the crack.
Mc.xc = 0.5*Mc.Lx;% the coupling location to the conduit.
Mc.yc = 0.5*Mc.Ly;% the coupling location to the conduit.

% if the crack is rigid.
Mc.isrigid = true;
Mc.use_fft = false;
Mc.npad_fft = 256;

% fluid properties, continuous with the conduit bottom.
Mc.rho  = Mp.rho(1);
Mc.c     = Mp.c(1);
Mc.K     = Mc.rho*Mc.c^2;
Mc.mu = 0;

% solid properties.
Mc.cp    = 5e3;
Mc.cs    = 2.7e3;
Mc.rhos = 3e3;
Mc.Gs   = Mc.cs^2*Mc.rhos;
Mc.nu  = ((Mc.cp/Mc.cs)^2-2)/((Mc.cp/Mc.cs)^2-1)/2;

% crack centroid location.
Mc.Xc = 0;
Mc.Yc = 0;
Mc.Zc = 1000 ;
Mc.strike = 0;
Mc.dip = 0;

% building crack model.
% East, North coordinate of observation w.r.t to the centroid.
x_obs = 0;
y_obs = 1000;

cracks.(Mc.name) = frac2dxy(Mc);
%
tic; model = pipeCracks(pipe, cracks); toc
% ------------------------------------Time integration-------------------------------------
CFL  = 0.5;
skip   = 10;
T       = 5;
plot_simu = true;

% time stepping
[cmax, hmin] = model.getCFL();
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

nkeep = floor(nt/skip);
pzs       = cell(nkeep, 1);
vzs       = cell(nkeep, 1);
pxys     = cell(nkeep, 1);

% [L,U,p,q,B] = imex_ark4_get_lu(model.Ai, dt);
fun = @(u, t) model.fun_integrate(u, t) + model.Ai*u;
tic
Us = zeros(3, nt);
model=model.init();
u = model.u;
cnt = 0;

for i=1:nt
    t = (i-1)*dt;
    u=lsrk4(fun, u,t,dt);
%     u = imex_ark4_lu(u,t,dt,fun, model.Ai, L, U, p, q);
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i, skip) == 0
            cnt = cnt+1;
            vz           = u(model.indu.pipe.vz);
            pz          = u(model.indu.pipe.pz);
            Xp_bt     = model.cracks.bt.grd.p.X;
            Yp_bt     = model.cracks.bt.grd.p.Y;
            pxy_bt    = model.cracks.bt.grd.p.grd(u(model.indu.bt.p));
            pzs{cnt, 1} = pz;
            vzs{cnt, 1} = vz;
            pxys{cnt, 1}= pxy_bt;
            
            z = model.pipe.grd.z;
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
            drawnow;
            
            figure(2)
            pcolor(Xp_bt,Yp_bt, pxy_bt);
            shading interp;
            cmap;
            colorbar;
            caxis([-Mp.pT.A, Mp.pT.A]/2);
            axis off
            title('Bottom Crack');
           drawnow();
        end
    end
end
d1.pz = pzs;
d1.vz = vzs;
d1.pxy = pxys;
d1.z     = z;
d1.X    = Xp_bt;
d1.Y    = Yp_bt;
%% frac3d+conduit solution
clear Mc
% conduit parameters
Mc.R   = 10;
Mc.L   = 1000;
Mc.nz = 64;
Mc.nr  = 6;
Mc.order = 2;
Mc.order_r = 2;
Mc.S = pi*Mc.R^2;
Mc.g = 0;
Mc.mu = 0;
Mc.interface_split = false;
Mc.with_exsolution = false;
Mc.rho = 1800 * ones(Mc.nz+1, 1);
Mc.c    = 1200 * ones(Mc.nz+1, 1);
Mc.K    = Mc.rho.*Mc.c.^2;

Mc.pT.A = 5e3; % pressure perturbation amplitude
Mc.pT.T = 0.25; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
% Mc.G = @(t) 0; % external force from the conduit.

% fluid-filled fracture parameters.
Mf.w0 = 2;
Mf.Lx   = 2e3;
Mf.Ly   = 2e3;
Mf.nx = 32;
Mf.ny = 32;
Mf.nz = 16;
Mf.order = 2;
Mf.order_z = 2;
Mf.interp_order = 2;
Mf.xs = 0.75*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) 0;
Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = true;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.15; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho = Mc.rho(1);
Mf.c    = Mc.c(1);
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 0;

Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu   = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
Mf.Xc   = 0;
Mf.Yc   = 0;
Mf.Zc   = 1000 ;
Mf.strike = 0;
Mf.dip = 0;

x_obs = 0;
y_obs = 1000;
% construct coupled models.
cond = conduit(Mc);
frac  = frac3d_o(Mf);
Model=coupledModel(cond, frac);
% time stepping:
CFL = 0.5;
skip = 10;
T = 5;
use_imex = false;
plot_simu = true;

% time stepping
hmin = min(Model.conduit.geom.dz, Model.frac.geom.p.hx);
cmax = max(max(Model.conduit.M.c), Model.frac.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

nkeep = floor(nt/skip);
pzs       = cell(nkeep, 1);
vzs       = cell(nkeep, 1);
pxys     = cell(nkeep, 1);

Model=Model.update(zeros(sum(Model.dimensions()), 1));
% surface displacement.
Us = zeros(3, nt);

if ~ use_imex
    A = Model.Ae + Model.Ai;
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end
fun = @(u,t) A*u + Model.Fp(:,1)*Model.conduit.M.G(t) + Model.Fp(:,2)*Model.frac.M.G(t);

cnt = 0;

tic
for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
    end
    
    if plot_simu
        if mod(i,skip) == 0
            cnt = cnt + 1;
            %plot solution.
            if Mc.with_exsolution
                [vz, pz, ~,~, pxy, vxy] = Model.fields(Model.u);
            else
                [vz, pz,~, pxy, ~] = Model.fields(Model.u);
            end
            uz = Model.conduit.op.W2*Model.field(Model.u, [1, 1]);
            
            pzs{cnt, 1} = pz;
            vzs{cnt, 1} = uz;
            pxys{cnt, 1}= pxy;
            
            z = Model.conduit.geom.z;
            X = Model.frac.geom.p.X;
            Y = Model.frac.geom.p.Y;

            figure(1)
            subplot(2,1,1)
            plot(z, pz);
            ylim(2*[-Mc.pT.A, Mc.pT.A]);
            xlim([0 sum(Mc.L)])
            title(sprintf('t = %8.2f', t));
            subplot(2,1,2)
            plot(z, uz);
            ylim([-0.01 0.01]);
            xlim([0 sum(Mc.L)])
            title(sprintf('t = %8.2f', t));
            drawnow;
            
            figure(2)
            pcolor(X, Y, pxy);
            shading interp;
            cmap;
            colorbar;
            caxis([-Mc.pT.A, Mc.pT.A]/2);
            axis off
            title('Bottom Crack');
           drawnow();
        end
        
    end
end

d2.pz = pzs;
d2.vz = vzs;
d2.pxy = pxys;
d2.z     = z;
d2.X    = X;
d2.Y    = Y;
%%
% save('compare_coupledmodel_inviscid_rigid_crack','d1','d2','nkeep');
%%
% load('compare_coupledmodel_inviscid_rigid_crack','d1','d2','nkeep');
nkeep = length(d1.pz);
% compare conduit pressure and velocity.
figure(1);
set(gcf,'position',[0 0 1200, 1200])
for i=1: nkeep
    subplot(4,2,[1,2])
    plot(d1.z, d1.pz{i,1},'k-',d2.z, d2.pz{i,1},'r-');
    ylim([-5e3, 5e3]);
    title('pressure');
    legend('pipeCrack','coupledModel');
    
    subplot(4,2,[3,4])
    plot(d1.z, d1.vz{i,1},'k-',d2.z, d2.vz{i,1},'r-');
     ylim([-0.01 0.01]);
    title('velocity');
    legend('pipeCrack','coupledModel');
   
    subplot(4,2,[5,7])
    pcolor(d1.X, d1.Y, d1.pxy{i});
    shading interp;
    cmap;
    colorbar;
    caxis([-5e3, 5e3]/2);
    title('pipeCrack');
    
    subplot(4,2,[6,8])
    pcolor(d2.X, d2.Y, d2.pxy{i});
    shading interp;
    cmap;
    colorbar;
    caxis([-5e3, 5e3]/2);
    title('coupledModel');
    pause(0.05);
end



