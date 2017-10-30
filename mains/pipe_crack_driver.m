% pipe_crack_flow_driver, multi-section pipe model coupled to multiple cracks.
source_dir = '../source';
addpath(genpath(source_dir));
%%
%-----------------------------------pipe parameter-----------------------------------
clear
% three section pipe. [0 200], [200, 400], [400, 1000]
Mp.name  = 'pipe';
Mp.L         = [300, 300, 400];        % length

% radius does not vary in depth for now.
Mp.R         = [20];                                         % radius

Mp.rho      = [2800, 2800, 1200];  % density
Mp.c         = [1800, 1800, 1000];     % wave speed 
Mp.K         = Mp.rho.*Mp.c.^2;                       % fluid bulk modulus.
Mp.mu      = [20];                                           %  viscosity.
Mp.S         = pi*Mp.R.^2;                               % pipe surface area.
Mp.nz       = [8, 8, 8];               % number grid points in each section.
Mp.g         = 10;                            % gravitational acceleration
Mp.order  = 4;                             % order of spatial discretization.

% interface_type: 
% 'c': crack (v+ - v- = -Q/A, p+ = p- = pc)
% 'm': material properties (v+ = v-, p+ = p-)
% if the interface is a crack, then must specify the name of the crack that it links to.

% Mp.interfaces = [];
Mp.interfaces{1}.name  = 'interface_1';
Mp.interfaces{1}.type    = 'c';
Mp.interfaces{1}.link     = 'interface_1'; 

% Mp.interfaces{1}.name  = 'interface_2';
Mp.interfaces{2}.type    =  'm';

% boundary conditions for x=0 and x=L 
% 'c': crack.
% 'h': moving surface, with extra pressure perturbation.
% 'p': pressure boundary condition, p = 0
% 'v' : velocity boundary condition, v = 0

% Note I only implement moving surface on top bc.
% bt bc can be p=0, v=0 or coupled to crack.

Mp.bc.bt.type = 'c';
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

pipe = pipe1d(Mp);
%% -----------------------------------crack parameters-----------------------------------
 % define the properties of crack.
% fluid-filled fracture parameters.

% the name of the interface of the pipe that this crack is coupled to. 
Mc.name = 'interface_1';
Mc.w0      = 2;
Mc.Lx       = 2e3;
Mc.Ly       = 2e3;
Mc.nx       = 32;
Mc.ny       = 32;
Mc.order = 4;
Mc.interp_order = 4;

% Injecting source location at the crack.
Mc.xs = 0.75*Mc.Lx;
Mc.ys = 0.5*Mc.Ly;
Mc.G  = @(t) 0;

% coupling location at the crack.
Mc.xc = 0.5*Mc.Lx;% the coupling location to the conduit.
Mc.yc = 0.5*Mc.Ly;% the coupling location to the conduit.

% if the crack is rigid.
Mc.isrigid = false;
Mc.use_fft = false;
Mc.npad_fft = 256;

% fluid properties, continuous with the conduit bottom.
Mc.rho  = Mp.rho(2);
Mc.c     = Mp.c(2);
Mc.K     = Mc.rho*Mc.c^2;
Mc.mu = 20;

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

% properties of a bottom crack.
Mc_bt = Mc;
Mc_bt.name = 'bt';
Mc_bt.Zc = 1000;
Mc_bt.Lx = 2000;
Mc_bt.Ly = 2000;
Mc_bt.rho  = Mp.rho(1);
Mc_bt.c     = Mp.c(1);
Mc_bt.K     = Mc_bt.rho*Mc_bt.c^2;
Mc_bt.xc = 0.75*Mc_bt.Lx;% the coupling location to the conduit.
Mc_bt.yc = 0.5*Mc_bt.Ly;% the coupling location to the conduit.

% building crack model.
% East, North coordinate of observation w.r.t to the centroid.
x_obs = 0;
y_obs = 1000;

cracks.(Mc.name)     = frac2dxy(Mc);
cracks.(Mc_bt.name) = frac2dxy(Mc_bt);
%%
tic; model = pipeCracks(pipe, cracks); toc
%% ------------------------------------Time integration-------------------------------------
CFL = 0.5;
skip = 20;
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
            vz           = u(model.indu.pipe.vz);
            pz          = u(model.indu.pipe.pz);
            Xp_bt     = model.cracks.bt.grd.p.X;
            Yp_bt     = model.cracks.bt.grd.p.Y;
            pxy_bt    = model.cracks.bt.grd.p.grd(u(model.indu.bt.p));
            name = 'interface_1';
            
            Xp    = model.cracks.(name).grd.p.X;
            Yp    = model.cracks.(name).grd.p.Y;
            pxy   = model.cracks.(name).grd.p.grd(u(model.indu.(name).p));
            
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
            view([1 1 1])
            drawnow();
            
            figure(3)
            pcolor(Xp,Yp, pxy);
            shading interp;
            cmap;
            colorbar;
            caxis([-Mp.pT.A, Mp.pT.A]/2);
            title('Middle Crack');
            axis off
            view([1 1 1])
           drawnow();
        end
    end
    
end






