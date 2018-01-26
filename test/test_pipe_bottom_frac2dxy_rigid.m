% test_pipe_bottom_crack_rigid
%% test the coupling of pipe and bottom crack use analytical solution.
% crack is radial symmetric, inifinite and rigid in the analytical solution.
% gravity in the turn-off.
clear;
source_dir = '../source';
addpath(genpath(source_dir));
%% pipe_crack model
%-----------------------------------pipe parameter-----------------------------------
clear
% three section pipe. [0 200], [200, 400], [400, 1000]
Mp.name  = 'pipe';
Mp.L         = 2000;        % length

% radius does not vary in depth for now.
Mp.R         = 1;                                         % radius
Mp.rho      = 1000;  % density
Mp.c         =  2000;     % wave speed 
Mp.K         = Mp.rho.*Mp.c.^2;                       % fluid bulk modulus.
Mp.mu      = 1e-3;                                           %  viscosity.
Mp.S         = pi*Mp.R.^2;                               % pipe surface area.
Mp.nz       = 128;               % number grid points in each section.
Mp.g         = 0;                            % gravitational acceleration
Mp.order  = 8;                             % order of spatial discretization.

% interface_type: 
% 'c': crack (v+ - v- = -Q/A, p+ = p- = pc)
% 'm': material properties (v+ = v-, p+ = p-)
% if the interface is a crack, then must specify the name of the crack that it links to.

Mp.interfaces = [];
Mp.bc.bt.type = 'c';
Mp.bc.bt.link  = 'bt';

Mp.pT.A = 1;  % pressure perturbation amplitude Pa.
Mp.pT.T = 0.1;   % pressure perturbation duration
Mp.pT.t  = 2;      % pressure perturbation center time

Mp.bc.tp.type = 'h'; 
Mp.bc.tp.f       = @(t) Mp.pT.A*exp(-0.5*((t-Mp.pT.t)/Mp.pT.T).^2); % external force from the conduit.

pipe = pipe1d(Mp);
% -----------------------------------crack parameters-----------------------------------
 % define the properties of crack.
% fluid-filled fracture parameters.

Mc.name = 'bt';
Mc.w0      = 1e-1;
Mc.Lx       = 2e3;
Mc.Ly       = 2e3;
Mc.nx       = 128;
Mc.ny       = 128;
Mc.order  = 8;
Mc.interp_order = 20;

% Injecting source location at the crack.
Mc.xs = 0.75*Mc.Lx;
Mc.ys = 0.5*Mc.Ly;
Mc.G  = @(t) 0;

% coupling location at the crack.
Mc.xc = 0.5*Mc.Lx;% the coupling location to the conduit.
Mc.yc = 0.5*Mc.Ly;% the coupling location to the conduit.

% if the crack is rigid.
Mc.isrigid = true;

% fluid properties, continuous with the conduit bottom.
Mc.rho  = 1000;
Mc.c      = 250;
Mc.K     = Mc.rho*Mc.c^2;
Mc.mu  = 1e-3;

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
cracks.(Mc.name)     = frac2dxy(Mc);
%
tic; model = pipeCracks(pipe, cracks); toc

% query points.
nq = 10;
% points in the pipe:
zq_idx = round(linspace(1,Mp.nz, nq))';
% analytica solution position downward with z=0 at the bottom.
% numerical solution position upward with z=0 at the bottom.
zq       = - full(model.pipe.grd.z(zq_idx));

% query points in the crack:
xq_idx = round(linspace(Mc.nx/2 + 0.1 *  Mc.nx/2, Mc.nx/2 + 1 + 0.75 *  Mc.nx/2, nq))';
xq        = model.cracks.bt.grd.p.x(xq_idx)' - Mc.Lx/2;
yq_idx = (Mc.nx/2 + 1)*ones(nq, 1);
yq       =  model.cracks.bt.grd.p.y(yq_idx)' - Mc.Ly/2;
rq        = sqrt(xq.^2 + yq.^2);
rq_idx = sub2ind(size(model.cracks.bt.grd.p.X), yq_idx, xq_idx);

% ------------------------------------Time integration-------------------------------------
CFL = 0.5;
skip = 20;
T     =  Mp.L/Mp.c +  0.5*Mc.Lx/Mc.c + Mp.pT.t;
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

pz_q   = zeros(nq, nt+1); 
vz_q   = zeros(nq, nt+1); 
pxy_q = zeros(nq, nt+1);
t_n     = [0:dt:nt*dt];

for i=1:nt
    t = (i-1)*dt;
%      u = lsrk4(fun,u,t,dt);
    u = imex_ark4_lu(u,t,dt,fun, model.Ai, L, U, p, q);
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    vz      = u(model.indu.pipe.vz);
    pz      = u(model.indu.pipe.pz);
    Xp     = model.cracks.bt.grd.p.X;
    Yp     = model.cracks.bt.grd.p.Y;
    pxy    = model.cracks.bt.grd.p.grd(u(model.indu.bt.p));
    z = model.pipe.grd.z;
    
    
    % record the solutions at the query points.
    pz_q(:, i+1)   = pz(zq_idx);
    vz_q(:, i+1)   = vz(zq_idx); 
    pxy_q(:, i+1) = pxy(rq_idx);
    
    if plot_simu
        if mod(i, skip) == 0   
            z = model.pipe.grd.z;
            figure(1)
            subplot(2,1,1)
            plot(z, pz,'k-');
            hold on;
            plot(-zq, zeros(nq, 1),'k*');
            hold off;
            ylim(2*[-Mp.pT.A, Mp.pT.A]);
            xlim([0 sum(Mp.L)])
            title(sprintf('t = %8.2f', t));
            subplot(2,1,2)
            plot(z, vz);
            ylim(2*[-Mp.pT.A, Mp.pT.A]/(Mp.rho*Mp.c));
            xlim([0 sum(Mp.L)])
            title(sprintf('t = %8.2f', t));
            drawnow;
            
            figure(2)
            pcolor(Xp,Yp, pxy);
            hold on;
            plot(xq+Mc.Lx/2, yq+Mc.Ly/2,'k*');
            hold off;
            shading interp;
            cmap;
            colorbar;
            caxis([-Mp.pT.A, Mp.pT.A]/4);
            axis off
            title('Bottom Crack');
            drawnow();
        end
    end
    
end
%% analytical solution:
dt_a = 0.01;
t_a   = [0:dt_a:2000];
g  = Mp.bc.tp.f(t_a);

pz_a   = cell(nq, 1);
vz_a   = cell(nq, 1);
pxy_a = cell(nq, 1);

% solution in the pipe.
for i = 1: nq
    z = zq(i);
    r = 0;
    mu = Mp.mu;
    [pz_a{i,1}, vz_a{i,1}, t_a, ~, ~] = pipe_crack_inf(Mp.L, Mp.R, Mc.w0, Mp.rho, Mp.c, Mc.c, g, dt_a, z, r, mu);
end

% solution in the crack.
for i = 1: nq
    z = 0;
    r = rq(i);
    mu = Mp.mu;
    [pxy_a{i,1}, ~, t_a, ~, ~] = pipe_crack_inf(Mp.L, Mp.R, Mc.w0, Mp.rho, Mp.c, Mc.c, g, dt_a, z, r, mu);
end

%% visualize the compare the solutions.
figure;
tn = [0:dt:nt*dt];
for i = 1: nq
    plot(t_a,real(pxy_a{i})+ (i-1)*0.01,'r-');
    hold on;
    plot(t_n,pxy_q(i,:)+(i-1)*0.01,'k-');
    hold on;
end
hold off;
xlim([1.5, T]);
%
figure;
for i = 1: nq
    plot(t_a,real(pz_a{i}),'r-');
    hold on;
    plot(t_n,pz_q(i,:),'k-');
    hold on;
end
hold off;
xlim([1.5, T]);
%
figure;
for i = 1: nq
    plot(t_a,real(-vz_a{i}),'r-');
    hold on;
    plot(t_n,vz_q(i,:),'k-');
    hold on;
end
hold off;
xlim([1.5, T]);

%% difference









