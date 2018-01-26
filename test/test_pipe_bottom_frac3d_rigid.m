% test_pipe_bottom_frac3d_rigid
%% test the coupling of pipe and bottom crack use analytical solution.
% crack is radial symmetric, inifinite and rigid in the analytical solution.
% gravity in the turn-off.
clear;
source_dir = '../source';
addpath(genpath(source_dir));
%% frac3d+conduit solution
% conduit parameters
R = 1;
L = 2000;
rho = 1000;
cp  = 2000;
cc  = 250;
w0 = 1e-1;
mu = 0;
Lx = 2000;
Ly = 2000;

nz = 128;
nx = 256;
ny = 256;
order = 4;

Mp.R   = R;
Mp.L   = L;
Mp.nz = nz;
Mp.nr  = 8;
Mp.order = order;
Mp.S = pi*Mp.R^2;
Mp.g = 0;
Mp.mu = mu;
Mp.interface_split = false;
Mp.with_exsolution = false;
Mp.rho = rho * ones(Mp.nz+1, 1);
Mp.c    = cp * ones(Mp.nz+1, 1);
Mp.K    = Mp.rho.*Mp.c.^2;

Mp.pT.A = 1; % pressure perturbation amplitude
Mp.pT.T = 0.1; % pressure perturbation duration
Mp.pT.t = 2; % pressure perturbation center time
Mp.G = @(t) Mp.pT.A*exp(-0.5*((t-Mp.pT.t)/Mp.pT.T).^2); % external force from the conduit.

% fluid-filled fracture parameters.
Mc.w0 = w0;
Mc.Lx   = Lx;
Mc.Ly   = Ly;
Mc.nx   = nx;
Mc.ny   = ny;
Mc.nz   = 8;
Mc.order = order;
Mc.interp_order = 4;
Mc.xs = 0.75*Mc.Lx;
Mc.ys = 0.5*Mc.Ly;
Mc.G  = @(t) 0;
Mc.xc = 0.5*Mc.Lx;% the coupling location to the conduit.
Mc.yc = 0.5*Mc.Ly;% the coupling location to the conduit.

Mc.isrigid = true;
Mc.r_g  =  0.3;% ratio of grid points in boundary layer.
Mc.r_bl =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
Mc.rho = rho;
Mc.c    =  cc;
Mc.K    = Mc.rho*Mc.c^2;
Mc.mu = 0;

Mc.cp    = 5e3;
Mc.cs    = 2.7e3;
Mc.rhos = 3e3;
Mc.Gs   = Mc.cs^2*Mc.rhos;
Mc.nu   = ((Mc.cp/Mc.cs)^2-2)/((Mc.cp/Mc.cs)^2-1)/2;
Mc.Xc   = 0;
Mc.Yc   = 0;
Mc.Zc   = 1000 ;
Mc.strike = 0;
Mc.dip = 0;

x_obs = 0;
y_obs = 1000;
% construct coupled models.
cond = conduit(Mp);
frac  = frac3d_o(Mc);
model=coupledModel(cond, frac);
%%
% query points.
nq = 10;
% points in the pipe:
zq_idx = round(linspace(1,Mp.nz, nq))';
% analytica solution position downward with z=0 at the bottom.
% numerical solution position upward with z=0 at the bottom.
zq       = - full(model.conduit.geom.z(zq_idx));

% query points in the crack:
xq_idx = round(linspace(Mc.nx/2 + 0.1 *  Mc.nx/2, Mc.nx/2 + 1 + 0.75 *  Mc.nx/2, nq))';
xq        = model.frac.geom.p.x(xq_idx)' - Mc.Lx/2;
yq_idx = (Mc.nx/2 + 1)*ones(nq, 1);
yq       =  model.frac.geom.p.y(yq_idx)' - Mc.Ly/2;
rq        = sqrt(xq.^2 + yq.^2);
rq_idx = sub2ind(size(model.frac.geom.p.X), yq_idx, xq_idx);
%%
% ------------------------------------Time integration-------------------------------------
CFL = 0.5;
skip = 20;
T     =  L/cp +  0.5*Mc.Lx/cc + Mp.pT.t;
use_imex = false;
plot_simu = true;

% time stepping
hmin = min(model.conduit.geom.dz, model.frac.geom.p.hx);
cmax = max(max(model.conduit.M.c), model.frac.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

pz_q   = zeros(nq, nt+1); 
vz_q   = zeros(nq, nt+1); 
pxy_q = zeros(nq, nt+1);
t_n     = [0:dt:nt*dt];

model=model.update(zeros(sum(model.dimensions()), 1));

if ~ use_imex
    A = model.Ae + model.Ai;
else
    [L,U,p,q,B] = imex_ark4_get_lu(model.Ai,dt);
    A = model.Ae;
end
fun = @(u,t) A*u + model.Fp(:,1)*model.conduit.M.G(t) + model.Fp(:,2)*model.frac.M.G(t);
tic

z = model.conduit.geom.z;
X = model.frac.geom.p.X;
Y = model.frac.geom.p.Y;

for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        model=model.update(lsrk4(fun,model.u,t,dt));
    else
        model=model.update(imex_ark4_lu(model.u,t,dt,fun, model.Ai,L,U,p,q));
    end
    
    %plot solution.
    if isfield(Mc, 'with_exsolution') &&  Mc.with_exsolution
        [vz, pz, ~,~, pxy, vxy] = model.fields(model.u);
    else
        [vz, pz,~, pxy, vx, vy] = model.fields(model.u);
    end
    uz = model.conduit.op.W2*model.field(model.u, [1, 1]);
    
    % record the solutions at the query points.
    pz_q(:, i+1)   = pz(zq_idx);
    vz_q(:, i+1)   = uz(zq_idx); 
    pxy_q(:, i+1) = pxy(rq_idx);
    
    if plot_simu
        if mod(i,skip) == 0
            figure(1)
            subplot(2,1,1)
            plot(z, pz,'k-');
            hold on;
            plot(-zq, zeros(nq, 1),'k*');
            hold off;
            ylim(2*[-Mp.pT.A, Mp.pT.A]);
            xlim([0, Mp.L])
            title(sprintf('t = %8.2f', t));
            subplot(2,1,2)
            plot(z, uz,'k-');
            hold on;
            plot(-zq, zeros(nq, 1),'k*');
            hold off;
            ylim(2*[-Mp.pT.A, Mp.pT.A]/(cp*rho));
            xlim([0, Mp.L])
            title(sprintf('t = %8.2f', t));
            drawnow;
            
            figure(2)
            pcolor(X, Y, pxy);
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
t_a   = [0:dt_a:1000];
g  = Mp.G(t_a);

pz_a   = cell(nq, 1);
vz_a   = cell(nq, 1);
pxy_a = cell(nq, 1);

% solution in the pipe.
for i = 1: nq
    z = zq(i);
    r = 0;
    mu = Mp.mu;
    [pz_a{i,1}, vz_a{i,1}, t_a, ~, ~] = pipe_crack_inf(Mp.L, Mp.R, Mc.w0, rho, cp, cc, g, dt_a, z, r, mu);
end

% solution in the crack.
for i = 1: nq
    z = 0;
    r = rq(i);
    mu = Mp.mu;
    [pxy_a{i,1}, ~, t_a, ~, ~] = pipe_crack_inf(Mp.L, Mp.R, Mc.w0, rho, cp, cc, g, dt_a, z, r, mu);
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
xlim([2, T]);
%
figure;
for i = 1: nq
    plot(t_a,real(pz_a{i})- (i-1)*0.1,'r-');
    hold on;
    plot(t_n,pz_q(i,:)-(i-1)*0.1,'k-');
    hold on;
end
hold off;
xlim([2, T]);
%%
figure;
for i = 1: nq
    plot(t_a,real(-vz_a{i})+ (i-1)*0.1/(rho*cp),'r-');
    hold on;
    plot(t_n,vz_q(i,:)+(i-1)*0.1/(rho*cp),'k-');
    hold on;
end
hold off;
xlim([2, T]);








