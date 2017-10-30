source_dir = '../source';
addpath(genpath(source_dir));
%% Simulation settings
clear; close all;
skip            = 10;
order           = 6;
show_grid       = false;
plot_source = false;

nx              = 100;
ny              = 50;
use_imex       = true;
test_stability  = true;
T               = 400;

% Fluid properties
% TODO: w0 is constant for now.
L   = 1e3;
w0  = 1;
% grid stretching parameters in y direction.
r_g = 0.3;% ratio of grid points in boundary layer.
r_bl = 0.3;% ratio of boundary layer thickness to fracture width.

rho = 1e3;
c    = 1e3;
K   = rho*c^2;
mu  = 1e-1*1e3;
cp    = 1.5e3;
cs    = 1e3;
rhos = 2e3;
G   = cs^2*rhos;
nu  = ((cp/cs)^2-2)/((cp/cs)^2-1)/2;

% rho = 2.5;
% K   = 2.25;
% mu  = 1e-7;
%  G   = 10;
% nu  = 0.25;

% Source term
xs = 0;
f0 = 1; t0 = 2/f0;
interp_order = 2;
% g = @(t) exp(-a^2*(t-t0).^2/2);% the unit of the source is the same as the velocity.
g = @(t) ricker(t, f0, t0);% the unit of the source is the same as the velocity.
% if plot source:
if plot_source
    dts = 1/(5*f0);
    ts = 0: dts:100;
    gs = g(ts);
    figure;
    plot(ts, gs);
    xlabel('t');
    ylabel('g');
    figure;
    [ghats,fs] = fft_dim(gs,dts);
    plot(fs, abs(ghats));
    xlabel('f (Hz)');
    ylabel('abs(ghat)');
end

% Initialize rectangular grid and construct staggered grid operators
grid_types = {'pp','mm','pm','mp','p','m'};  

% Remove the end points of the m grid
truncate = true;
for i=1:length(grid_types)
  gt = grid_types{i};
  geom.(gt)       = stretched_grid(gt,order,'strong',nx,ny,L,w0,r_g,r_bl,truncate);
%   geom.(gt)       = grids(gt, nx,ny,operator_type,L,w0);
  op.(gt)         = operators(gt,order,geom.(gt),'strong');
end
metrics = struct();
f = fluid(geom,metrics,op,w0,rho, K,mu, G, nu);
e = f.source(xs,interp_order);
[f.Ae f.Ai] = f.interior();

% Test stability
if test_stability && nx*ny <= 500
  E = f.energy_norm();
  A = f.Ae + f.Ai;
  [is_stable, is_energy_stable,eig_s,eig_es] = test_energy_stability(A,E);
  is_stable
  is_energy_stable
end
%%
% c = sqrt(K/rho);
dt = 0.5*op.p.hx/c;
nt = ceil(T/dt);

if ~ use_imex
  f.A = f.Ae + f.Ai;
else
  [L,U,p,q,B] = imex_ark4_get_lu(f.Ai,dt);
  f.A = f.Ae;
end
fun = @(u,t) f.A*u + e*g(t);

tic
for i=1:nt
   t = (i-1)*dt;
   if ~use_imex
    f.u = lsrk4(fun,f.u,t,dt);
   else
    f.u = imex_ark4_lu(f.u,t,dt,fun,f.Ai,L,U,p,q); 
   end
   
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
% 
   if mod(skip,10) == 0
    sol.p = f.field(f.u,1);
    sol.v = f.field(f.u,2);
    sol_v = geom.mm.grd(sol.v);
    % Pressure
    %plot(geom.p.x,sol.p,'o-');
%     ylim([-1 1]);
    % Velocity
    subplot(3,1,1);
    pcolor(geom.mm.X,geom.mm.Y,geom.mm.grd(sol.v));   
    xlabel('x'); ylabel('y');
    cmap;
    caxis([-0.5 0.5]);
%     ylim([0 1e-2]);
    shading INTERP;
%     colorbar;
    title(sprintf('t=%f',t));
   
    subplot(3,1,2);
    plot(geom.m.x, sol_v(round(ny/2), :)');
    xlabel('x');
    ylabel('v');
    ylim([-0.5 0.5]);
    subplot(3,1,3);
    plot(geom.p.x, sol.p);
    xlabel('x');
    ylabel('p');
    ylim([-0.5 0.5]*rho*c);
    
%     sol_v = geom.mm.grd(sol.v);
    
%     plot(geom.mm.x,sol_v(10,:), geom.mm.x, zeros(size(geom.mm.x)),'o');
%     ylim([-0.1, 0.1]);
%     plot(geom.p.x, f.Ks_inv*sol.p);
    drawnow;
    if t >=2-dt/2 && t<=2+dt/2
        pause();
    end
    if t >=4-dt/2 && t<=4+dt/2
        pause();
    end
    
   end
end
toc;


