% test_compare_frac2dxy_frac3d inviscid
%
% test if the frac3d+conduit and frac2dxy+pipe1d produces the same
% solution in time domain under the limit of zero viscosity.
%
clear
source_dir = '../source';
addpath(genpath(source_dir));
%% frac3d
% fluid-filled fracture parameters.
Mf.w0  = 2;
Mf.Lx   = 2e3;
Mf.Ly   = 2e3;
Mf.nx   = 64;
Mf.ny   = 64;
Mf.nz   = 20;
Mf.order = 4;
Mf.order_z = Mf.order;
Mf.interp_order = 6;
Mf.xs = 0.5*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) ricker(t, 2, 2/2);
Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = true;
Mf.r_g     =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl    =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
% fluid and solid properties.
Mf.rho = 1800;
Mf.c    = 1200;
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
model = frac2dxy(Mf);
%
% ------------------------------------Time integration-------------------------------------
CFL  = 0.5;
T       = 5;
plot_simu = true;

% time stepping
[cmax, hmin] = model.getCFL();
dt = CFL*hmin/cmax;
nt = ceil(T/dt);
skip    =ceil(nt/20);
nkeep = floor(nt/skip);
pxys    = cell(nkeep, 1);
tkeep  = zeros(nkeep, 1);

A = model.getA();
fun = @(u, t) A*u + model.Fp*model.M.G(t);
model=model.init();
u = model.u;
cnt = 0;

for i=1:nt
    t = (i-1)*dt;
    u = lsrk4(fun,u,t,dt);
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i, skip) == 0
            cnt = cnt+1;
            pxy         = model.grd.p.grd(u(model.indu.indp));
            
            X            = model.grd.p.X;
            Y            = model.grd.p.Y;
            pxys{cnt, 1}= pxy;
            tkeep(cnt) = t+dt;
            
            pcolor(X,Y, pxy);
            shading interp;
            cmap;
            colorbar;
            caxis([-5000, 5000]);
            %             axis off
            title('Bottom Crack');
            drawnow();
        end
    end
end

d1.name = 'frac2dxy';
d1.isrigid = false;
d1.pxy    = pxys;
d1.X       = X;
d1.Y       = Y;
d1.n       = nkeep;
%%  frac3d_o.
clear model;
model = frac3d_o(Mf);

%  ------------------------------------Time integration-------------------------------------
CFL  = 0.5;
T       = 5;
plot_simu = true;

% time stepping
hmin = min(model.geom.p.hx, model.geom.p.hy);
cmax = Mf.c;
dt = CFL*hmin/cmax;

nt = ceil(T/dt);
skip    =ceil(nt/20);
nkeep = floor(nt/skip);
pxys    = cell(nkeep, 1);
tkeep  = zeros(nkeep, 1);

fun = @(u, t) model.A*u + model.Fp*model.M.G(t);
u = zeros(sum(model.dimensions()),1);
cnt = 0;
for i=1:nt
    t = (i-1)*dt;
    u = lsrk4(fun,u,t,dt);
    
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i, skip) == 0
            cnt = cnt+1;
            pxy         = model.geom.p.grd(model.field(u,1));
            X            = model.geom.p.X;
            Y            = model.geom.p.Y;
            pxys{cnt, 1}= pxy;
            tkeep(cnt) = t+dt;
            
            pcolor(X,Y, pxy);
            shading interp;
            cmap;
            colorbar;
            caxis([-5000, 5000]);
            %             axis off
            title('Bottom Crack');
            drawnow();
        end
    end
end

d2.name = 'frac3d_o';
d2.isrigid = false;
d2.pxy    = pxys;
d2.X       = X;
d2.Y       = Y;
d2.n       = nkeep;

save('test_compare_frac2dxy_frac3d_inviscid_rigid','d1','d2','nkeep');
%%
load('test_compare_frac2dxy_frac3d_inviscid_rigid','d1','d2','nkeep');
for i = 1: nkeep
    subplot(1,3,1);
    pcolor(d1.X,d1.Y, d1.pxy{i});
    shading interp;
    cmap;
    colorbar;
    caxis([-5000, 5000]);
    title('frac2dxy');
    subplot(1,3,2);
    pcolor(d2.X,d2.Y, d2.pxy{i});
    shading interp;
    cmap;
    colorbar;
    caxis([-5000, 5000]);
    title('frac 3d');
    
    subplot(1,3,3);
    pcolor(d2.X,d2.Y, d2.pxy{i}-d1.pxy{i});
    shading interp;
    cmap;
    colorbar;
%     caxis([-1e-6, 1e-6]);
    title('frac2dxy-frac3d difference');
   drawnow
end

