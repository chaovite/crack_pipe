% crack2d subject to distributed external loading.
clear
source_dir = '../source';
addpath(genpath(source_dir));
%%

% fluid
Mf.rho = 1.8e3;
Mf.c    = 2.5e3;
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 0; % viscosity

%solid
Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;

% crack geometry and grid
Mf.w0 = 0.01;
Mf.L   = 50;
Mf.nx = 100;
Mf.ny = 5;
Mf.order = 2;
Mf.interp_order = 4;

% external load distribution in time and space
Mf.pT.Ts = 0.2;% Source period.
Mf.pT.f0  = 20;
Mf.pT.t0  = 2/Mf.pT.f0;
Mf.pT.sigma = @(t) sin(2*pi/Mf.pT.Ts*t);
Mf.G  = @(t) 2*pi/Mf.pT.Ts*cos(2*pi/Mf.pT.Ts*t); % source time function of external forcing. d_sigma/d_t.
Mf.pX.A = 1e3;% external load amplitude
Mf.pX.L = Mf.L/5;% period of ex_load in space.
Mf.Gx = @(x) Mf.pX.A*sin(2*pi*(x-Mf.L/2)/Mf.pX.L);%a sin function center at x=L/2

Mf.xs = 0; % location of point source (not used in this script but keep it for the code to run)
Mf.xc = 0.5*Mf.L;% the coupling location to the conduit. (not used in this script but keep it for the code to run.)

Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.15; % estimated ration of boundary layer.

Mf.Xc   = 0;
Mf.Yc   = 0;
Mf.Zc   = 1000 ;
Mf.strike = 0;
Mf.dip = 0;
Mf.Xc   = 0;
Mf.Yc   = 0;
Mf.Zc   = 1000 ;
Mf.strike = 0;
Mf.dip = 0;
Mf.Lz  = Mf.L;

%
Model = frac2d(Mf);

%%
CFL = 0.5;
skip = 10;
T = 5;
use_imex = false;
plot_simu = true;

% time stepping
hmin = min(Model.geom.p.hx);
cmax = Mf.c;
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

if ~ use_imex
    A = Model.Ae + Model.Ai;
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end

% external load source:
dim = Model.dimensions;
Fp_ex = zeros(sum(dim), 1);
Fp_ex(1: dim(1)) =  -Model.M.K_t*1/Model.M.w0*Model.M.Ks_inv*Model.M.Gx(Model.geom.p.x);

fun = @(u,t) A*u + Fp_ex*Model.M.G(t); %Model.Fp*Model.M.G(t)

tic
for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
    end
    
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i,skip) == 0
            %plot solution.
            sol.p = Model.field(Model.u,1);
            sol.v = Model.field(Model.u,2);
            sol_v=Model.geom.mm.grd(sol.v);
            X = Model.geom.mm.X;
            Y= Model.geom.mm.Y;
            
            subplot(5,1,1);%v map.
            pcolor(X, Y, sol_v);
            xlabel('x'); ylabel('y');
            cmap;
            caxis([-5e-3 5e-3]);
            shading INTERP;
            %     colorbar;
            title(sprintf('t=%f',t));
            
            subplot(5,1,2);%v at y = w/2
            plot(Model.geom.m.x, sol_v(round(Model.geom.mm.ny/2), :)');
            xlabel('x');
            ylabel('v');
            ylim([-5e-3 5e-3]);
            
            subplot(5,1,3);% p
            plot(Model.geom.p.x, sol.p);
            xlabel('x');
            ylabel('p');
            ylim([-1, 1]*Mf.pX.A);
            
            subplot(5,1,4);% sigma: ex_load
            sigma = Model.M.Gx(Model.geom.p.x)*Mf.pT.sigma(t+dt);
            plot(Model.geom.p.x, sigma);
            xlabel('x');
            ylabel('sigma');
            ylim([-1, 1]*Mf.pX.A);
            
            subplot(5,1,5);% opening
            opening = Model.M.Ks_inv*(sigma + sol.p);
            plot(Model.geom.p.x, opening);
            xlabel('x');
            ylabel('w');
            ylim([-2 2]*1e-6)

            drawnow;
        end
        
    end
end
toc;
