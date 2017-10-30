% 3D crack wave simulation
addpath('../source')
%%
clear
% fluid-filled fracture parameters.
Mf.w0 = 2;
Mf.Lx   = 1000;
Mf.nx = 100;
Mf.order = 4;
Mf.interp_order = 4;

Mf.xs = 0.5*Mf.Lx;

% point source.
Mf.source.pT.A = 1;
Mf.source.pT.f0 = 5;
Mf.source.pT.t0 = 2/Mf.source.pT.f0;
% Mf.source.G  = @(t) Mf.source.pT.A*ricker(t, 5, 2/5);
Mf.source.G  = @(t) Mf.source.pT.A*ricker(t, 5, 2/5);
% Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.

% boundary condition. pressure at x=0, zero velocity at x=L;
% Mf.bc.p0.pT.A = 1e5;
% Mf.bc.p0.pT.f0 = 2;
% Mf.bc.p0.pT.t0 = 2/Mf.bc.p0.pT.f0;
% Mf.bc.p0.G = @(t) Mf.bc.p0.pT.A*ricker(t, Mf.bc.p0.pT.f0, Mf.bc.p0.pT.t0);

Mf.bc.v0.G = @(t) 0;

% Mf.bc.vn.pT.A = 0;
% Mf.bc.vn.pT.f0 = 5;
% Mf.bc.vn.pT.t0 = 2/Mf.bc.vn.pT.f0;
Mf.bc.vn.G = @(t) 0;

% fluid and solid properties.
Mf.rho = 1.8e3;
Mf.c    = 2.5e3;
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 50;
Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
Mf.isrigid = false;
%% construct coupled models.
Model = frac1d(Mf);
test_stability = true;

if test_stability
  E = Model.energy_norm();
  [is_stable, is_energy_stable,eig_s,eig_es] = test_energy_stability(Model.A,E);
  is_stable
  is_energy_stable
end
%% time stepping:
CFL = 0.5;
skip = 5;
T = 10;
use_imex = true;
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
fun = @(u,t) A*u + Model.Fp.ps*Model.M.source.G(t) + Model.M.bc.v0.G(t)*Model.Fp.bc.v0...
                           + Model.M.bc.vn.G(t)*Model.Fp.bc.vn;

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
            sol.p = Model.field(Model.u, 1);
            sol.v = Model.field(Model.u, 2);
            x_p = Model.geom.p.x;
            x_v = Model.geom.v.x;
            
            figure(1);
            subplot(2,1,1);
            plot(x_p, sol.p);
            xlabel('x'); ylabel('p');
            title(sprintf('t=%f',t));
            ylim([-2 2]*Model.M.source.pT.A*(Model.M.rho*Model.M.c));
            xlim([0 Mf.Lx]);
            
            subplot(2,1,2);
            plot(x_v, sol.v);
            xlabel('x'); ylabel('v');
            title(sprintf('t=%f',t));
            xlim([0 Mf.Lx]);
            ylim([-2 2]*Model.M.source.pT.A);
            drawnow;
        end
        
    end
end
toc;







