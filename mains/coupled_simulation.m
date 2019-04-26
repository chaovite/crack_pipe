% coupled conduit and bottom crack wave simulation
clear
source_dir = '../source';
addpath(genpath(source_dir));
%%
% conduit parameters
Mc.R = 20;
Mc.L = 1000;
Mc.nz = 100;
Mc.nr = 50;
Mc.order = 6;
Mc.order_r = 6;
Mc.S = pi*Mc.R^2;
Mc.g = 10;
Mc.mu = 50;
dz = Mc.L/Mc.nz;
z = dz*[0: Mc.nz]';
Mc.tau = 1;

% conduit wave speed, density for the upper and lower section.
Mc.rho_upper = 1200;
Mc.rho_lower  = 2800;
Mc.c_upper     = 800;
Mc.c_lower     = 1800;
Mc.L_upper    = 0.5*Mc.L;
Mc.L_lower    = Mc.L - Mc.L_upper;

Mc.interface_split = true;
Mc.with_exsolution = false;

Mc.split_index = round(Mc.L_lower/dz) + 1;
Mc.rho = [Mc.rho_lower*ones(Mc.split_index, 1); Mc.rho_upper*ones(Mc.nz+2 - Mc.split_index, 1)];
Mc.c    = [Mc.c_lower*ones(Mc.split_index, 1); Mc.c_upper*ones(Mc.nz+2 - Mc.split_index, 1)];
Mc.K   = Mc.rho.*(Mc.c).^2;

% [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
if Mc.interface_split && Mc.with_exsolution
    [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
    Mc.split_index = find(Mc.ex, 1) - 1;
    z = Mc.L/Mc.nz*[0:(Mc.split_index -1), Mc.split_index, Mc.split_index:Mc.nz]';
    [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
end

Mc.pT.A = 1e4; % pressure perturbation amplitude
Mc.pT.T = 0.5; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit (external pressure)

% fluid-filled fracture parameters.
Mf.w0 = 1.5;
Mf.L   = 1500;
Mf.nx = 50;
Mf.ny = 30;
Mf.order = 6;
Mf.interp_order = 4;
Mf.xs = 0*Mf.L;
% Mf.G  = @(t) 5e4/(Mc.c(1)*Mc.rho(1))*ricker(t, 1, 2); % no external forcing in the crack.
Mf.G  = @(t) 0; % external forcing in the crack (point source injection)
Mf.xc = 0.5*Mf.L;% the coupling location to the conduit.

Mf.lc  = 50;%coupling length scale.
Mf.Lz = Mf.L; % length of crack in the out-of-plane direction. 2D plane strain crack.
Mf.Xc = 0;
Mf.Yc = 0;
Mf.Zc = 1000;
Mf.strike = 0;
Mf.dip = 0;

% the cross-section surface area of the crack coupled to the conduit.
Mf.S = Mf.lc * Mf.w0; 
Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.15; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho = Mc.rho(1);
Mf.c    = Mc.c(1);
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 20;
Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;

% observation points.
x_obs = 0;
y_obs = 1000;
%% construct coupled models.
frac = frac2d(Mf); % 2d crack with 1 length dimension, 1 width dimension, resolve viscous boundary layer.
cond = conduit(Mc); % 2d radially symmetric conduit.
%%
Model = coupledModel(cond, frac);
%% time stepping:
CFL = 0.5;
skip = 50;
T = 400;
use_imex = true;
plot_simu = true;

% time stepping
hmin = min(Model.conduit.geom.dz, Model.frac.geom.p.hx);
cmax = max(max(Model.conduit.M.c), Model.frac.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

% surface displacement.
Us = zeros(3, nt);

if ~ use_imex
    A = Model.Ae + Model.Ai;
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end
fun = @(u,t) A*u + Model.Fp(:,1)*Model.conduit.M.G(t) + Model.Fp(:,2)*Model.frac.M.G(t);

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
    
    Us(:,i) = Model.frac.disloc2d.eval_disp_p(Model.field(Model.u, [2 1]), x_obs, y_obs);
    
    if plot_simu
        if mod(i,skip) == 0
            %plot solution.
            if Mc.with_exsolution
                [vz, pz, ~,~, px, vx] = Model.fields(Model.u);
            else
                [vz, pz,~, px, vx] = Model.fields(Model.u);
            end
            uz = Model.conduit.op.W2*Model.field(Model.u, [1, 1]);
            
            rm = Model.conduit.geom.rm;
            z = Model.conduit.geom.z;
            
            xm = Model.frac.geom.m.x;
            ym = Model.frac.geom.m.y;
            xp  = Model.frac.geom.p.x;
            
            %fracture
            figure(1);
            subplot(3,1,1);
            pcolor(xm, ym, vx);
            xlabel('x'); ylabel('y');
            cmap;
            caxis([-2.5e-1 2.5e-1]);
            shading INTERP;
            title(sprintf('t=%f',t));
            xlim([0 Model.frac.M.L]);
            
            subplot(3,1,2);
            plot(xm, vx(round(size(vx,1)/2), :)');
            xlabel('x');
            ylabel('v');
            ylim([-2.5e-1 2.5e-1]);
            xlim([0 Model.frac.M.L]);
            
            subplot(3,1,3);
            plot(xp, px);
            xlabel('x');
            ylabel('p');
            ylim([-10e4, 10e4]);
            xlim([0 Model.frac.M.L]);
            
            % plot the surface displacement vertical component.
            figure(2);
            plot([1:i]*dt, Us(3,1:i),'-k');
            xlim([0 T]);xlabel('time (s)');ylabel('Uz');
            ylim([-10e-5, 10e-5]);
            grid on;
            drawnow;
            
%             conduit.
            figure(3);
            subplot(1,3,1)
            pcolor(rm, z, vz');
            shading INTERP;
            cmap;
            caxis([-5e-2 5e-2])
            xlabel('r (m)'), ylabel('z (m)')
            
            subplot(1,3,2)
            plot(uz, z);
            xlim([-5e-2, 5e-2])
            ylim([0 Model.conduit.M.L])
            xlabel('u'), ylabel('z (m)')
            
            subplot(1,3,3)
            plot(pz, z);
            xlim([-5e4,5e4]);
            ylim([0 Model.conduit.M.L])
            xlabel('p'), ylabel('z (m)')
            drawnow;
            
        end
        
    end
end
toc;







