%coupled_simulation_frac3d

%% Model parameters.
% conduit parameters
% change to two section pipe model.
% 
source = '../source';
addpath(genpath(source));
%%
% baseline model parameters, two section pipe model.
para.L                  = 300; % m
para.upper_ratio  = 0.5; % [0, 1]
para.R                  = 5;% m

para.rho_upper = 1800;
para.rho_lower  = 2800;
para.c_upper     = 1000;
para.c_lower      = 1800;

para.mu      = 10; % Pa.s
para.Lx       = 1000; % m, size in down dip direction.
para.Ly       = 1000; % m, size along the strike
para.w0      = 4;% m   crack full width
para.PTA    = 1e5;%

para.xc_ratio = 0.5;
para.yc_ratio = 0.5;

para.strike  = 0;%
para.dip      = 0;%

para.target  = 0.033; %
para.xq       = 0; %
para.yq       = 1000; %
para.nx       = 32;%
para.nz       = 16/min(para.upper_ratio, 1 - para.upper_ratio);%
para.nz_f    = 16;%
para.nr       = 16;%
para.order  = 4;
para.g        = 10;
para.Xc      = 0;
para.Yc      = 0;
para.Zc      = 1000;

%% conduit and crack.

Mc.R = para.R;
Mc.L = para.L;
Mc.nz = para.nz;
Mc.nr = para.nr;
Mc.order = para.order;
Mc.S = pi*Mc.R^2;
Mc.g = para.g;
Mc.mu = para.mu;
Mc.interface_split = true;
Mc.with_exsolution = false;

% conduit parameters
dz = Mc.L/Mc.nz;
% conduit wave speed, density for the upper and lower section.
Mc.rho_upper = para.rho_upper;
Mc.rho_lower  = para.rho_lower;
Mc.c_upper     = para.c_upper;
Mc.c_lower     = para.c_lower;
Mc.L_upper    = para.L*para.upper_ratio;
Mc.L_lower     = Mc.L - Mc.L_upper;

Mc.split_index = round(Mc.L_lower/dz) + 1;
Mc.rho = [Mc.rho_lower*ones(Mc.split_index, 1); Mc.rho_upper*ones(Mc.nz+2 - Mc.split_index, 1)];
Mc.c    = [Mc.c_lower*ones(Mc.split_index, 1); Mc.c_upper*ones(Mc.nz+2 - Mc.split_index, 1)];
Mc.K   = Mc.rho.*(Mc.c).^2;

Mc.pT.A = para.PTA; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.

%% fracture parameters.
Mf.w0 = para.w0;
Mf.Lx   = para.Lx;
Mf.Ly   = para.Ly;
Mf.nx   = para.nx;
Mf.ny   = para.nx;
Mf.nz   = para.nz_f;
Mf.order = para.order;
Mf.interp_order = para.order;
Mf.xs = 0.5*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) 0;
Mf.xc = para.xc_ratio*Mf.Lx;% the coupling location to the conduit.
Mf.yc = para.yc_ratio*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho  = Mc.rho(1);
Mf.c     = Mc.c(1);
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = para.mu;

% solid parameters are fixed for now.
Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu   = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;

% crack centroid location, fixed!
Mf.Xc = para.Xc;
Mf.Yc = para.Yc;
Mf.Zc = para.Zc;
Mf.strike = para.strike;
Mf.dip = para.dip;
%% construct coupled models.
cond = conduit(Mc);
frac = frac3d_o(Mf);
Model = coupledModel(cond, frac);
%% time domain simulation.

CFL = 0.5;
skip = 50;
T = 200;
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
Hp = Model.get_Hp(para.xq,para.yq);
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
    
    % surface displacement
    Us(:,i) = Hp*Model.u;
    
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
            
%             xm = Model.frac.geom.m.x;
%             ym = Model.frac.geom.m.y;
%             xp  = Model.frac.geom.p.x;
%             
            %fracture
%             figure(1);
%             subplot(3,1,1);
%             pcolor(xm, ym, vx);
%             xlabel('x'); ylabel('y');
%             cmap;
%             caxis([-2.5e-1 2.5e-1]);
%             shading INTERP;
%             title(sprintf('t=%f',t));
%             xlim([0 Model.frac.M.L]);
%             
%             subplot(3,1,2);
%             plot(xm, vx(round(size(vx,1)/2), :)');
%             xlabel('x');
%             ylabel('v');
%             ylim([-2.5e-1 2.5e-1]);
%             xlim([0 Model.frac.M.L]);
%             
%             subplot(3,1,3);
%             plot(xp, px);
%             xlabel('x');
%             ylabel('p');
%             ylim([-10e4, 10e4]);
%             xlim([0 Model.frac.M.L]);
            
            % plot the surface displacement vertical component.
%             figure(2);
            plot([1:i]*dt, Us(3,1:i),'-k');
            xlim([0 T]);xlabel('time (s)');ylabel('Uz');
            ylim([-10e-5, 10e-5]);
            grid on;
            drawnow;
            
            %conduit.
%             figure(3);
%             subplot(1,3,1)
%             pcolor(rm, z, vz');
%             shading INTERP;
%             cmap;
%             caxis([-5e-2 5e-2])
%             xlabel('r (m)'), ylabel('z (m)')
%             
%             subplot(1,3,2)
%             plot(uz, z);
%             xlim([-5e-2, 5e-2])
%             ylim([0 Model.conduit.M.L])
%             xlabel('u'), ylabel('z (m)')
%             
%             subplot(1,3,3)
%             plot(pz, z);
%             xlim([-5e4,5e4]);
%             ylim([0 Model.conduit.M.L])
%             xlabel('p'), ylabel('z (m)')
            
%             subplot(1,4,4)
%             plot(nz, z)
%             xlim([-1e-5, 1e-5]);
%             ylim([0 Model.conduit.M.L])
%             xlabel('n'), ylabel('z (m)');
            drawnow;
        end
        
    end
end
toc;
t = [1:nt]*dt; 
save('Us_nx_16_nr_8_200s','Us','t');

%%
d1 = load('Us_nx_16_nr_8_200s');
d2 = load('Us_nx_8_nr_8_400s');

plot(d1.t, d1.Us(3,:),'-k');
hold on;
plot(d2.t, d2.Us(3,:),'-k');
xlim([0 200]);xlabel('time (s)');ylabel('Uz');
ylim([-10e-5, 10e-5]);







