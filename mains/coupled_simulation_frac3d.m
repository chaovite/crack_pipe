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

para.strike  = 90;%
para.dip      = 60;%

para.target  = 0.033; %
para.xq       = 0; %
para.yq       = 1000; %
para.nx       = 8;%
para.nz       = 8/min(para.upper_ratio, 1 - para.upper_ratio);%
para.nz_f    = 8;%
para.nr       = 8;%
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
Mc.order_r = para.order;
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
Mc.pT.t = 5; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T).^2); % external force from the conduit.

%% fracture parameters.
Mf.w0 = para.w0;
Mf.Lx   = para.Lx;
Mf.Ly   = para.Ly;
Mf.nx   = para.nx;
Mf.ny   = para.nx;
Mf.nz   = para.nz_f;
Mf.order = para.order;
Mf.order_z = para.order;
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
T = 100;
use_imex = true;
plot_simu = true;

% time stepping
hmin = min([Model.conduit.geom.dz, Model.frac.geom.p.hx, Model.frac.geom.p.hy]);
cmax = max(max(Model.conduit.M.c), Model.frac.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt);
%% fields to be stored.
nkeep   = floor(nt/skip);
time      = [skip:skip:nt]*dt;
% pzs       = cell(nkeep, 1);
% vzs       = cell(nkeep, 1);
% uzs       = cell(nkeep, 1);
% pxys     = cell(nkeep, 1);
% vx_ms  = cell(nkeep, 1);
% uxs       = cell(nkeep, 1);
Us         = zeros(3*length(para.xq), nt);
%%
if ~ use_imex
    A = Model.Ae + Model.Ai;
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end

fun = @(u,t) A*u + Model.Fp(:,1)*Model.conduit.M.G(t) + Model.Fp(:,2)*Model.frac.M.G(t);
tic
Hp = Model.get_Hp(para.xq, para.yq);

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
    
    if mod(i,skip) == 0
        iter = i/skip;
        
        %plot solution.
        % unknowns in the conduit.
        nr  = Model.conduit.geom.nr;
        nz = Model.conduit.geom.nz;
        vz  = reshape(Model.field(Model.u, [1, 1]),[nr,nz]); % [nr, nz]
        pz  = Model.field(Model.u, [1, 2]); % [nz, 1]
        h    = Model.field(Model.u, [1, 3]); % [1,  1]
        uz  = (Model.conduit.op.W1*vz)';%   [nz, 1]
        
        z   = Model.conduit.geom.z;
        rm  = Model.conduit.geom.rm;
        
        % unknowns in the crack:
        pxy      = Model.frac.geom.p.grd(Model.field(Model.u, [2, 1])); %[nxp, nxp]
        X_pxy  = Model.frac.geom.p.X;
        Y_pxy  = Model.frac.geom.p.Y;
        
        vx   = Model.field(Model.u, [2, 2]);
        ux   = Model.frac.geom.ux.grd(Model.frac.op.vx.Wz3*vx);
        X_ux = Model.frac.geom.ux.X;
        Y_ux = Model.frac.geom.ux.Y;
        
        vy   = Model.field(Model.u, [2, 3]);
        uy   = Model.frac.geom.uy.grd(Model.frac.op.vy.Wz3*vy);
        X_uy = Model.frac.geom.uy.X;
        Y_uy = Model.frac.geom.uy.Y;
        
        % reshape vx, vy and keep the middle slice.
        vx = Model.frac.geom.vx.grd(vx);
        vy = Model.frac.geom.vy.grd(vy);
        
        vx_ymid =  round(length(Model.frac.geom.vx.y)/2);
        vy_xmid =  round(length(Model.frac.geom.vy.x)/2);
        
        vx_m = squeeze(vx(:,vx_ymid,:)); %[nz, nx]
        [Z_vxm, X_vxm] = meshgrid(Model.frac.geom.vx.x, Model.frac.geom.vx.z); %[nz, nx]
        vy_m = squeeze(vy(:,:,vy_xmid)); %[nz, ny]
        [Z_vym, Y_vym] = meshgrid(Model.frac.geom.vy.y, Model.frac.geom.vy.z);%[nz, ny]
        
%         pzs{iter,1}       = pz;
%         vzs{iter,1}       = vz;
%         uzs{iter,1}       = uz;
%         pxys{iter,1}     = pxy;
%         vx_ms{iter,1}  = vx_ms;
%         uxs{iter,1}       = ux;
        
        %% save the data in d structure.
%         d.time = time;
% 
%         % fields in the pipe.
%         d.nr = nr;
%         d.nz = nz;
%         d.z   = z;
%         d.rm = rm;
% 
%         d.pzs = pzs;
%         d.vzs = vzs;
%         d.uzs = uzs;
% 
%         % fields in the crack.
%         % save pxys
%         d.X_pxy = X_pxy;
%         d.Y_pxy = Y_pxy;
%         d.pxys = pxys;
% 
%         % save uxs
%         d.X_ux = X_ux;
%         d.Y_ux = Y_ux;
%         d.uxs   = uxs;
% 
%         % save vxm
%         d.Z_vxm = Z_vxm;
%         d.X_vxm = X_vxm;
%         d.vx_ms = vx_ms;
% 
%         d.Mc = Mc;
%         d.Mf  = Mf;
%         save(['data',num2str(iter)],'d');
        
        if plot_simu
            % plot the crack.
            figure(1);
            subplot(4,1,1)
            pcolor(Z_vym, Y_vym, vy_m); cmap;shading interp
            caxis([-0.03, 0.03]);
            colorbar
            subplot(4,1,2)
            pcolor(Z_vxm, X_vxm, vx_m); cmap;shading interp
            caxis([-0.03, 0.03]);
            colorbar
            subplot(4,1,[3,4])
            pcolor(X_pxy, Y_pxy, pxy); cmap;shading interp
            caxis([-1e4, 1e4]);
            colorbar
            drawnow;
            
            % plot the conduit.
            figure(2);
            subplot(1,3,1)
            pcolor(rm, z, vz');
            shading INTERP;
            cmap;
            caxis([-5e-1 5e-1])
            colorbar
            xlabel('r (m)'), ylabel('z (m)')
            
            subplot(1,3,2)
            plot(uz, z);
            xlim([-5e-1, 5e-1])
            ylim([0 Model.conduit.M.L])
            xlabel('u'), ylabel('z (m)')
            
            subplot(1,3,3)
            plot(pz, z);
            xlim([-1e5,1e5]);
            ylim([0 Model.conduit.M.L])
            xlabel('p'), ylabel('z (m)')
            drawnow;
            
            % plot the surface displacement vertical component.
            figure(3);
            plot([1:i]*dt, Us(3,1:i),'-k');
            xlim([0 T]);xlabel('time (s)');ylabel('Uz');
            ylim([-10e-5, 10e-5]);
            grid on;
            drawnow;
        end
    end
end

%% save the data in d structure.
% d.time = time;
% 
% % fields in the pipe.
% d.nr = nr;
% d.nz = nz;
% d.z   = z;
% d.rm = rm;
% 
% d.pzs = pzs;
% d.vzs = vzs;
% d.uzs = uzs;
% 
% % fields in the crack.
% % save pxys
% d.X_pxy = X_pxy;
% d.Y_pxy = Y_pxy;
% d.pxys = pxys;
% 
% % save uxs
% d.X_ux = X_ux;
% d.Y_ux = Y_ux;
% d.uxs   = uxs;
% 
% % save vxm
% d.Z_vxm = Z_vxm;
% d.X_vxm = X_vxm;
% d.vx_ms = vx_ms;
% 
% d.Mc = Mc;
% d.Mf  = Mf;
% save('data','d');














