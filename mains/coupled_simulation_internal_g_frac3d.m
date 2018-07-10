%coupled_simulation_frac3d

%% Model parameters.
% conduit parameters
% change to two section pipe model.
%
source = '../source';
addpath(genpath(source));
%%

sourcedir = '../source';
addpath(genpath(sourcedir));

Mc.R  = 5; 
Mc.L   = 300; 
Mc.nz = 20; 
Mc.nr = 16; 
Mc.order = 6;
Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1); 
Mc.g = 10;  
Mc.with_exsolution=false;
Mc.interface_split=false;
Mc.mu = 50*ones(Mc.nz+1, 1);
z = Mc.L/Mc.nz*[0:Mc.nz]';

% density profile
rho0     = 800;
rho1     = 2000;
alpha    = (log(rho1)-log(rho0))/Mc.L;

Mc.rho = rho0*exp(alpha*(Mc.L-z));
Mc.c    = 1000*ones(Mc.nz+1,1);
Mc.K    = Mc.rho.*Mc.c.^2;

Mc.Mg = alpha - Mc.rho*Mc.g./Mc.K;

%%
Mc.epsilon = 0;
Mc.pT.A = 2e5; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t  = 5; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
model = conduit_internal_g(Mc);
%% fracture parameters.
Mf.w0 =  4;
Mf.Lx   = 2000;
Mf.Ly   = 2000;
Mf.nx   = 32;
Mf.ny   = 32;
Mf.nz   = 16;
Mf.order = 6;
Mf.interp_order = 6;
Mf.xs = 0.5*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) 0;
Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho  = Mc.rho(1);
Mf.c     = Mc.c(1);
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = Mc.mu(1);

% solid parameters are fixed for now.
Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu   = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;

% crack centroid location, fixed!
Mf.Xc = 0;
Mf.Yc = 500;
Mf.Zc = 1000;
Mf.strike = 0;
Mf.dip     = 0;
%% construct coupled models.
cond = conduit_internal_g(Mc);
frac = frac3d_o(Mf);
Model = coupledModel(cond, frac);
%% time domain simulation.
CFL = 0.5;
skip = 50;
T = 2;
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
        
        % save the data in d structure.
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
        time_i = i*dt;
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
            suptitle(sprintf('time = %f sec',time_i));
            drawnow;
            
            % plot the conduit.
            figure(2);
            subplot(1,3,1)
            pcolor(rm, z, vz');
            shading INTERP;
            cmap;
            caxis([-1 1]);
            colorbar
            xlabel('r (m)'), ylabel('z (m)')
            
            subplot(1,3,2)
            plot(uz, z);
            xlim([-1, 1])
            ylim([0 Model.conduit.M.L])
            xlabel('u'), ylabel('z (m)')
            
            subplot(1,3,3)
            plot(pz, z);
            xlim([-2.5e5,2.5e5]);
            ylim([0 Model.conduit.M.L])
            xlabel('p'), ylabel('z (m)')
            drawnow;
            
            % plot the surface displacement vertical component.
            figure(3);
            plot([1:i]*dt, Us(3,1:i),'-k');
            xlim([0 T]);xlabel('time (s)');ylabel('Uz');
            ylim([-1e-4, 1e-4]);
            grid on;
            drawnow;
        end
    end
end