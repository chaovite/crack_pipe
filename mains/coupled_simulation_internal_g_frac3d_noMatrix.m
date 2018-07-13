% coupled conduit and bottom crack wave simulation

% use absolute path if working on the cluster.
source_dir = '../source';
addpath(genpath(source_dir));
%%
Mc.R  = 5;
Mc.L   = 300;
Mc.nz = 20;
Mc.nr = 30;
Mc.order = 6;
Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1);
Mc.g = 10;
Mc.with_exsolution=false;
Mc.interface_split=false;
Mc.mu = 50*ones(Mc.nz+1, 1);
z = Mc.L/Mc.nz*[0:Mc.nz]';
Mc.variables = {'vz','pz','h','hL'};

% density profile
rho0     = 800;
rho1     = 2000;
alpha    = (log(rho1)-log(rho0))/Mc.L;

Mc.rho = rho0*exp(alpha*(Mc.L-z));
Mc.c    = 1000*ones(Mc.nz+1,1);
Mc.K    = Mc.rho.*Mc.c.^2;
Mc.Mg = alpha - Mc.rho*Mc.g./Mc.K;
%%
Mc.epsilon  = 0;
Mc.pT.A     = 2e5; % pressure perturbation amplitude
Mc.pT.T     = 0.25; % pressure perturbation duration
Mc.pT.t      = 5; % pressure perturbation center time
Mc.G         = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
%% fracture parameters.
Mf.w0   =  4;
Mf.Lx   = 2000;
Mf.Ly   = 2000;
Mf.nx   = 50;
Mf.ny   = 50;
Mf.nz   = 50;
Mf.order = 6;
Mf.interp_order = 6;
Mf.xs = 0.5*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) 0;
Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.
Mf.variables = {'pxy','vx','vy'};

Mf.isrigid = false;
Mf.r_g   =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl  =  0.1; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho  = Mc.rho(1);
Mf.c      = Mc.c(1);
Mf.K     = Mf.rho*Mf.c^2;
Mf.mu  = Mc.mu(1);

% solid parameters are fixed for now.
Mf.cp     = 5e3;
Mf.cs     = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs    = Mf.cs^2*Mf.rhos;
Mf.nu    = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;

% crack centroid location, fixed!
Mf.Xc      = 0;
Mf.Yc      = 500;
Mf.Zc      = 1000;
Mf.strike = 0;
Mf.dip     = 0;
x_obs     = 0;
y_obs     = 1000;
%% construct coupled models.
Model = coupledModel_NoMatrix(Mc, Mf);
%% time stepping:
CFL = 0.5;
% skip = 40;
T = 10;
use_imex = true;
plot_simu = false;

% time stepping
[cmax, hmin] = Model.getCFL();
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

% number of frames to keep.
nkeep  = 5;
skip     = floor(nt/nkeep);

% put the absolute path of the savefolder if working on the cluster.
save_folder = '/Users/Chaovite/Documents/data/time_simulation3d/time_simu_internal_g'; 

if ~ use_imex
    fun = @(u, t) Model.ft(u,t) + Model.Ai*u;
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
    fun = @(u,t) Model.ft(u,t);
end

disp('Time integration starts here:');
Us = zeros(3, nt);
itr = 0;
ta=tic;

for i=1:nt
    t = (i-1)*dt;
    
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt, fun, Model.Ai,L,U,p,q));
    end
    
    if mod(i, round(nt/100))==0
        t_elapse = toc(ta);
        fprintf( '%% %f  finished, %8.5f sec / step \n', round( i*100/nt), t_elapse/round(nt/100));
        ta=tic;
    end
    
    Us(:,i) = Model.frac.disloc3d.eval_disp_p(Model.field(Model.u, [2 1]), x_obs, y_obs);
    
    t_i = i*dt;
    
    % save the data for later plotting.
    
    if mod(i, skip) == 0
        itr = itr+1;
        % extract solution.
        [vz, pz, ~, ~, p_mat, vx_mat, ~] = Model.fields(Model.u);
        uz = Model.cond.op.W2*Model.field(Model.u, [1, 1]);
        rm = Model.cond.geom.rm;
        z = Model.cond.geom.z;
        Xp = Model.frac.geom.p.X;
        Yp = Model.frac.geom.p.Y;
        
        v_slice = squeeze(vx_mat(:, round(Model.frac.geom.vx.ny/2), :));
        [X_vx_slice, Z_vx_slice] = meshgrid(Model.frac.geom.vx.x, Model.frac.geom.vx.z);
        px_slice        = squeeze(p_mat(round(Model.frac.geom.p.ny/2), :));
        x_px_slice    = Model.frac.geom.p.x;
        
        % save solution into a structure.
        d=struct();
        d.Mf  = Mf;
        d.Mc = Mc;
        d.time = t_i;
        d.z = z;
        d.rm = rm;
        d.vz = vz;
        d.pz = pz;
        d.uz = uz;
        
        d.Xp = Xp;
        d.Yp = Yp;
        d.pxy = p_mat;
        
        d.X_vx_slice = X_vx_slice;
        d.Z_vx_slice = Z_vx_slice;
        d.vx_slice = v_slice;
        
        d.x_px_slice=x_px_slice;
        d.px_slice =px_slice;
        
        d.Lx = Model.frac.M.Lx;
        d.Ly = Model.frac.M.Ly;
        d.w0 = Model.frac.M.w0;
        d.L   = Model.cond.M.L;
        d.R  = Model.cond.M.R;
        save_path = sprintf('%s/data_frame_%d.mat', save_folder, itr);
        save(save_path, 'd');
        
        % plot solution.
        if plot_simu
            h1=figure(1);
            pcolor(Xp, Yp, p_mat);
            xlabel('x'); ylabel('y');
            cmap;
            caxis([-1e4, 1e4]);
            colorbar('location','eastoutside');
            shading INTERP;
            title(sprintf('t=%f',t));
            axis([0 Model.frac.M.Lx, 0 Model.frac.M.Ly]);
            daspect([1 1 1]);
            set(gcf,'PaperPositionMode','auto');
            set(gca, 'fontsize',18)
            
            h2 = figure(2);
            % plot the surface displacement vertical component.
            plot([1:i]*dt, Us(3,1:i),'-k');
            xlim([0 T]);xlabel('time (s)');ylabel('Uz');
            ylim([-1e-4, 1e-4]);
            set(gcf,'PaperPositionMode','auto');
            set(gca, 'fontsize',18);
            
            h3 = figure(3); %vx
            
            subplot(2,1,1);
            pcolor(X_vx_slice, Z_vx_slice, v_slice);
            xlabel('x'); ylabel('z');
            cmap;
            caxis([-3e-2 3e-2]);
            shading INTERP;
            title(sprintf('t=%f',t));
            set(gca,'fontsize',16);
            axis([0 Model.frac.M.Lx, 0 Model.frac.M.w0]);
            colorbar('location','southoutside');
            
            subplot(2,1,2);
            plot(Model.frac.geom.p.x, px_slice);
            xlabel('x'); ylabel('p');
            xlim([0 Model.frac.M.Lx]);
            ylim([-1e4, 1e4]);
            set(gca,'fontsize',16);
            set(gcf,'PaperPositionMode','auto');
            
            %conduit.
            h4 = figure(4);
            subplot(1,2,1);
            plot(uz, z);
            xlim([-1 1])
            ylim([0 Model.cond.M.L])
            xlabel('u'), ylabel('z (m)');
            
            subplot(1,2,2)
            plot(pz, z);
            xlim([-2.5e5,2.5e5]);
            ylim([0 Model.cond.M.L])
            xlabel('p'), ylabel('z (m)');
            suptitle(sprintf('t=%f',t));
            set(gcf,'PaperPositionMode','auto');
            drawnow
            
        end
    end
end

% save the time series of surface displacements
t_Us = [0:nt]*dt;
save_path = sprintf('%s/data_Us.mat', save_folder);
save(save_path, 'Us','t_Us');