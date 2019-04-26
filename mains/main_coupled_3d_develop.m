% coupled conduit and bottom crack wave simulation

% This script is used do some development, can be messy.
clear
source_dir = '../source';
addpath(genpath(source_dir));
%%
clear
% conduit parameters
Mc.R = 10;
Mc.L = 1000;
Mc.nz = 40;
Mc.nr = 20;
Mc.order = 2;
Mc.order_r = 2;
Mc.S = pi*Mc.R^2;
Mc.g = 10;
Mc.mu = 0;
Mc.interface_split = false;
z = Mc.L/Mc.nz*[0:Mc.nz]';
Mc.tau = 1;
% [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
% 
% if Mc.interface_split
%     Mc.split_index = find(Mc.ex, 1) - 1;
%     z = Mc.L/Mc.nz*[0:(Mc.split_index -1),Mc.split_index, Mc.split_index:Mc.nz]';
%     [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
% end

Mc.with_exsolution = false;
Mc.rho = 1200 * ones(Mc.nz+1, 1);
Mc.c    = 1800 * ones(Mc.nz+1, 1);
Mc.K    = Mc.rho.*Mc.c.^2;
Mc.variables = {'vz','pz','h'};

Mc.pT.A = 1e4; % pressure perturbation amplitude
Mc.pT.T = 0.25; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
% Mc.G = @(t) 0; % external force from the conduit.

% fluid-filled fracture parameters.
Mf.w0 = 2;
Mf.Lx   = 2e3;
Mf.Ly   = 2e3;
Mf.nx = 50;
Mf.ny = 50;
Mf.nz = 20;
Mf.order = 2;
Mf.order_z = 2;
Mf.interp_order = 4;
Mf.xs = 0.75*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) 0;
Mf.xc = 0.75*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.15; % estimated ration of boundary layer.

Mf.variables = {'pxy','vx','vy'};

% fluid and solid properties.
Mf.rho = Mc.rho(1);
Mf.c    = Mc.c(1);
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 20;

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

x_obs = 0;
y_obs = 1000;
%% construct coupled models.
Model = coupledModel_NoMatrix(Mc, Mf);
%% time stepping:
CFL = 0.5;
skip = 40;
T = 100;
use_imex = false;
plot_simu = true;

% data_folder = '../../../../data/Videos/coupled_inviscid_L2000/';

% time stepping
hmin = min(Model.cond.geom.dz, Model.frac.geom.p.hx);
cmax = max(max(Model.cond.M.c), Model.frac.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt);

if ~ use_imex
    fun = @(u, t) Model.ft(u,t) + Model.Ai*u;
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
    fun = @(u,t) Model.ft(u,t);
end

%
disp('Time integration starts here:');
tic
Us = zeros(3, nt);
itr = 0;
for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt, fun, Model.Ai,L,U,p,q));
    end
%     
%    if i==100
%         t_elapse = toc;
%         fprintf( '%8.5f sec / step \n', t_elapse/100);
%         return
%    end
%     
    if mod(i, round(nt/100))==0
        toc;
        t_elapse = toc - tic;
        fprintf( '%% %f  finished, %8.5f sec / step', round( i*100/nt), t_elapse/round(nt/100));
        tic;
    end
    
    if t>Model.cond.M.pT.t
        
        Us(:,i) = Model.frac.disloc3d.eval_disp_p(Model.field(Model.u, [2 1]), x_obs, y_obs);
        
        if plot_simu
            if mod(i,skip) == 0
                itr = itr+1;
                %plot solution.
%                 [vz, pz, nz, ~, p_mat, vx_mat, ~] = Model.fields(Model.u);
                 [vz, pz, ~, p_mat, vx_mat, ~] = Model.fields(Model.u);
                uz = Model.cond.op.W2*Model.field(Model.u, [1, 1]);
                
                rm = Model.cond.geom.rm;
                z = Model.cond.geom.z;
                
                %fracture
                
                %plot solution.
                X = Model.frac.geom.p.X;
                Y = Model.frac.geom.p.Y;
                
                h1=figure(1);
                pcolor(X, Y, p_mat);
                xlabel('x'); ylabel('y');
                cmap;
                caxis([-1e3 1e3]);
                colorbar('location','eastoutside');
                shading INTERP;
                title(sprintf('t=%f',t));
                axis([0 Model.frac.M.Lx, 0 Model.frac.M.Ly]);
                daspect([1 1 1]);
%                 set(gcf, 'Position',[0 0 800 800]);
                set(gcf,'PaperPositionMode','auto');
                set(gca, 'fontsize',18)
%                 filename_px = [data_folder,'/px/px_',num2str(itr)];
%                 print(filename_px,'-dpng','-r0');
                
                h2 = figure(2);
                % plot the surface displacement vertical component.
                plot([1:i]*dt, Us(3,1:i),'-k');
                xlim([0 T]);xlabel('time (s)');ylabel('Uz');
                ylim([-1e-5, 1e-5]);
%                 filename_U = [data_folder,'/U/U_',num2str(itr)];
%                 set(gcf, 'Position',[0 0 800 300]);
                set(gcf,'PaperPositionMode','auto');
                set(gca, 'fontsize',18);
%                 print(filename_U,'-dpng');
                
                h3 = figure(3); %vx
                v_slice = squeeze(vx_mat(:, round(Model.frac.geom.vx.ny/2), :));
                [X_vx, Z_vx] = meshgrid(Model.frac.geom.vx.x, Model.frac.geom.vx.z);
                subplot(2,1,1);
                pcolor(X_vx, Z_vx, v_slice);
                xlabel('x'); ylabel('z');
                cmap;
                caxis([-2e-2 2e-2]);
                shading INTERP;
                title(sprintf('t=%f',t));
                set(gca,'fontsize',16);
                axis([0 Model.frac.M.Lx, 0 Model.frac.M.w0]);
                colorbar('location','southoutside');
                
                subplot(2,1,2);
                px_slice = squeeze(p_mat(round(Model.frac.geom.p.ny/2), :));
                plot(Model.frac.geom.p.x, px_slice);
                xlabel('x'); ylabel('p');
                xlim([0 Model.frac.M.Lx]);
                ylim([-1e3, 1e3]);
                set(gca,'fontsize',16);
                
%                 filename_vx = [data_folder,'/vx/vx_',num2str(itr)];
%                 set(gcf, 'Position',[0 0 800 600]);
                set(gcf,'PaperPositionMode','auto');
%                 print(filename_vx,'-dpng');

                %conduit.
                h4 = figure(4);
                subplot(1,2,1);
                plot(uz, z);
                xlim([-2e-1, 2e-1])
                ylim([0 Model.cond.M.L])
                xlabel('u'), ylabel('z (m)');
                
                subplot(1,2,2)
                plot(pz, z);
                xlim([-0.5e4,0.5e4]);
                ylim([0 Model.cond.M.L])
                xlabel('p'), ylabel('z (m)');
                suptitle(sprintf('t=%f',t));
%                 filename_pz_uz = [data_folder,'/pz_uz/pz_uz_',num2str(itr)];
%                 set(gcf, 'Position',[0 0 600 600]);
                set(gcf,'PaperPositionMode','auto');
%                 print(filename_pz_uz,'-dpng');
                drawnow
                
            end
            
        end
    end
end
toc;







