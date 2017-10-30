% test rigid wall limit with viscosity.
source_dir = '../source';
addpath(genpath(source_dir));
%% Simulation settings
clear; close all;

skip            = 10;
order           = 6;
show_grid       = false;
nx              = 100;
ny              = 50;
use_imex       = true;
test_stability  = false;
plot_source   = false;

% grid stretching parameters in y direction.
r_g  = 0.3;% ratio of grid points in boundary layer.
r_bl = 0.3;% ratio of boundary layer thickness to fracture width.


% Fluid properties
% TODO: w0 is constant for now.

% Lt = [0.5, 0.8, 1, 2];
% wt =[1, 2, 4]*1e-3;

Lt = [1];
wt =[1];
rho = 1;
c    = 1;
K   = rho*c^2;
mu  = 1e-1;%kPa.s
cp    = 5;
cs    = 2.7;
rhos = 3;
G   = cs^2*rhos*1e-3;
nu  = ((cp/cs)^2-2)/((cp/cs)^2-1)/2;
isrigid = true;

%%

for iL = 1: length(Lt)
    for iw= 1: length(wt)
        L = Lt(iL);
        T = max(L/1*400,200);
        w0= wt(iw);
        %         L   = 1e3;
        %         w0  = 2;
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
        
        %% Initialize rectangular grid and construct staggered grid operators
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
        f = fluid(geom,metrics,op,w0,rho, K,mu, G, nu, isrigid);
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
        %% time stepping
        fprintf('L=%f km, w0=%f m \n', L, w0);
        c =  sqrt(K/rho);
        dt = 0.5*op.p.hx/c;
        nt = ceil(T/dt);
        t_vec = [0:dt:nt*dt];
        p0 = zeros(1, nt+1);
        
        if ~ use_imex
            f.A = f.Ae + f.Ai;
        else
            [L_lu,U_lu,p_lu,q_lu,~] = imex_ark4_get_lu(f.Ai,dt);
            f.A = f.Ae;
        end
        fun = @(u,t) f.A*u + e*g(t);
        
        tic
        for i=1:nt
            t = (i-1)*dt;
            if ~use_imex
                f.u = lsrk4(fun,f.u,t,dt);
            else
                f.u = imex_ark4_lu(f.u,t,dt,fun,f.Ai,L_lu,U_lu,p_lu,q_lu);
            end
            if mod(i, round(nt/100))==0
                fprintf( '%% %f  finished', round(i*100/nt));
                toc;
            end
            sol.p = f.field(f.u,1);
            p0(i+1)     = sol.p(1);
            
            %             if mod(skip,1) == 0
            %                 sol.p = f.field(f.u,1);
            %                 sol.v = f.field(f.u,2);
            %                 sol_v = geom.mm.grd(sol.v);
            %                 Pressure
            %                 plot(geom.p.x,sol.p,'o-');
            %                 ylim([-1 1]);
            %                 Velocity
            %                 subplot(4,1,1:2);
            %                 pcolor(geom.mm.X,geom.mm.Y,geom.mm.grd(sol.v));
            %                 xlabel('x'); ylabel('y');
            %                 cmap;
            %                 caxis([-1 1]);
            %                 shading INTERP;
            %                 colorbar('location','southoutside');
            %                 title(sprintf('t=%f',t));
            %
            %                 subplot(4,1,3);
            %                 plot(geom.m.x, sol_v(round(ny/2), :)');
            %                 xlabel('x');
            %                 ylabel('v');
            %                 ylim([-1 1]);
            %                 subplot(4,1,4);
            %                 plot(geom.p.x, sol.p);
            %                 xlabel('x');
            %                 ylabel('p');
            %                 ylim([-1 1]);
            %
            %                             %     sol_v = geom.mm.grd(sol.v);
            %                             %
            %                             %     plot(geom.mm.x,sol_v(10,:), geom.mm.x, zeros(size(geom.mm.x)),'o');
            %                             %     ylim([-1, 1]);
            %                             %     plot(geom.p.x, f.Ks_inv*sol.p);
            %                 drawnow;
            %                 if t >=2-dt/2 && t<=2+dt/2
            %                     pause();
            %                 end
            %                 if t >=3-dt/2 && t<=3+dt/2
            %                     pause();
            %                 end
            %             end
        end
        toc;
        
        %% construct data for transfer function calculation
        data_root= '/Users/Chaovite/Documents/data/basalt_rock/';
        test_name = sprintf('test_rigid_nx%dny%d_no_stretch',nx, ny);
        name = ['L',num2str(L*1000), 'w', num2str(w0),test_name,'.mat'];
        data.t = t_vec;
        data.u = g(t_vec);
        data.p = p0;
        data.Z = rho*c;
        save([data_root,name], 'data');
    end
end

%% now construct fracture transfer function.
%
for iL = 1: length(Lt)
    for iw= 1: length(wt)
        L = Lt(iL);
        w0= wt(iw);
        name = ['L',num2str(L*1000), 'w', num2str(w0),test_name,'.mat'];
        data_load = load([data_root, name]);
        data_in = data_load.data;
        data_in.p = data_in.p;
        T = compute_transfer_function(data_in, 5*data_in.t(end),[],1, true);
        save([data_root, name], 'T');
    end
end

%% compare to analytical solution
L =1;
w0=1;
data_root= '/Users/Chaovite/Documents/data/basalt_rock/';
test_name1 = sprintf('test_rigid_nx%dny%d',50, 30);
test_name2 = sprintf('test_rigid_nx%dny%d',50,100);
test_name3 = sprintf('test_rigid_nx%dny%d',200,30);
test_name4 = sprintf('test_rigid_nx%dny%d',100,200);
test_name5 = sprintf('test_rigid_nx%dny%d_no_stretch',100,200);
test_name6 = sprintf('test_rigid_nx%dny%d_no_stretch',100,50);
data_test1= load([data_root,'L',num2str(L*1000), 'w', num2str(w0),test_name1,'.mat']);
data_test2= load([data_root,'L',num2str(L*1000), 'w', num2str(w0),test_name2,'.mat']);
data_test3= load([data_root,'L',num2str(L*1000), 'w', num2str(w0),test_name3,'.mat']);
data_test4= load([data_root,'L',num2str(L*1000), 'w', num2str(w0),test_name4,'.mat']);
data_test5= load([data_root,'L',num2str(L*1000), 'w', num2str(w0),test_name5,'.mat']);
data_test6= load([data_root,'L',num2str(L*1000), 'w', num2str(w0),test_name6,'.mat']);
omega = 2*pi*data_test1.T.f;
kesi = sqrt(-1i*(w0/2)^2*omega/(mu/rho));
K_dispersion = (1-tanh(kesi)./kesi).^(-1/2);
k_dispersion = K_dispersion.*omega/c;
k_dispersion(1) = 0;
Fa = -1i*1./K_dispersion.*tan(k_dispersion*L);

%low frequency limit
f_slope = [0: 0.01: 0.5];
F_slope = -1i*2*pi*f_slope/c*L;

plot(data_test1.T.f, abs(data_test1.T.F),'b');hold on;
plot(data_test2.T.f, abs(data_test2.T.F),'g');hold on;
plot(data_test3.T.f, abs(data_test3.T.F),'m-');hold on;
plot(data_test4.T.f, abs(data_test4.T.F),'k-');hold on;
plot(data_test5.T.f, abs(data_test5.T.F),'y*');hold on;
plot(data_test6.T.f, abs(data_test6.T.F),'c*');hold on;
plot(data_test1.T.f, abs(Fa),'r');hold on;
plot(f_slope, abs(F_slope),'r--');hold off;

% legend({'numerical 50, 30','numerical 50, 100','numerical 200, 30','analytical'})
legend({'numerical 50, 30','numerical 50, 100','numerical 200, 30',...
    'numerical 100, 200','numerical 100, 200 no stretch','numerical 100, 50 no stretch','analytical','low freq limit'})
% legend({'numerical 50, 30','numerical 50, 100','numerical 200, 30','numerical 100, 200','analytical','low freq limit'})
xlim([0 6]);shg;
set(gca,'fontsize',18);
title(sprintf('L= %d km, w= %d m', L, w0))
