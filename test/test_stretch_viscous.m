% test grid stretch, rigid wall limit with viscosity.
source_dir = '../source';
addpath(genpath(source_dir));
%% Simulation settings
skip                = 10;
order              = 2;
show_grid      = false;
nx                   = 50;
ny                   = 32;
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
wt =[2];
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
data_names = {'data_s','data_ns'};
stretch_paras = {[0.3 0.15],[0.3, 0.3]};

for s_i = 1: 2
    data_name = data_names{s_i};
    r_g = stretch_paras{s_i}(1);
    r_bl= stretch_paras{s_i}(2);
    disp(data_name);
for iL = 1: length(Lt)
    for iw= 1: length(wt)
        L = Lt(iL);
        T = max(L/1*400,200);
        w0= wt(iw);
        % Source term
        xs = 0;
        f0 = 1; t0 = 2/f0;
        interp_order = 2;
        % g = @(t) exp(-a^2*(t-t0).^2/2);% the unit of the source is the same as the velocity.
        g = @(t) ricker(t, f0, t0);% the unit of the source is the same as the velocity.
        
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
            
        end
        toc;
        
        %% construct data for transfer function calculation
        data.t = t_vec;
        data.u = g(t_vec);
        data.p = p0;
        data.Z = rho*c;
        switch data_name
        case 'data_s'
           data_s = compute_transfer_function(data, 5*data.t(end),[],100, true);
        case 'data_ns'
           data_ns = compute_transfer_function(data, 5*data.t(end),[],100, true);
        end
    end
end
end

%% compare to analytical solution
L =1;
w0 =2;
omega = 2*pi*data_s.f;
kesi = sqrt(-1i*(w0/2)^2*omega/(mu/rho));
K_dispersion = (1-tanh(kesi)./kesi).^(-1/2);
k_dispersion = K_dispersion.*omega/c;
k_dispersion(1) = 0;
Fa = -1i*1./K_dispersion.*tan(k_dispersion*L);

%low frequency limit
f_slope = [0: 0.01: 0.5];
F_slope = -1i*2*pi*f_slope/c*L;

plot(data_s.f, abs(data_s.F),'b-');hold on;
plot(data_ns.f, abs(data_ns.F),'k-');hold on;
plot(data_s.f, abs(Fa),'r');hold on;
plot(f_slope, abs(F_slope),'r--');hold off;
legend({'stretch','no stretch','analytical','low freq limit'})

xlim([0 6]);shg;
set(gca,'fontsize',18);
title(sprintf('L= %d km, w= %d m', L, w0))
