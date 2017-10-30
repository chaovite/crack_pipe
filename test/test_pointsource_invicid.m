% test rigid wall limit with no viscosity and with a point source.
source_dir = '../source';
addpath(genpath(source_dir));
%% Simulation settings
clear; 

skip            = 10;
order           = 6;
show_grid       = false;
nx              = 100;
ny              = 30;
use_imex       = true;
test_stability  = false;
plot_source   = false;
plot_simu      = true;
isrigid = true;

% grid stretching parameters in y direction.
r_g  = 0.3;% ratio of grid points in boundary layer.
r_bl = 0.15;% ratio of boundary layer thickness to fracture width.


% Fluid properties
% TODO: w0 is constant for now.
rho = 1;
c    = 1;
K   = rho*c^2;
mu  = 0;%kPa.s
cp    = 5;
cs    = 2.7;
rhos = 3;
G   = cs^2*rhos*1e-3;
nu  = ((cp/cs)^2-2)/((cp/cs)^2-1)/2;

% Lt = [0.5, 0.8, 1, 2];
% wt =[1, 2, 4]*1e-3;

Lt = [3];
wt =[1];

for iL = 1: length(Lt)
    for iw= 1: length(wt)
        L = Lt(iL);
        T = max(L/1*400,200);
        w0= wt(iw);
        
        % Source term
        xs = 0.5*L;
        f0 = 2; t0 = 2/f0;
        interp_order = 2;
%         g = @(t) exp(-f0^2*(t-t0).^2/2);% the unit of the source is the same as the velocity.
        g = @(t) ricker(t, f0, t0);% the unit of the source is the same as the velocity.
%         g = @(t) 0;
        
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
        
        f = fluid(geom,metrics,op,w0,K,rho,mu, G, nu, isrigid);
        
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
             if plot_simu
                if mod(skip,1) == 0
                    sol.p = f.field(f.u,1);
                    sol.v = f.field(f.u,2);
                    sol_v = geom.mm.grd(sol.v);

                    subplot(4,1,1:2);
                    pcolor(geom.mm.X,geom.mm.Y,geom.mm.grd(sol.v));
                    xlabel('x'); ylabel('y');
                    cmap;
                    caxis([-0.5 0.5]*K/c);
                    shading INTERP;
                    colorbar('location','southoutside');
                    title(sprintf('t=%f',t));
                    
                    subplot(4,1,3);
                    xv = geom.m.x;
                    v_analytic = zeros(size(xv));
                    ind_vp = (xv-xs)/c<t & (xv-xs)/c>=0;
                    ind_vm = (xv-xs)/c>-t & (xv-xs)/c<=0;
                    v_analytic(ind_vp) = g(t-(xv(ind_vp)-xs))/2;
                    v_analytic(ind_vm) = -g(t+(xv(ind_vm)-xs))/2;
                    plot(xv, sol_v(round(ny/2), :)','k-',xv, v_analytic, 'r*');
                    xlabel('x');
                    ylabel('v');
                    ylim([-0.5 0.5]*K/c/rho/c);
                    
                    subplot(4,1,4);
                    xp = geom.p.x;
                    p_analytic = zeros(size(xp));
                    ind_pp = (xp-xs)/c<t & (xp-xs)/c>=0;
                    ind_pm = (xp-xs)/c>-t & (xp-xs)/c<=0;
                    p_analytic(ind_pp) = K/(2*c)*g(t-(xp(ind_pp)-xs)/c);
                    p_analytic(ind_pm) = K/(2*c)*g(t+(xp(ind_pm)-xs)/c);
                    plot(xp, sol.p, 'k-',xp, p_analytic, 'r*');  
                    xlabel('x');
                    ylabel('p');
                    ylim([-0.5 0.5]*K/c);
                    suptitle('Test point source, inviscid fluid, rigid wall')
                    drawnow;
%                     pause();
                    if t >=t0-dt/2 && t<=t0+dt/2
                        pause();
                    end
                    if t >=t0+L/3/c-dt/2 && t<=t0+L/3/c+dt/2
                        pause();
                    end
                end
            end
            if ~use_imex
                f.u = lsrk4(fun,f.u,t,dt);
            else
                f.u = imex_ark4_lu(f.u,t,dt,fun,f.Ai,L_lu,U_lu,p_lu,q_lu);
            end
            if mod(i, round(nt/100))==0
                fprintf( '%% %f  finished', round(i*100/nt));
                toc;
            end

        end
        toc;
    end
end
