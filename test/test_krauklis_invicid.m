% test inviscid krauklis wave with no zero initial condition.
source_dir = '../source';
addpath(genpath(source_dir));
%% Simulation settings
clear; 

skip            = 2;
order           = 4;
show_grid       = false;
nx              = 256;
ny              = 20;
use_imex       = true;
test_stability  = false;
plot_source   = false;
plot_simu      = true;
isrigid = false;

% grid stretching parameters in y direction.
r_g  = 0.3;% ratio of grid points in boundary layer.
r_bl = 0.3;% ratio of boundary layer thickness to fracture width.


% Fluid properties


Lt = [8];
wt =[0.5];
rho = 2.5;
c    = 1.8; 
K   = rho*c^2;
mu  = 0; %kPa.s
cp    = 5;
cs    = 2.7;
rhos = 3;
G   = cs^2*rhos*1e-3;
nu  = ((cp/cs)^2-2)/((cp/cs)^2-1)/2;
Gs = G/(1-nu);

for iL = 1: length(Lt)
    for iw= 1: length(wt)
        L = Lt(iL);
        T = 10;
        w0= wt(iw);
        
        % Source term
        xs = 0;
        f0 = 1; t0 = 2/f0;
        interp_order = 2;
        % g = @(t) exp(-a^2*(t-t0).^2/2);% the unit of the source is the same as the velocity.
        %         g = @(t) ricker(t, f0, t0);% the unit of the source is the same as the velocity.
        g = @(t) 0;
        
        % specify initial condition;
        x_c = L/2; x_dist = L/200;
        fp = @(x) exp(-(x-x_c).^2/(2*x_dist^2));
        
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
            op.(gt)         = operators(gt,order,geom.(gt),'strong');
        end
        metrics = struct();
        
        f = fluid(geom,metrics,op,w0,rho,K,mu, G, nu, isrigid);
        
        %initial condition:
        p_t0 = fp(geom.p.x);
        dim  = f.dimensions();
        v_t0 = zeros(dim(2), 1);
        u_t0 = [p_t0; v_t0];
        f = f.init(u_t0);
        
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
                if mod(i,skip) == 0
                    sol.p = f.field(f.u,1);
                    sol.v = f.field(f.u,2);
                    sol_v = geom.mm.grd(sol.v);
                    subplot(4,1,1:2);
                    pcolor(geom.mm.X,geom.mm.Y,geom.mm.grd(sol.v));
                    xlabel('x'); ylabel('y');
                    cmap;
                    caxis([-1 1]);
                    shading INTERP;
                    colorbar('location','southoutside');
                    title(sprintf('t=%f',t));
                    xp = geom.p.x;
                    xv = geom.m.x;
                    
                    N = nx;
                    xa =[-N/2+1:N/2]*L/N;
                    b0 = fp(xa+L/2);
                    [p_a, v_a] = krauklis_wave(xa,t,b0,c,rho,Gs,w0);
                    
                    subplot(4,1,3);
                    plot(geom.m.x, sol_v(round(ny/2), :),'k-', xa+L/2, real(v_a),'r-');
                    legend('FD','FFT');
                    xlabel('x');
                    ylabel('v');
                    ylim([-1 1]);
                    
                    subplot(4,1,4);
                    plot(geom.p.x, sol.p, 'k-',xa+L/2,real(p_a), 'r-');
                    legend('FD','FFT');
                    xlabel('x');
                    ylabel('p');
                    ylim([-1 1]);
                    suptitle('Test invicid Krauklis wave');
                    drawnow;
%                     pause();
                    pause_time = [0.5: 0.5: 5];
                    if any(t > (pause_time-dt/2) & t <= (pause_time+dt/2))
                        pause(2)
                    end
                end
              end
            
            if ~use_imex
                f.u = lsrk4(fun,f.u,t,dt);
            else
                f.u = imex_ark4_lu(f.u,t,dt,fun,f.Ai,L_lu,U_lu,p_lu,q_lu);
            end
            
        end
        toc;
    end
end
