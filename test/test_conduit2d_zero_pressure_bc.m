% test conduit2D zero pressure bc with gravity
% compare to analytical solution.

source_dir = '../source';
addpath(genpath(source_dir));
%-----------------------------------pipe parameter-----------------------------------
% one section pipe.
Mp.name  = 'pipe';
Mp.nz       = 50;                % number grid points in th pipe.
Mp.nr        = 50;    % number of grid points in R direction.
Mp.L         = 300;                           % length
Mp.R         = 5;                                     % radius
Mp.rho      = 1000*ones(Mp.nz+1,1);    % density
Mp.c         = 1000*ones(Mp.nz+1,1);     % wave speed
Mp.K         =  Mp.c.^2.*Mp.rho;     % fluid bulk modulus.
Mp.mu      = 50;                                     %  viscosity.
Mp.S         = pi*Mp.R.^2;                     % pipe surface area.
Mp.g        = 10;                                     % gravitational acceleration
Mp.order  = 4;                % order of spatial discretization.
Mp.order_r = 4;
Mp.interface_split=false;
Mp.with_exsolution = false;
Mp.interfaces = [];
Mp.bc.bt.type = 'p';
Mp.bc.bt.f       = @(t) 0;
Mp.bc.tp.type = 'h';
Mp.bc.tp.f       = @(t) 1e5*exp(-(t-3).^2);
model = conduit(Mp);
%% ------------------------------------Time integration--------------------------------------
CFL = 0.5;
T     = 100;
plot_simu = true;

% time stepping
cmax = max(Mp.c);
hmin  = Mp.L/Mp.nz;
dt = CFL*hmin/cmax;
nt = ceil(T/dt);
skip = 10;

fun = @(u, t) model.A*u + model.Fp*Mp.bc.tp.f(t);
tic
Us = zeros(3, nt);
itr = 0;
u = model.u;

%% visualize in time domain.

for i=1:nt
    t = (i-1)*dt;
    % invicid, integrate explicitly
    u = lsrk4(fun,u,t,dt);
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i, skip) == 0
            indv = field_indices(model.dimensions(),1);
            indp = field_indices(model.dimensions(),2);
            indh = field_indices(model.dimensions(),3);
            vz          = u(indv);
            uz          = model.op.W2*vz;
            pz          = u(indp);
            h            = u(indh);
            z = model.geom.z;
            
            figure(1)
            subplot(2,1,1)
            plot(z, pz,'k-');
            ylim([-1, 1]*1e5);
            xlim([0 sum(Mp.L)]);
            ylabel('pressure');
            title(sprintf('t = %8.2f', t));
            
            subplot(2,1,2)
            plot(z, uz,'k-');
            ylim([-1, 1]);
            xlim([0 sum(Mp.L)]);
            ylabel('velocity');
            title(sprintf('t = %8.2f', t));
            
            drawnow();
        end
    end
end


