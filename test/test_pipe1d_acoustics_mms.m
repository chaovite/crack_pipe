% test pipe1d, no gravity, acoustics
% p=0 at z=1 and z=0
%
% mms solutions:
%           pz = sin(2*pi*z)*sin(2*pi*t)
%           vz = cos(2*pi*z)*cos(2*pi*t)
%

source_dir = '../source';
addpath(genpath(source_dir));

order = [2, 4, 6];
Nz     = [16, 32, 64, 128, 256, 512, 1024];
err     = zeros(length(order), length(Nz));
for i_order = 1: length(order)
    
    for i_nz  =  1: length(Nz)
        %%
        %-----------------------------------pipe parameter-----------------------------------
        % one section pipe.
        Mp.name  = 'pipe';
        Mp.L         = 1;        % length
        Mp.R         = 1;                                         % radius
        Mp.rho      = 1;  % density
        Mp.c         = 1;     % wave speed
        Mp.K         = Mp.rho.*Mp.c.^2;                       % fluid bulk modulus.
        Mp.mu      = 0;                                           %  viscosity.
        Mp.S         = pi*Mp.R.^2;                               % pipe surface area.
        Mp.nz       = Nz(i_nz);               % number grid points in each section.
        Mp.g         = 0;                            % gravitational acceleration
        Mp.order  = order(i_order);                             % order of spatial discretization.
        Mp.interfaces = [];
        Mp.bc.bt.type = 'p';
        Mp.bc.bt.f       = @(t) 0;
        Mp.bc.tp.type = 'h';
        Mp.bc.tp.f       = @(t) 0;
        model = pipe1d(Mp);
        %% ------------------------------------Create MMS -------------------------------------------
        z               = model.grd.z;
        v_mms_z = cos(2*pi*z);
        p_mms_z = sin(2*pi*z);
        h_mms_z = 1;
        u_mms_zt     = @(t) [v_mms_z*cos(2*pi*t); p_mms_z*sin(2*pi*t); h_mms_z*sin(2*pi*t)/(2*pi)];
        
        %% ------------------------------------Time integration--------------------------------------
        CFL = 0.5;
        T     = 1.5;
        plot_simu = true;
        
        % time stepping
        [cmax, hmin] = model.getCFL();
        dt = CFL*hmin/cmax;
        nt = ceil(T/dt);
        skip = ceil(nt/5);
        
        fun = @(u, t) model.fun_integrate(u, t);
        tic
        Us = zeros(3, nt);
        itr = 0;
        u = u_mms_zt(0);
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
                    vz          = u(model.indu.indv);
                    pz          = u(model.indu.indp);
                    u_mms  = u_mms_zt(t + dt);
                    vz_mms = u_mms(model.indu.indv);
                    pz_mms = u_mms(model.indu.indp);
                    z = model.grd.z;
                    figure(1)
                    subplot(2,1,1)
                    plot(z, pz,'k-',z, pz_mms,'r--');
                    legend('code','mms');
                    ylim([-1, 1]);
                    xlim([0 sum(Mp.L)])
                    title(sprintf('t = %8.2f', t));
                    subplot(2,1,2)
                    plot(z, vz,'k-',z, vz_mms,'r--');
                    legend('code','mms');
                    ylim([-1, 1]);
                    xlim([0 sum(Mp.L)])
                    title(sprintf('t = %8.2f', t));
                    suptitle(sprintf('mms test pipe1d inviscid acoustics, order=%d, nz=%d', Mp.order, Mp.nz));
                    drawnow();
                end
            end
        end
        
        t = t + dt; 
        % finish up so that the simulation ends exactly at t=T
        if t<T
            dt = T - t;
            t = t + dt;
            u = lsrk4(fun, u,t,dt);
        end
        
        % calculate the error
        E = model.energy_norm(); % energy norm
        u_mms = u_mms_zt(t);
        err(i_order, i_nz) = sqrt((u_mms -u)'*E*(u_mms - u));
    end
end
%%
figure;
for i = 1:length(order)
    loglog(Nz, err(i,:),'-d','linew',2);
    hold on;
end
hold off;
legend('2','4','6');
xlabel('nz','fontsize',16);
ylabel('err','fontsize',16);
set(gca,'fontsize',16);
save('pipe1d_inviscid_no_gravity_mms','order','Nz','err');
%%
dp = load('pipe1d_inviscid_no_gravity_mms','order','Nz','err');
dc = load('conduit_inviscid_no_gravity_mms','order','Nz','err');

figure;
for i = 1:length(dp.order)
    loglog(dp.Nz, dp.err(i,:),'-d','linew',2);
    hold on;
end

for i = 1:length(dc.order)
    loglog(dc.Nz, dc.err(i,:),'--d','linew',2);
    hold on;
end
hold off;
xlim([0, 1024])
legend('2-pipe1d','4-pipe1d','6-pipe1d','2-conduit2d','4-conduit2d','6-conduit2d');
title('convergence plot comparison');
xlabel('nz','fontsize',16);
ylabel('err','fontsize',16);
set(gca,'fontsize',16);









