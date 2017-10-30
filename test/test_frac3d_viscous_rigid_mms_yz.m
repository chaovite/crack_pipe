% test frac3d with viscosity and rigid wall using manufactured solution.
% the mms solution is checked in script mms_frac3d_rigid_viscous.m

% this script use solution that is constant in y direction and varies only
% in y z direction.
source_dir = '../source';
addpath(genpath(source_dir));

%% convergence test.
order = [2, 4, 6];
Nx     = [6, 8, 12];
Nyz   = [8, 16, 32, 64];
err    = zeros(numel(order), numel(Nyz));

for i_order = 1:length(order)
    Mf.order = order(i_order);
    Mf.nx = Nx(i_order);
    for i_n = 1: length(Nyz)
        fprintf('order %d, n  %d \n', order(i_order), Nyz(i_n));
        Mf.ny = Nyz(i_n);
        Mf.nz = Nyz(i_n);
        
        if Mf.order==6 && Mf.ny==8
            err(i_order, i_n) = nan;
            continue
        end
        
        % fluid-filled fracture parameters.
        Mf.w0 = 1;
        Mf.Lx   = 1;
        Mf.Ly   = 1;
        Mf.interp_order = 2;
        % Mf.xs = 0.5;
        % Mf.ys = 0.5;
        Mf.G  = @(t) 0; % no external forcing in the crack.
        % Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
        % Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.
        
        % Mf.lc  = 2*pi*Mc.R;%coupling length scale.
        % Mf.S = Mf.lc * Mf.w0;
        Mf.isrigid = true;
        Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
        Mf.r_bl =  0.15; % estimated ration of boundary layer.
        
        % fluid and solid properties.
        Mf.rho = 1;
        Mf.c    = 1;
        Mf.K    = Mf.rho*Mf.c^2;
        Mf.mu = 1e-2;
        Mf.cp   = 5e3;
        Mf.cs    = 2.7e3;
        Mf.rhos = 3e3;
        Mf.Gs   = Mf.cs^2*Mf.rhos;
        Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
        %%
        Model = frac3d(Mf);
        %% time stepping:
        CFL = 0.5;
        mu = Mf.mu;
        T = 1;
        use_imex = true;
        plot_simu = true;
        dim = Model.dimensions;
        %% creating mms solution.
        k = pi;
        y_p = Model.geom.p.vec(Model.geom.p.Y);
        % source.
        s_mms_t = @(t) exp(-mu*k^2*t);
        s_mms_xyz = block_matrix(dim,1, 0);
        s_mms_xyz = block_matrix_insert(s_mms_xyz, dim,1,1,1, ...
            (1-cos(k)).*cos(k*y_p));
        s_mms = @(t)  s_mms_xyz * s_mms_t(t);
        
        %p
        
        p_mms_xyz = block_matrix(dim,1,0);
        p_mms_t      = @(t) 0;
        p_mms = @(t) p_mms_xyz * p_mms_t(t);
        
        [~, y_vx, z_vx] = Model.geom.vx.mesh();
        y_vx = Model.geom.vx.vec(y_vx);
        z_vx = Model.geom.vx.vec(z_vx);
        
        [~, y_vy, z_vy] = Model.geom.vy.mesh();
        y_vy = Model.geom.vy.vec(y_vy);
        z_vy = Model.geom.vy.vec(z_vy);
        
        % vy
        vy_mms_xyz = block_matrix(dim,1, 0);
        vy_mms_xyz = block_matrix_insert(vy_mms_xyz, dim,1,3,1,...
            sin(k*y_vy).*sin(k*z_vy));
        vy_mms_t     = @(t) exp(-mu*k^2*t);
        vy_mms          = @(t) vy_mms_xyz*vy_mms_t(t);
        
        % vx
        vx_mms_xyz1 = block_matrix(dim,1, 0);
        vx_mms_t1     = @(t) 0;
        vx_mms          = @(t) vx_mms_xyz1*vx_mms_t1(t);
        
        % uy
        y_uy = Model.geom.uy.vec(Model.geom.uy.Y);
        uy_mms_xz    = (1-cos(k))/k*sin(k*y_uy);
        uy_mms_t       = exp(-mu*k^2*t);
        uy_mms          = @(t)  uy_mms_xz*uy_mms_t;
        
        U_mms  = @(t)  p_mms(t) + vx_mms(t) + vy_mms(t);
        %%
        % time stepping
        hmin = min(Model.geom.p.hx, Model.geom.p.hy);
        cmax = Mf.c;
        dt = CFL*hmin/cmax;
        nt = ceil(T/dt);
        skip = ceil(nt/5);
        
        if ~ use_imex
            A = Model.Ae + Model.Ai;
        else
            [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
            A = Model.Ae;
        end
        fun = @(u,t) A*u + Model.Fp*Model.M.G(t) + s_mms(t);
        tic;
        
        % initiation
        Model = Model.init(U_mms(0));
        
        for i=1:nt
            t = (i-1)*dt;
            
            if ~use_imex
                Model=Model.update(lsrk4(fun,Model.u,t,dt));
            else
                Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
            end
            
            if mod(i, round(nt/10))==0
                fprintf( '%% %f  finished', round(i*10/nt)*10);
                toc;
            end
            
            if plot_simu
                if mod(i,skip) == 0
                    %plot solution ux.                    
                    X_uy = Model.geom.uy.X;
                    Y_uy = Model.geom.uy.Y;
                    uy_vec =  Model.op.vy.Wz3*Model.field(Model.u, 3);
                    
                    uy_mat = reshape(uy_vec, size(X_uy));
                    uy_mms_mat = reshape(uy_mms(t+dt), size(X_uy));
                    
                    subplot(1,2,1);
                    pcolor(X_uy, Y_uy, uy_mat);
                    xlabel('x'); ylabel('y');
                    cmap;
                    %             caxis([-20 20]);
                    shading INTERP;
                    title(sprintf('p t=%f',t));
                    axis([0 Model.M.Lx, 0 Model.M.Ly]);
                    colorbar;
                    subplot(1,2,2);
                    pcolor(X_uy, Y_uy, uy_mms_mat);
                    xlabel('x'); ylabel('y');
                    cmap;
                    %             caxis([-20 20]);
                    shading INTERP;
                    title(sprintf('p mms t=%f',t));
                    axis([0 Model.M.Lx, 0 Model.M.Ly]);
                    colorbar;
                    drawnow;
                end
                
            end
        end
        toc;
        
        t = t + dt; 
        % finish up so that the simulation ends exactly at t=T
        if t<T
            dt = T - t;
            t = t + dt;
            if ~ use_imex
                A = Model.Ae + Model.Ai;
            else
                [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
                A = Model.Ae;
            end
            fun = @(u,t) A*u + Model.Fp*Model.M.G(t) + s_mms(t);
            if ~use_imex
                Model=Model.update(lsrk4(fun,Model.u,t,dt));
            else
                Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
            end
        end
        u_mms = U_mms(t);
        err(i_order, i_n) = sqrt((u_mms - Model.u)'*Model.E*(u_mms - Model.u));
    end
end

%% plot convergence:
figure;
loglog(Nyz, err(1,:), 'k-*', Nyz, err(2,:), 'r-o', Nyz, err(3,:), 'b-d');
legend({'2','4','6'});
xlabel('n');
ylabel('||err||');
set(gca,'fontsize',16);


