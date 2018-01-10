% test frac3d with viscosity and rigid wall using manufactured solution.
% the mms solution is checked in script mms_frac3d_rigid_viscous.m
source_dir = '../source';
addpath(genpath(source_dir));
%% convergence test.
order = [2, 4, 6];
Nxyz   = [8, 16, 32];
err    = zeros(numel(order), numel(Nxyz));

for i_order = 1:length(order)
    Mf.order = order(i_order);
    for i_n = 1: length(Nxyz)
        fprintf('order %d, n  %d \n', order(i_order), Nxyz(i_n));
        Mf.nx = Nxyz(i_n);
        Mf.ny = Nxyz(i_n);
        Mf.nz = Nxyz(i_n);
        
        if Mf.order==6 && Mf.nx==8
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
        Mf.r_bl =  0.3; % estimated ration of boundary layer.
        
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
        Model = frac3d_o(Mf);
        %% time stepping:
        CFL = 0.5;
        T = 1;
        use_imex = true;
        plot_simu = true;
        
        dim = Model.dimensions;
        %% creating mms solution.
        k = pi;
        x_p = Model.geom.p.vec(Model.geom.p.X);
        y_p = Model.geom.p.vec(Model.geom.p.Y);
        
        mu = Mf.mu;
        
        % source.
        s_mms_t = @(t) 0;
        s_mms_xyz = block_matrix(dim,1, 0);
        s_mms = @(t)  s_mms_xyz * s_mms_t(t);
        
        %p
        
        p_mms_xyz = block_matrix(dim,1,0);
        p_mms_t      = @(t) 0;
        p_mms = @(t) p_mms_xyz * p_mms_t(t);
        
        [x_vx, y_vx, z_vx] = Model.geom.vx.mesh();
        x_vx = Model.geom.vx.vec(x_vx);
        y_vx = Model.geom.vx.vec(y_vx);
        z_vx = Model.geom.vx.vec(z_vx);
        
        [x_vy, y_vy, z_vy] = Model.geom.vy.mesh();
        x_vy = Model.geom.vy.vec(x_vy);
        y_vy = Model.geom.vy.vec(y_vy);
        z_vy = Model.geom.vy.vec(z_vy);
        
        % vx
        vx_mms_xyz = block_matrix(dim,1, 0);
        vx_mms_xyz = block_matrix_insert(vx_mms_xyz, dim,1,2,1,...
            sin(k*x_vx).*cos(k*y_vx).*sin(k*z_vx));
        
        vx_mms_t     = @(t) exp(-mu*k^2*t);
        vx_mms          = @(t) vx_mms_xyz*vx_mms_t(t);
        
        % vy
        vy_mms_xyz = block_matrix(dim,1, 0);
        vy_mms_xyz = block_matrix_insert(vy_mms_xyz, dim,1,3,1,...
            -cos(k*x_vy).*sin(k*y_vy).*sin(k*z_vy));
        
        vy_mms_t     = @(t) exp(-mu*k^2*t);
        vy_mms          = @(t) vy_mms_xyz*vy_mms_t(t);
        
        % ux
        x_ux = Model.geom.ux.vec(Model.geom.ux.X);
        y_ux = Model.geom.ux.vec(Model.geom.ux.Y);
        ux_mms_xy    = (1-cos(k))/k*sin(k*x_ux).*cos(k*y_ux);
        ux_mms_t       =@(t) exp(-mu*k^2*t);
        ux_mms          = @(t)  ux_mms_xy*ux_mms_t(t);
        
        % uy
        x_uy = Model.geom.ux.vec(Model.geom.uy.X);
        y_uy = Model.geom.uy.vec(Model.geom.uy.Y);
        uy_mms_xy    = (1-cos(k))/k*cos(k*x_uy).*sin(k*y_uy);
        uy_mms_t       =@(t) exp(-mu*k^2*t);
        uy_mms          = @(t)  uy_mms_xy*uy_mms_t(t);
        
        U_mms  = @(t)  p_mms(t) + vx_mms(t) + vy_mms(t);
        %%
        % time stepping
        hmin = min(Model.geom.p.hx, Model.geom.p.hy);
        cmax = Mf.c;
        dt = CFL*hmin/cmax;
        nt = ceil(T/dt);
        skip = ceil(nt/10);
        
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
                    X_ux = Model.geom.ux.X;
                    Y_ux = Model.geom.ux.Y;
                    ux_vec =  Model.op.vx.Wz3*Model.field(Model.u, 2);
                    
                    ux_mat = reshape(ux_vec, size(X_ux));
                    ux_mms_mat = reshape(ux_mms(t+dt), size(X_ux));
                    
                    subplot(1,2,1);
                    pcolor(X_ux, Y_ux, ux_mat);
                    xlabel('x'); ylabel('y');
                    cmap;
                    %             caxis([-20 20]);
                    shading INTERP;
                    title(sprintf('p t=%f',t));
                    axis([0 Model.M.Lx, 0 Model.M.Ly]);
                    colorbar;
                    subplot(1,2,2);
                    pcolor(X_ux, Y_ux, ux_mms_mat);
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
loglog(Nxyz, err(1,:), 'k-*', Nxyz, err(2,:), 'r-o', Nxyz, err(3,:), 'b-d');
legend({'2','4','6'});
xlabel('n');
ylabel('||err||');
set(gca,'fontsize',16);


