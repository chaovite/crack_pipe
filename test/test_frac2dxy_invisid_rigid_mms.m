% test frac2dxy, inviscid fluid and rigid wall using manufactured solution.
clear;
source_dir = '../source';
addpath(genpath(source_dir));
%%
% Mf.order = 4;
% Mf.nx = 32;
% Mf.ny = 32;
% % fluid-filled fracture parameters.
% Mf.w0 = 1;
% Mf.Lx   = 1;
% Mf.Ly   = 1;
% 
% Mf.interp_order = 2;
% % Mf.xs = 0.5;
% % Mf.ys = 0.5;
% Mf.G  = @(t) 0; % no external forcing in the crack.
% % Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
% % Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.
% 
% % fluid and solid properties.
% Mf.rho = 1;
% Mf.c    = 1;
% Mf.K    = Mf.rho*Mf.c^2;
% Mf.mu = 0;
% Mf.isrigid = true;
% Mf.cp   = 5e3;
% Mf.cs    = 2.7e3;
% Mf.rhos = 3e3;
% Mf.Gs   = Mf.cs^2*Mf.rhos;
% Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
% %
% Model = frac2dxy(Mf);
%%
order = [2, 4, 6];
Nxy   = [8, 16, 32, 64, 128, 256];
% Nxy   = [8, 16, 32, 64];
err    = zeros(numel(order), numel(Nxy));

for i_order = 1:length(order)
    Mf.order = order(i_order);
    for i_n = 1: length(Nxy)
        fprintf('order %d, n  %d \n', order(i_order), Nxy(i_n));
        Mf.nx = Nxy(i_n);
        Mf.ny = Nxy(i_n);
        
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
        
        % fluid and solid properties.
        Mf.rho = 1;
        Mf.c    = 1;
        Mf.K    = Mf.rho*Mf.c^2;
        Mf.mu = 0;
        Mf.isrigid = true;
        Mf.cp   = 5e3;
        Mf.cs    = 2.7e3;
        Mf.rhos = 3e3;
        Mf.Gs   = Mf.cs^2*Mf.rhos;
        Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
        %%
        Model = frac2dxy(Mf);
        %% time stepping:
        CFL = 0.1;
        T = 1.5;
        use_imex = false;
        plot_simu = true;
        dim = Model.dimensions;
        %% creating mms solution.
        k = pi;
        x_p = Model.grd.p.vec(Model.grd.p.X);
        y_p = Model.grd.p.vec(Model.grd.p.Y);
        mu = Mf.mu;
        %source, no source.
        s_mms = @(t)  0;
        %p
        p_mms_xyz = block_matrix(dim,1,0);
        p_mms_xyz = block_matrix_insert(p_mms_xyz, dim,1,1,1, ...
            cos(k*x_p).*cos(k*y_p));
        p_mms_t      = @(t) sin(k*sqrt(2)*t);
        p_mms = @(t) p_mms_xyz * p_mms_t(t);
        
%         [x_vx, y_vx, ~] = Model.grd.vx.grd();
        x_vx = Model.grd.vx.vec(Model.grd.vx.X);
        y_vx = Model.grd.vx.vec(Model.grd.vx.Y);
        
%         [x_vy, y_vy, ~] = Model.grd.vy.grd();
        x_vy = Model.grd.vy.vec(Model.grd.vy.X);
        y_vy = Model.grd.vy.vec(Model.grd.vy.Y);
        
        % vx
        vx_mms_xyz = block_matrix(dim,1, 0);
        vx_mms_xyz = block_matrix_insert(vx_mms_xyz, dim,1,2,1,...
            - 1/sqrt(2)*sin(k*x_vx).*cos(k*y_vx));
        vx_mms_t     = @(t) cos(k*sqrt(2)*t);
        vx_mms          = @(t)  vx_mms_xyz*vx_mms_t(t);
        
        % vy
        vy_mms_xyz = block_matrix(dim,1, 0);
        vy_mms_xyz = block_matrix_insert(vy_mms_xyz, dim,1,3,1,...
            - 1/sqrt(2)*cos(k*x_vy).*sin(k*y_vy));
        vy_mms_t     = @(t) cos(k*sqrt(2)*t);
        vy_mms          = @(t)  vy_mms_xyz*vy_mms_t(t);
        
        U_mms  = @(t)  p_mms(t) + vx_mms(t) + vy_mms(t);
        %%
        % time stepping
%         hmin = min(Model.grd.p.hx, Model.grd.p.hy);
%         cmax = Mf.c;
        [cmax, hmin] = Model.getCFL();
        dt = CFL*hmin/cmax;
        nt = ceil(T/dt);
        skip = ceil(nt/5);
        
        if ~ use_imex
            A = Model.getA();
        else
            [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
            A = Model.Ae;
        end
        fun = @(u,t) A*u + Model.Fp*Model.M.G(t) + s_mms(t);
        tic;
        
        % initiation
        Model = Model.update(U_mms(0));
        
        for i=1:nt
            t = (i-1)*dt;
            
            if ~use_imex
                Model=Model.update(lsrk4(fun,Model.u,t,dt));
            else
                Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
            end
            
            if mod(i, round(nt/10))==0
                fprintf( '%% %f  finished', round(i*10/nt)*10);
                toc
            end
            
            if plot_simu
                if mod(i,skip) == 0
                    %plot solution.
                    X = Model.grd.p.X;
                    Y = Model.grd.p.Y;
                    p_vec = Model.field(Model.u, 1);
                    p_mat = reshape(p_vec, size(X));
                    
%                     X_ux = Model.grd.ux.X;
%                     Y_ux = Model.grd.ux.Y;
%                     vx_vec = Model.field(Model.u, 2);
                    
                    u_mms = U_mms(t + dt);
                    p_mms = reshape(Model.field(u_mms, 1), size(X));
                    
                    subplot(1,2,1);
                    pcolor(X, Y, p_mat);
                    xlabel('x'); ylabel('y');
                    cmap;
                    %             caxis([-20 20]);
                    shading INTERP;
                    title(sprintf('p t=%f',t));
                    axis([0 Model.M.Lx, 0 Model.M.Ly]);
                    colorbar;
                    
                    subplot(1,2,2);
                    pcolor(X, Y, p_mms);
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
loglog(Nxy, err(1,:), 'k-*', Nxy, err(2,:), 'r-o', Nxy, err(3,:), 'b-d');
legend({'2','4','6'});
xlabel('n');
ylabel('||err||');
set(gca,'fontsize',16);
save('test_frac2dxy_inviscid_rigid_convergence','Nxy', 'err');
%% compare with frac3d
d1=load('test_frac2dxy_inviscid_rigid_convergence','Nxy', 'err');
d2=load('inviscid_mms_convergence_test','Nxy', 'err');
figure;
loglog(d1.Nxy, d1.err(1,:), 'k-*', d1.Nxy, d1.err(2,:), 'r-o', d1.Nxy, d1.err(3,:), 'b-d');
hold on;
loglog(d2.Nxy, d2.err(1,:), 'k--*', d2.Nxy, d2.err(2,:), 'r--o', d2.Nxy, d2.err(3,:), 'b--d');
hold off;
legend({'2dxy-2','2dxy-4','2dxy-6','3d-2','3d-4','3d-6'});
xlabel('n');
ylabel('||err||');
set(gca,'fontsize',16);
% save('test_frac2dxy_inviscid_rigid_convergence','Nxy', 'err')





