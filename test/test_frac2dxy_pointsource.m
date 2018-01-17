% test frac2dxy, inviscid 2D acoustics from point source.
clear;
source_dir = '../source';
addpath(genpath(source_dir));
%% material properties.
%%
% fluid-filled fracture parameters.
Mf.w0 = 2;
Mf.Lx   = 4;
Mf.Ly   = 4;
Mf.nx = 200;
Mf.ny = 200;
Mf.order = 4;
Mf.interp_order = 6;
Mf.xs = 0.5*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
% Mf.G  = @(t) 5e4/(Mc.c(1)*Mc.rho(1))*ricker(t, 1, 2); % no external forcing in the crack.
% Mf.G  = @(t) 1e5/(Mf.c*Mf.rho)*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2)
Mf.G  = @(t) ricker(t, 5, 2/5);
Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = true;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho = 2;
Mf.c    = 1;
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 0;
%% construct coupled models.
Model = frac2dxy(Mf);
%% time stepping:
CFL = 0.25;
skip = 10;
T = Mf.Lx/2/Mf.c;
use_imex = false;
plot_simu = false;

% time stepping
[cmax, hmin] = Model.getCFL();
dt = CFL*hmin/cmax;
nt = ceil(T/dt);
t_n = (0:nt)*dt;

%% query points to compare the solutions;
x = Model.grd.p.x;
y = Model.grd.p.y;

[~, indx] =  find(x>Mf.xs & x<0.75*Mf.Lx);
[~, indy] =  find(y>Mf.ys & y<0.75*Mf.Ly);

indr = indy + (indx - 1)*Model.grd.p.ny;

xq = x(indx);
yq = y(indy);

rq = sqrt((xq-Mf.xs).^2 + (yq-Mf.ys).^2);
g = Mf.G(0:dt:nt*dt)*Mf.rho/Mf.w0;
p_a= zeros(length(rq), nt);
Nrq = length(rq);
tic
for i = 1 : Nrq
    [p_a(i,:), t_a] = acoustics2D_pointsource(rq(i), Mf.c, g, dt);
end
assert(max(abs(imag(p_a(:)./real(p_a(:)))))<1e-6);
p_a = real(p_a);
toc

%%
p_n= zeros(length(rq), nt+1);

if ~ use_imex
    A = Model.getA();
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end
fun = @(u,t) A*u + Model.Fp*Model.M.G(t);

tic
for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
    end
    
    p_vec = Model.field(Model.u, 1);
    
    p_n(:, i +1) = p_vec(indr);
    
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if plot_simu
        if mod(i,skip) == 0
            %plot solution.
            X = Model.grd.p.X;
            Y = Model.grd.p.Y;
            p_mat = reshape(p_vec, size(X));
            pcolor(X, Y, p_mat);
            xlabel('x'); ylabel('y');
            cmap;
            caxis([-6 6]);
            shading INTERP;
            title(sprintf('t=%f',t));
            axis([0 Model.M.Lx, 0 Model.M.Ly]);
            hold on;
            %query points
            plot(xq, yq,'*','markersize',10);
            hold off;
            colorbar;
            drawnow;
        end
    end
end
toc;

%% compare the solutions
for i = 1:2:length(rq); plot(t_a, p_a(i,:) + 5*rq(i),'k-'); hold on;end
for i = 1:2:length(rq); plot(t_n, p_n(i,:) + 5*rq(i),'r--'); hold on;end
hold off;
xlabel('time');
ylabel('pressure');
title('Red: analytical (FFT); Black: numerical');
set(gca,'fontsize',18);
