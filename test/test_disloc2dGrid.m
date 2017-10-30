% test disloc2dGrid
% script the validate the implementation of disloc2dGrid, it must reduce to
% the solution of KernelBEM2D when the length of out-of-plane dimension to
% large number. This should also reduce to solution given by disloc3dGrid
% after massaging the grid a little bit.

source_dir = '../source';
addpath(genpath(source_dir));
%%
% fault length.
xc = 0;
yc = 0;
zc = 600;

L   = 100;
W  = 100;

nx  = 200;
ny  = 1;

dx  = W/nx;
dy = L/ny;

x    = [1/2:1: (nx-1/2)]*dx -1/2*dx;
x = x';

mu = 10e9;
nu  = 0.25;

% uniform opening.
w    =  1*ones(size(x));

% observation points.
xq = linspace(0,1,20)*W;
yq = zeros(1, 20)+W;
zq = zeros(1, 20);

%% check for different strikes and dips.
strike = 0; dip = 90;
% reference point.
depth = zc + W/2*sind(dip);
east   = xc + W/2*cosd(dip)*cosd(strike);
north = yc - W/2*cosd(dip)*sind(strike);

tic;
dsl_grid = disloc2dGrid(xc, yc, zc, strike, dip, x, L, mu, nu);
toc;
% compare displacement.
U_w   = dsl_grid.eval_disp_w(w, xq, yq);
U_p    = dsl_grid.eval_disp_p(dsl_grid.K*w(:), xq, yq);
%%
% calculation from disloc3d.  mdl = [length width depth dip strike east north ss ds op]';
mdl = [L, W, depth, dip, strike, east, north, 0, 0, 1]';
obs =[xq; yq; zq];
[U, ~, ~, ~] = disloc3d(mdl,obs,mu,nu);
figure(1); plot(xq,U_p,'-',xq,U_w,xq, U,'*');xlabel('x'); ylabel('U');
legend({'x','y','z'});
title(sprintf('y=0, strike = %d, dip = %d',strike, dip));

U_w   = dsl_grid.eval_disp_w(w, yq, xq);
U_p    = dsl_grid.eval_disp_p(dsl_grid.K*w(:), yq, xq);

% calculation from disloc3d.  mdl = [length width depth dip strike east north ss ds op]';
mdl = [L, W, depth, dip, strike, east, north, 0, 0, 1]';
obs =[yq; xq; zq];
[U, ~, ~, ~] = disloc3d(mdl,obs,mu,nu);
figure(2); plot(xq, U_p,'-',xq,U_w,'--',xq, U,'*');xlabel('y'); ylabel('U');
legend({'x','y','z'});
title(sprintf('x=0, strike = %d, dip = %d',strike, dip));

%% compare the elastic DEM kernel with 
nx  = 200;
dx  = W/nx;
x    = [1/2:1: (nx-1/2)]*dx -1/2*dx;
x = x';
p    = 1e6;       % uniform pressure 1 MPa.
N = length(x);
p_vec = ones(length(x),1)*p; 
[K, K_inv] = KernelBEM2D(x, mu, nu);
dsl_grid   = disloc2dGrid(xc, yc, zc, strike, dip, x, L, mu, nu);
K_inv_dsl = dsl_grid.K_inv;

w_n = K_inv*p_vec; 
w_n_dsl = K_inv_dsl*p_vec;

% Analytical solution;
c    = W/2;
w_a = 2*c*p*(1-nu)*sqrt(1-((x-c)/c).^2)/mu;

% compare results:
figure(3);
plot(x, w_n, 'b--',x, w_n_dsl, 'r--',x, w_a,'k-','linew',2);
legend({'DDM2D','disloc3d','Analytical solution'},'Location','South');
xlabel('x along the crack (m)');
ylabel('Opening (m)');
title(['Test DDM2D, N = ',num2str(N)]);
set(gca, 'fontsize',20);
