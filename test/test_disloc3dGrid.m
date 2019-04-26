%test disloc3dGrid.
clear;close all;
source_dir = '../source';
addpath(genpath(source_dir));
%%
% fault length.
xc = 0;
yc = 0;
zc = 600;

L   = 100;
W  = 100;

nx  = 20;
ny  = 20;

dx  = W/nx;
dy = L/ny;

x    = [1/2:1: (nx-1/2)]*dx -1/2*dx;
y    = [1/2:1: (ny-1/2)]*dy -1/2*dy;

[X, Y] = meshgrid(x, y);

mu = 10e9;
nu  = 0.25;

% uniform opening.
w    =  1*ones(size(X));

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
dsl_grid = disloc3dGrid(xc, yc, zc, strike, dip, X, Y, mu, nu);
toc;
% compare displacement.
U_w   = dsl_grid.eval_disp_w(w, xq, yq);
U_p    = dsl_grid.eval_disp_p(dsl_grid.K*w(:), xq, yq);
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




