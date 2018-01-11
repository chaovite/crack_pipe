function [g, op] = grids_frac2dxy(nx, ny, Lx, Ly, order, operator_type)
% construct grid and operators for 2d fracture, grouped by unknowns, p, vx, vy
%
% vx, vy are width averaged velocities, assume Poiseuille flow for viscosity.
%
% strong boundary conditions are imposed by default.
%
% p: x,y standard grid (pp)
% vx: x staggered, y standard.  (mp)
% vy: x standard,   y staggered (pm)
%

if nargin < 6
    operator_type = 'strong';
end

%% ***************************************grid and operators********************************************

% Now, generate grid and operators for each of the unknow field, p, vx, vy
switch operator_type
    case 'weak'
        [xp, xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_weak(order,nx,Lx/nx);
        [yp, ym, Pyp, Pym, Qyp, Qym] = sbp_staggered_weak(order,ny,Ly/ny);
    case 'strong'
        [xp, xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_strong(order,nx,Lx/nx,true);
        [yp, ym, Pyp, Pym, Qyp, Qym] = sbp_staggered_strong(order,ny,Ly/ny,true);
end
xp = xp'; xm = xm';
yp = yp'; ym = ym';

% Difference operators
Dxp      = inv(Pxp)*Qxp;
Dxm      = inv(Pxm)*Qxm;
Dyp      = inv(Pyp)*Qyp;
Dym      = inv(Pym)*Qym;

nxp = length(xp);
nxm = length(xm);
nyp = length(yp);
nym = length(ym);

% Identity matrices
Ixp = speye(nxp);
Iyp = speye(nyp);
hx = Lx/nx;
hy = Ly/ny;
%% p grid:
g.p.x = xp;
g.p.y = yp;
g.p.hx = hx;
g.p.hy = hy;
g.p.nx = nxp;
g.p.ny = nyp;
[g.p.X, g.p.Y] =  meshgrid(xp,yp);
g.p.vec = @(u) reshape(u, g.p.ny*g.p.nx, 1);
g.p.grd = @(u) reshape(u, g.p.ny, g.p.nx);

% p op:
op.p.Dx1    = Dxm;% 1D.
op.p.Dy1    = Dym;
op.p.Dx2      = kron(Dxm, Iyp);% dp/dx for vx (mp), 2D
op.p.Dy2      = kron(Ixp,Dym); % dp/dy for vy (pm), 2D
op.p.Px1    = Pxp; % 1D
op.p.Py1    = Pyp; % 1D
op.p.Pxy2  = kron(Pxp, Pyp); % for energy norm dp/dt.
op.p.restrictions = restrictions(nxp, nyp); % restriction operator.

%% vx grid (x, y) and op.
g.vx.x = xm;
g.vx.y = yp;
g.vx.hx = hx;
g.vx.hy = hy;
g.vx.nx = nxm;
g.vx.ny = nyp;
[g.vx.X, g.vx.Y] =  meshgrid(xm,yp);
g.vx.vec = @(u) reshape(u, g.vx.ny*g.vx.nx, 1);
g.vx.grd = @(u) reshape(u, g.vx.ny, g.vx.nx);

% p op:
op.vx.Dx1    = Dxp;% 1D.
op.vx.Dx2    = kron(Dxp, Iyp);%dvx/dx for dp/dt
op.vx.Px1    = Pxm; % 1D
op.vx.Py1    = Pyp; % 1D
op.vx.Pxy2  = kron(Pxm, Pyp); % for energy norm dp/dt.
op.vx.restrictions = restrictions(nxm, nyp); 
%% vy grid (x, y) and op.
g.vy.x = xp;
g.vy.y = ym;
g.vy.hx = hx;
g.vy.hy = hy;
g.vy.nx = nxp;
g.vy.ny = nym;
[g.vy.X, g.vy.Y] =  meshgrid(xp,ym);
g.vy.vec = @(u) reshape(u, g.vy.ny*g.vy.nx, 1);
g.vy.grd = @(u) reshape(u, g.vy.ny, g.vy.nx);

% p op:
op.vy.Dy1    = Dyp;% 1D.
op.vy.Dy2    = kron(Ixp, Dyp);% dvy/dy for dp/dt. 2D.
op.vy.Px1    = Pxp; % 1D
op.vy.Py1    = Pym; % 1D
op.vy.Pxy2  = kron(Pxp, Pym); % for energy norm dp/dt.
op.vy.restrictions = restrictions(nxp, nym); 
end