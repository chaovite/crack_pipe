function [g, op] = grids_frac3d(nx,ny,nz, Lx, Ly, Lz, order, operator_type, r_g, r_bl)
% construct grid and operators for 3d fracture, grouped by unknowns, p, vx, vy, ux, uy.
% r_g, r_bl are grid stretching parameters in z direction.
%
% p: x,y standard grid, z, staggered grid. (ppm)
% vx: x staggered, y,z standard. (mpp)
% vy: y staggered, x,z standard. (pmp)
% ux: x staggered, y standard.(mp)
% uy: y staggered, x standard.(pm)

if nargin < 10
    % no grid stretching.
    r_g = 0.3;
    r_bl = 0.3;
end

if nargin < 8
    operator_type = 'weak';
end

switch operator_type
    case 'weak'
        truncate = 0;
    case 'strong'
        truncate = 1;
    otherwise
        error('Incorrect operator type.');
end

%% ******************************construct Jacobian in z direction*************************************

% weak operators are used to construct the Jacobian.
[zp, zm, Pzp, Pzm, Qzp, Qzm] = sbp_staggered_weak(order,nz,Lz/nz);
Dzp = inv(Pzp)*Qzp;
Dzm = inv(Pzm)*Qzm;

% stretch the grid only in z direction.
zp = boundary_layer_thickness(zp/Lz, r_g, r_bl);% stretched grid
zm = boundary_layer_thickness(zm/Lz, r_g, r_bl);% stretched grid.
zp = Lz*zp'; zm = Lz*zm';

nzp = length(zp);
nzm = length(zm);

% Jacobian. d_eta/dz. (eta is the stretched grid and z is the uniform grid)
Jzp = spdiags(Dzp*zm,0, nzp, nzp); 
Jzm = spdiags(Dzm*zp,0, nzm, nzm);

% Jacobian. dz/d_eta.
Jzpi  = inv(Jzp);
Jzmi  = inv(Jzm);

% if the operator type is 'strong', then truncate the end points.
Jzm  = trunc_mat(Jzm,truncate);
Jzmi = trunc_mat(Jzmi,truncate);
zm   = trunc_vec(zm,truncate);
zp   = trunc_vec(zp,truncate);
nzp = length(zp);
nzm = length(zm);

%% ***************************************grid and operators********************************************

% Now, generate grid and operators for each of the unknow field, p, vx, vy, ux, uy
switch operator_type
    case 'weak'
        [xp, xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_weak(order,nx,Lx/nx);
        [yp, ym, Pyp, Pym, Qyp, Qym] = sbp_staggered_weak(order,ny,Ly/ny);
        [~, ~, Pzp, Pzm, Qzp, Qzm] = sbp_staggered_weak(order,nz,Lz/nz);
    case 'strong'
        [xp, xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_strong(order,nx,Lx/nx,true);
        [yp, ym, Pyp, Pym, Qyp, Qym] = sbp_staggered_strong(order,ny,Ly/ny,true);
        [~, ~, Pzp, Pzm, Qzp, Qzm]   = sbp_staggered_strong(order,nz,Lz/nz, true);
end
xp = xp'; xm = xm';
yp = yp'; ym = ym';

% Difference operators
Dxp      = inv(Pxp)*Qxp;
Dxm      = inv(Pxm)*Qxm;
Dyp      = inv(Pyp)*Qyp;
Dym      = inv(Pym)*Qym;
Dzp      = inv(Pzp)*Qzp;
Dzm      = inv(Pzm)*Qzm;

nxp = length(xp);
nxm = length(xm);
nyp = length(yp);
nym = length(ym);

% Identity matrices
Ixp = speye(nxp);
Iyp = speye(nyp);
Izp = speye(nzp);
Ixm = speye(nxm);
Iym = speye(nym);
Izm = speye(nzm);

hx = Lx/nx;
hy = Ly/ny;
hz = Lz/nz;

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
op.p.Dx3    = kron(kron(Dxm, Iyp), Izm);% dp/dx for vx (mpm)
op.p.Dy3    = kron(kron(Ixp,Dym), Izm); % dp/dy for vy (pmm)
op.p.ez       = ones(nzm, 1); 
op.p.Ez      = kron(kron(Ixp, Iyp), op.p.ez); % extend p in z direction into 3D.
op.p.Px1      = Pxp; % 1D
op.p.Py1     = Pyp;
op.p.Pxy2  = kron(Pxp, Pyp); % for energy norm dp/dt.
op.p.restrictions = restrictions(nxp, nyp); % restriction operator.

%% vx, grid and op.
g.vx.x = xm;
g.vx.y = yp;
g.vx.z = zm;
g.vx.Jzp = Jzp;
g.vx.Jzm = Jzm;
g.vx.Jzpi = Jzpi;
g.vx.Jzmi = Jzmi;
g.vx.hx = hx;
g.vx.hy = hy;
g.vx.hz = hz;
g.vx.nx = nxm;
g.vx.ny = nyp;
g.vx.nz = nzm;
% store 1D grid and construct 3D when using it.
g.vx.mesh = @() meshgrid(xm, yp, zm);
% [g.vx.X, g.vx.Y, g.vx.Z] =  meshgrid(xm, yp, zm); this cost too much memory.
g.vx.vec = @(u) reshape(permute(u,[3,1,2]), g.vx.nz*g.vx.ny*g.vx.nx,1);
g.vx.grd = @(u) reshape(u, g.vx.nz, g.vx.ny, g.vx.nx);

% vx, op. Attention! Grid stretch in z direction!
op.vx.D2z1 = Jzmi*Dzm*Jzpi* Dzp;% D2z in 1D
op.vx.D2z3 = kron(kron(Ixm,Iyp), op.vx.D2z1); % D2z in 3D
op.vx.Wz1  = 1/Lz*ones(1, g.vx.nz)*Jzm*Pzm;% width-average operator in 1D
op.vx.Wz3  = kron(kron(Ixm, Iyp), op.vx.Wz1);% width-average operator in 3D
op.vx.Pxyz3 = kron(kron(Pxm,Pyp), Jzm*Pzm);% for the energy norm.
op.vx.restrictions = restrictions(nxm,nyp,nzm);

%% vy grid and op.
g.vy.x = xp;
g.vy.y = ym;
g.vy.z = zm;
g.vy.Jzp = Jzp;
g.vy.Jzm = Jzm;
g.vy.Jzpi = Jzpi;
g.vy.Jzmi = Jzmi;
g.vy.hx = hx;
g.vy.hy = hy;
g.vy.hz = hz;
g.vy.nx = nxp;
g.vy.ny = nym;
g.vy.nz = nzm;
g.vy.mesh = @() meshgrid(xp, ym, zm);
g.vy.vec = @(u) reshape(permute(u,[3,1,2]), g.vy.nz*g.vy.ny*g.vy.nx,1);
g.vy.grd = @(u) reshape(u, g.vy.nz, g.vy.ny, g.vy.nx);

op.vy.D2z1 = Jzmi*Dzm*Jzpi* Dzp;
op.vy.D2z3 = kron(kron(Ixp,Iym), op.vy.D2z1);
op.vy.Wz1  = 1/Lz*ones(1, g.vx.nz)*Jzm*Pzm;
op.vy.Wz3  = kron(kron(Ixp, Iym), op.vy.Wz1);
op.vy.Pxyz3 = kron(kron(Pxp,Pym), Jzm*Pzm);% for the energy norm.
op.vy.restrictions = restrictions(nxp,nym,nzm);
%% ux grid (x, y) and op.
g.ux.x = xm;
g.ux.y = yp;
g.ux.hx = hx;
g.ux.hy = hy;
g.ux.nx = nxm;
g.ux.ny = nyp;
[g.ux.X, g.ux.Y] =  meshgrid(xm,yp);
g.ux.vec = @(u) reshape(u, g.ux.ny*g.ux.nx, 1);
g.ux.grd = @(u) reshape(u, g.ux.ny, g.ux.nx);

% p op:
op.ux.Dx1    = Dxp;% 1D.
op.ux.Dx2    = kron(Dxp, Iyp);%dux/dx for dp/dt
op.ux.restrictions = restrictions(nxm, nyp); 
%% uy grid (x, y) and op.
g.uy.x = xp;
g.uy.y = ym;
g.uy.hx = hx;
g.uy.hy = hy;
g.uy.nx = nxp;
g.uy.ny = nym;
[g.uy.X, g.uy.Y] =  meshgrid(xp,ym);
g.uy.vec = @(u) reshape(u, g.uy.ny*g.uy.nx, 1);
g.uy.grd = @(u) reshape(u, g.uy.ny, g.uy.nx);

% p op:
op.uy.Dy1    = Dyp;% 1D.
op.uy.Dy2    = kron(Ixp, Dyp);% duy/dy for dp/dt. 2D.
op.uy.restrictions = restrictions(nxp, nym); 
end

function A = trunc_mat(A,truncate)
  if ~truncate
    return;
  end
  A = A(2:end-1,2:end-1);
end

function a = trunc_vec(a,truncate)
  if ~truncate
    return;
  end
  a = a(2:end-1);
end
