function g = stretched_grid(grid_type,order,operator_type,nx,ny,Lx,Ly,r_g,r_bl,truncate)

% Construct grids with np = n + 1, and nm = n + 2 grid points


order = min(order,6);

[xp,xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_weak(order,nx,Lx/nx);
[yp, ym, Pyp, Pym, Qyp, Qym] = sbp_staggered_weak(order,ny,Ly/ny);
Dxp = inv(Pxp)*Qxp;
Dyp = inv(Pyp)*Qyp;
Dxm = inv(Pxm)*Qxm;
Dym = inv(Pym)*Qym;

g.hx = xp(2) - xp(1);
g.hy = yp(2) - yp(1);
g.nx = nx;
g.ny = ny;

yp = boundary_layer_thickness(yp/Ly, r_g, r_bl);
ym = boundary_layer_thickness(ym/Ly, r_g, r_bl);

xp = xp'; xm = xm';
yp = Ly*yp'; ym = Ly*ym';

nxp = length(xp);
nyp = length(yp);
nxm = length(xm);
nym = length(ym);

Jxp = spdiags(Dxp*xm,0, nxp, nxp);
Jxm = spdiags(Dxm*xp,0, nxm, nxm);
Jyp = spdiags(Dyp*ym,0, nyp, nyp);
Jym = spdiags(Dym*yp,0, nym, nym);

Jxpi  = inv(Jxp);
Jxmi  = inv(Jxm);
Jypi  = inv(Jyp);
Jymi  = inv(Jym);

Jxm  = trunc_mat( Jxm,truncate);
Jxmi = trunc_mat(Jxmi,truncate);
Jym  = trunc_mat( Jym,truncate);
Jymi = trunc_mat(Jymi,truncate);
xm   = trunc_vec(xm,truncate);
ym   = trunc_vec(ym,truncate);

grid_types = select_grid(grid_type);
if grid_types.is_p
  g.n1      = length(xp);
  g.n2      = length(yp);
  g.x       = xp;
  g.y       = yp;
  g.Jx      = Jxp;
  g.Jy      = Jyp;
  g.Jxi     = Jxpi;
  g.Jyi     = Jypi;
end

if grid_types.is_m
    g.n1     = length(xm);
    g.x      = xm;
    g.n2     = length(ym);
    g.y      = ym;
  g.Jx      = Jxm;
  g.Jy      = Jym;
  g.Jxi     = Jxmi;
  g.Jyi     = Jymi;
end

if grid_types.is_pp
    g.n1      = length(xp);
    g.n2      = length(yp);
    g.x       = xp;
    g.y       = yp;
    [g.X g.Y] = meshgrid(xp,yp);
    g.Jx      = Jxp;
    g.Jy      = Jyp;
    g.Jxi     = Jxpi;
    g.Jyi     = Jypi;
end

if grid_types.is_mm
    g.n1      = length(xm);
    g.n2      = length(ym);
    g.x       = xm;
    g.y       = ym;
    [g.X g.Y] = meshgrid(xm,ym);
    g.Jx      = Jxm;
    g.Jy      = Jym;
    g.Jxi     = Jxmi;
    g.Jyi     = Jymi;
end

if grid_types.is_pm
    g.n1      = length(xp);
    g.n2      = length(ym);
    g.x       = xp;
    g.y       = ym;
    [g.X g.Y] = meshgrid(xp,ym);
    g.Jx      = Jxp;
    g.Jy      = Jym;
    g.Jxi     = Jxpi;
    g.Jyi     = Jymi;
end

if grid_types.is_mp
    g.n1       = length(xm);
    g.n2       = length(yp);
    g.x        = xm;
    g.y        = yp;
    [g.X, g.Y] = meshgrid(xm,yp);
    g.Jx      = Jxm;
    g.Jy      = Jyp;
    g.Jxi     = Jxmi;
    g.Jyi     = Jypi;
end


g.vec = @(u) reshape(u,g.n1*g.n2,1);
g.grd = @(u) reshape(u,g.n2,g.n1);

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

