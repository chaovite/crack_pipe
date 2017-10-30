function g = grids(grid_type,nx,ny,operator_type,Lx,Ly)

if nargin < 4
    operator_type = 'weak';
end

if nargin < 5
    Lx = 1;
end

if nargin < 6
    Ly = 1;
end

grid_types = select_grid(grid_type);

[xp xm] = grids_1d(nx,operator_type);
[yp ym] = grids_1d(ny,operator_type);


xp = Lx*xp; xm = Lx*xm;
yp = Ly*yp; ym = Ly*ym;


g.hx = xp(2) - xp(1);
g.hy = yp(2) - yp(1);

g.nx = nx;
g.ny = ny;

if grid_types.is_p
    g.n1      = length(xp);
    g.n2      = length(yp);
    g.x        = xp;
    g.y        = yp;
end

if grid_types.is_m
    g.n1     = length(xm);
    g.x       = xm;
    g.n2     = length(ym);
    g.y       = ym;
end


if grid_types.is_pp
    g.n1      = length(xp);
    g.n2      = length(yp);
    g.x       = xp;
    g.y       = yp;
    [g.X g.Y] = meshgrid(xp,yp);
end

if grid_types.is_mm
    g.n1      = length(xm);
    g.n2      = length(ym);
    g.x       = xm;
    g.y       = ym;
    [g.X g.Y] = meshgrid(xm,ym);
end

if grid_types.is_pm
    g.n1      = length(xp);
    g.n2      = length(ym);
    g.x       = xp;
    g.y       = ym;
    [g.X g.Y] = meshgrid(xp,ym);
end

if grid_types.is_mp
    g.n1      = length(xm);
    g.n2      = length(yp);
    g.x       = xm;
    g.y       = yp;
    [g.X, g.Y] = meshgrid(xm,yp);
end

g.vec = @(u) reshape(u,g.n1*g.n2,1);
g.grd = @(u) reshape(u,g.n2,g.n1);

end
