function g = operators(grid_type,order,geom,operator_type)

g.grid_type     = grid_type;
g.operator_type = operator_type;

nx = geom.nx;
ny = geom.ny;

hx = geom.hx;
hy = geom.hy;

nxp = nx + 1;
nxm = nx + 2;
nyp = ny + 1;
nym = ny + 2;

switch operator_type
    case 'weak'
        [xp xm Pxp Pxm Qxp Qxm] = sbp_staggered_weak(order,nx,hx);
        [yp ym Pyp Pym Qyp Qym] = sbp_staggered_weak(order,ny,hy);
    case 'strong'
        [xp xm Pxp Pxm Qxp Qxm] = sbp_staggered_strong(order,nx,hx,true);
        [yp ym Pyp Pym Qyp Qym] = sbp_staggered_strong(order,ny,hy,true);
        nxm = nxm - 2;
        nym = nym - 2;
end

% Difference operators
Dxp      = inv(Pxp)*Qxp;
Dxm      = inv(Pxm)*Qxm;
Dyp      = inv(Pyp)*Qyp;
Dym      = inv(Pym)*Qym;

% Identity matrices
Ixp = speye(nxp);
Iyp = speye(nyp);
Ixm = speye(nxm);
Iym = speye(nym);

grid_types = select_grid(grid_type);

if grid_types.is_p
    
    % Number of grid points
    g.n1 = nxp;
    g.n2 = nyp;
    
    % Difference operators
    g.Dx = Dxp;
    g.Dy = Dyp;
    
    % Quadrature operators
    g.Px = Pxp;
    g.Py = Pyp;
    
    g.Ix = Ixp;
    g.Iy = Iyp;
end

if grid_types.is_m
    
    % Number of grid points
    g.n1 = nxm;
    g.n2 = nym;
    
    % Difference operators
    g.Dx = Dxm;
    g.Dy = Dym;
    
    % Quadrature operators
    g.Px = Pxm;
    g.Py = Pym;

    % Identity operators
    g.I = Ixm;
end




if grid_types.is_pp
    % Number of grid points
    g.n1 = nxp;
    g.n2 = nyp;
    
    % Difference operators
    g.Dx = kron(Dxp,Iyp);
    g.Dy = kron(Ixp,Dyp);
    
    % 1D version
    g.Dx1 = Dxp;
    g.Dy1 = Dyp;
    
    % Quadrature operators
    g.Px = Pxp;
    g.Py = Pyp;
    
    % Restriction operators
    g.restrictions = restriction_operators(nxp,nyp);
end

if grid_types.is_mm
    % Number of grid points
    g.n1 = nxm;
    g.n2 = nym;
    
    % Difference operators
    g.Dx = kron(Dxm,Iym);
    g.Dy = kron(Ixm,Dym);
    
    % 1D version
    g.Dx1 = Dxm;
    g.Dy1 = Dym;
    
    % Quadrature operators
    g.Px = Pxm;
    g.Py = Pym;
    
    % Restriction operators
    g.restrictions = restriction_operators(nxm,nym);
end

if grid_types.is_pm
    % Number of grid points
    g.n1 = nxp;
    g.n2 = nym;
    
    % Difference operators
    g.Dx = kron(Dxp,Iym);
    g.Dy = kron(Ixp,Dym);
    
    % 1D version
    g.Dx1 = Dxp;
    g.Dy1 = Dym;
    
    % Quadrature operators
    g.Px = Pxp;
    g.Py = Pym;
    
    % Restriction operators
    g.restrictions = restriction_operators(nxp,nym);
end

if grid_types.is_mp
    % Number of grid points
    g.n1 = nxm;
    g.n2 = nyp;
    
    % Difference operators
    g.Dx = kron(Dxm,Iyp);
    g.Dy = kron(Ixm,Dyp);
    
    % 1D version
    g.Dx1 = Dxm;
    g.Dy1 = Dyp;
    
    % Quadrature operators
    g.Px = Pxm;
    g.Py = Pyp;
    
    % Restriction operators
    g.restrictions = restriction_operators(nxm,nyp);
end

% Grid spacing
g.hx = hx;
g.hy = hy;
g.nx = nx;
g.ny = ny;

end
