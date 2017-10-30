function [g, op] = grids_frac1d(nx, Lx, order, operator_type)
% construct grid and operators for 1d fracture, grouped by unknowns, p, u
%  weak b.c. treatment is used in this code.
%
% p: on staggered grid (m)
% u: on standard grid (p)

if nargin < 4
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

%% ***************************************grid and operators********************************************

% Now, generate grid and operators for each of the unknow field, p, u
switch operator_type
    case 'weak'
        [xp, xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_weak(order, nx, Lx/nx);
    case 'strong'
        [xp, xm, Pxp, Pxm, Qxp, Qxm] = sbp_staggered_strong(order,nx,Lx/nx,true);
end
xp = xp'; xm = xm';

% Difference operators
Dxp     = inv(Pxp)*Qxp;
Dxm    = inv(Pxm)*Qxm;

nxp = length(xp);
nxm = length(xm);

hx = Lx/nx;
%% p grid and op. staggered grid.
g.p.x = xp;
g.p.hx = hx;
g.p.nx = nxp;

% p op:
op.p.Dx    = Dxp;% 1D.
op.p.Px    = Pxp; % 1D
op.p.restrictions = restrictions(nxp); % restriction operator.

%% v grid and op. standard grid.
g.v.x = xm;
g.v.hx = hx;
g.v.nx = nxm;

% v op:
op.v.Dx    = Dxm;% 1D.
op.v.Px    = Pxm;% 1D.
op.v.restrictions = restrictions(nxm); 
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
