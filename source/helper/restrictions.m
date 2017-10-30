function E = restrictions(nx,ny,nz)
% E = restriction(nx)
% E = restriction(nx,ny)
% E = restriction(nx,ny,nz)
%
% Constructs sparse matrices that can be used to extract values along the
% boundaries
% 
% Input:
% nx,ny,nz : Number of grid points
%
% Output:s
% E  : struct containing vectors and matrices that extract the values along the
% a boundary by computing
% uX0 = E.X0*u, 
% or puts back values
% u^ = E.X0^T*uX0.

% 1D case
if nargin == 1
  E.e0 = spalloc(nx,1,1);
  E.e0(1) = 1;
  E.en = spalloc(nx,1,1);
  E.en(end) = 1;
end

% 2D case
if nargin == 2
  Ix = speye(nx);
  Iy = speye(ny);
  ex0 = spalloc(nx,1,1);
  ex0(1) = 1;
  exn = spalloc(nx,1,1);
  exn(end) = 1;
  ey0 = spalloc(ny,1,1);
  ey0(1) = 1;
  eyn = spalloc(ny,1,1);
  eyn(end) = 1;

  E.X0 = kron(ex0',Iy);
  E.XN = kron(exn',Iy);
  E.Y0 = kron(Ix,ey0');
  E.YN = kron(Ix,eyn');
  
  E.ex0 = ex0;
  E.exn = exn;
  E.ey0 = ey0;
  E.eyn = eyn;
  
end

% 3D case
if nargin == 3
  Ix = speye(nx);
  Iy = speye(ny);
  Iz = speye(nz);
  ex0 = spalloc(nx,1,1);
  ex0(1) = 1;
  exn = spalloc(nx,1,1);
  exn(end) = 1;
  ey0 = spalloc(ny,1,1);
  ey0(1) = 1;
  eyn = spalloc(ny,1,1);
  eyn(end) = 1;
  ez0 = spalloc(nz,1,1);
  ez0(1) = 1;
  ezn = spalloc(nz,1,1);
  ezn(end) = 1;
  E.X0 = kron(kron(ex0',Iy), Iz);
  E.XN = kron(kron(exn',Iy), Iz);
  E.Y0 = kron(kron(Ix,ey0'),Iz);
  E.YN = kron(kron(Ix,eyn'),Iz);
  E.Z0 = kron(kron(Ix,Iy),ez0');
  E.ZN = kron(kron(Ix,Iy),ezn');
  E.ex0 = ex0;
  E.exn = exn;
  E.ey0 = ey0;
  E.eyn = eyn;
  E.ez0 = ez0;
  E.ezn = ezn;
end

end
