function E = restriction_operators(nx,ny)
% E = restriction_operators(nx)
% E = restriction_operators(nx,ny)
% E = restriction_operators(nx,ny,nz)
%
% Constructs sparse matrices that can be used to extract values along the
% boundaries
% 
% Input:
% nx,ny,nz : Number of grid points
%
% Output:
% E  : struct containing vectors and matrices that extract the values along the
% a boundary by computing
% uw = E.W*u, 
% or puts back values
% u^ = E.W^T*u.

% 1D case
if nargin == 1
  E.e0 = spalloc(nx,1,1);
  E.e0(1) = 1;
  E.en = spalloc(nx,1,1);
  E.en(end) = 1;
end

% 2D case
if nargin == 2
  exp0 = spalloc(nx,1,1);
  exp0(1) = 1;
  exm0 = spalloc(nx,1,1);
  exm0(1) = 1;

  eyp0 = spalloc(ny,1,1);
  eyp0(1) = 1;
  eym0 = spalloc(ny,1,1);
  eym0(1) = 1;

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

  E.W = kron(ex0',Iy);
  E.E = kron(exn',Iy);
  E.N = kron(Ix,eyn');
  E.S = kron(Ix,ey0');

end
end
