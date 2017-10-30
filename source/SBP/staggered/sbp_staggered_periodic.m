function [Dp Dm Pm Pp xp xm] = sbp_staggered_periodic(acc,m,dx)
% [Dp Dm Pm Pp xp xm] = sbp_staggered_periodic(acc,m,dx)
% Periodic first derivative operator on staggered grids.
%
% Input arguments:
%   acc: Order of accuracy                                             (Even integer)
%     m: Number of grid points                
%    dx: Grid spacing 
%
% Output arguments:
%   Dp: Periodic first derivative operator on the nodal grid          (Sparse matrix m x m)
%   Dm: Periodic first derivative operator on the cell-centered grid  (Sparse matrix m x m)
%   Pp: Quadrature rule on the nodal grid                             (Sparse matrix m x m)
%   Pm: Quadrature rule on the cell-centered grid                     (Sparse matrix m x m)
%   xp: Nodal grid                                                    (Vector of size m)
%   xm: cell-centered grid                                            (Vector of size m)

  if nargin < 4
    test = true;
  end

  assert(mod(acc,2)==0,'The order of accuracy must be even');
  
  w = acc;

  % Construct the central finite difference stencil approximation and apply 
  % it for each row
  d  = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))'/dx;  
  Dp = spdiags(repmat(-d(end:-1:1)',[m 1]), -(w/2-1):w/2, m, m);  

  % The previous calculations needs to be adjusted near the boundaries
  % Since the operator should be periodic, we apply a circular shift
  % to each of the affected rows in D1.

  s       = full(Dp(acc/2,1:w));
  es      = spalloc(1,m,w);
  es(1:w) = s;
  for i = 1:acc/2
    Dp(i,:) = circshift(es,[0 (i-acc/2)]);
  end
  for i = 1:acc/2
    g             = circshift(es,[0 (i-acc/2-1)]);
    Dp(end-i+1,:) = -g(end:-1:1);
  end
  Dm = -Dp';

  % Quadrature rules
  Pp = dx*speye(m);
  Pm = dx*speye(m);

  % Grids
  xp = [0:m-1]'*dx;
  xm = ([0:m-1]-1/2)'*dx;

end                         
