function [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_weak(order,n,h,test,x)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_half(order,n,h)
%
% Construct SBP staggered grid operators
% that satisfy the property (Q+)_00  = (Q_-)_00,  (H+)_00 = (H-)_00
%
% Input:
% order : Order of accuracy
% n     : Number of grid points n+1 (nodal grid) n+2 (cell-centered grid)
% h     : Grid spacing

%
% Output
% xp,xm         : Grid vectors xp (nodal grid) xm (cell-centered grid)
% Pp,Pm,Qp,Qm,  : Staggered grid operators

if nargin < 4
  test = false;
end
if nargin < 5
  x = [];
end

    switch order
      case 2
      [xp,xm,Pp,Pm,Qp,Qm]   = sbp_staggered_weak_2nd(n,h,test);
      case 4
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_weak_4th(n,h,test,x);
      case 6
        [xp,xm,Pp,Pm,Qp,Qm] = sbp_staggered_weak_6th(n,h,test,x);
      case 8
       error('SBP staggered grid operator not implemented');   
        %[xp,xm,Pp,Pm,Qp,Qm] = sbp_half_8th(n,h,test);
      otherwise
       error('SBP staggered grid operator not implemented');   
      end

xp = xp';
xm = xm';


end
