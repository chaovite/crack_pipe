function [C,L,U,p,q,Ci,B_] = imex_pack(B,dt)
% y = imex_pack(B,dt,alpha)
% Packs the matrix B by removing all rows that are zero and then computes the LU
% decomposition of A = (I - dt*B_) where B_ is the compressed B. 
%
% Input arguments: 
% B  : The matrix to factorize
% dt : Scaling in the factorization (I - dt*B). The scaling coefficient must be
%      greater than zero.
% 
% Output arguments:
%   C : matrix that requires no decompostion
% L,U : Decomposition of A
% p,q : pivot vectors
% 
% 
% Split the matrix B into two parts, 
% a known part and an unknown part
% [ C | B_ ]*|u| = x 
%           |v|, 
% reduces to
% Bv = x - C*u, where u is known
%
% The solution to the problem Ax = b, where A = I - dt*B_ is given by
% LUx = b,
% set y = Ux,
% solve Ly = b
% solve Ux = y
% 
% y = L\b(p);
% x(q) = U\y;
% and then add -C*u to x

% Get non-zero rows
[Ci,cols] = find(B ~=0);
Ci = int32(unique(Ci));
r = false(size(B,1),1);
r(Ci) = true; 
rows = r;
B_ = B;
B_(~rows,:) = [];

% Remove B_ indices from B
C = B_;
C(:,r~=0) = 0;

% Remove zero columns from Bhat
B_(:,~rows) = [];

% Compute LU 
I = speye(size(B_));
[L,U,p,q]=lu(sparse(I-dt*B_),'vector');

p = int32(p);
q = int32(q);
end




