function R = ark_stability_function(F,S,dt,eA,iA,I,b,e)
% Compute the stability function of a problem discretized in time
% with the additive Runge-kutta scheme
%
% Input:
% F : Implicit part
% S : Explicit part
%
% Output:
% R (scalar), if R > 1 then the scheme is unstable

if nargin < 4
  [eA iA I b e] = ark4_matrices();
end
      I2 = speye(size(F));
      Q = dt*F + dt*S;
      M = kron(I,I2) - kron(iA,dt*F) - kron(eA,dt*S);
      Z = kron(I,I2) + kron(e*b',Q)*inv(M);  
      R = max(abs(eig(full(Z)))); 
end
