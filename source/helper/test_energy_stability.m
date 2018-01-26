function [is_stable, is_energy_stable,eig_s,eig_es] = ...
  test_energy_stability(A,H,J,halt,tol)
% [is_stable, is_energy_stable,eig_s,eig_es] = test_energy_stability(A,H,J,halt,tol)
%
% Input arguments:
%               A: Spatial discretization (Sparse matrix kN x kN)
%               H: Energy matrix          (Sparse matrix kN x kN)
%
% Optional arguments:
%               J: Jacobian               (Sparse matrix N x N)
%            halt: Halt execution if true (Default: False)
%       tolerance: Minimum tolerance used to test that 
%                  real part of the eigenvalues is non-positive.
%                                         (Default: 1e-8)
%
% Output arguments:
% is_stable,is_energy_stable: Flags indicating whether approximation is stable
%                             or not.
%               eig_s,eig_es: Eigenvalues of the spatial discretization and
%                             corresponding energy.
  if (nargin < 3)
    J = 1;
  end

  if (nargin < 4)
    halt = false;
  end

  if (nargin < 5)
    tol = 1e-8;
  end

  if size(J,1) > 1
    kN = size(A,1);
    N = size(J,1);
    J = kron(speye(kN/N),J); 
  end

  eig_s = eig(full(A));
  is_stable = max(real(eig_s)) < tol;
                   
  eig_es = eig(full(J*H*A + A'*H*J));
  is_energy_stable = max(real(eig_es)) < tol;

  if (halt)
    assert(is_stable,['Not stable!  max(real(lambda)) = ' ...
           num2str(max(real(eig_s)))]);
    assert(is_energy_stable,['Not energy stable!  max(real(lambda)) = ' ...
           num2str(max(real(eig_es)))]);
  end
end
