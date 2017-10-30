function [K, K_inv] = KernelBEM2D(x, G, nu)
%[K, K_inv] = KernelBEM2D(x, G, nu)
%   KernelBEM2d computes the BEM kernel K such that p  = K*w.
% Input:
%  x    : grid points coordinate.
%  G: rock shear modulus in Pa.
%  nu: Poisson's ratio
%
% Reference:
% Uenishi and Rice, 2003. Universal nucleation length for slip-weakening
% rupture instability under normal fault loading. JGR. 108-B1, 2042.

dx = x(2) - x(1);
N = length(x);
K = zeros(N, N);
G_star = G/(1 - nu);
for i = 1: N
    for j = 1: N
        K(i, j) = -G_star/(2*pi*dx*((i - j)^2 - 0.25));
    end
end
K = sparse(K);
K_inv = inv(K);

end

