function [d, V] = eigs_A(A, k, sigma, opt)
% this function take sparse matrix A and calculate k eigenvalues most near
% sigma with option specified by opt.
eps = 1e-4;
if nargin<4
    [V, D, ~] = eigs(A, k, sigma);
else
    [V, D, ~] = eigs(A, k, sigma, opt);
end

d=diag(D);
V=V(:,imag(d)>eps);
d = d(imag(d)>eps);

% sort the eigen-values in terms of the imaginary part in ascending order
[~, I] = sort(imag(d),1);

d=d(I);
V = V(:,I);% sort eigen vectors;

end
