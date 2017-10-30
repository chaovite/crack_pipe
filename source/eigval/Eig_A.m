function [d, V] = Eig_A(A, eps, Q)
%[D, V]=EigA(A, n) returns eigenvalues in the order of descending period. 
%Input:
%       A: the matrix on which eigen value decomposition is to be performed
%       Q: threshold quality factor.
%output:
%       D: a vector contain the eigen value of the first n modes with
%       longest period.
%       V: a matrix, each column contains the corresponding eigenvector
%           given each eigenvalue in D

if issparse(A)
    A=full(A);
end

if nargin<3
      Q = -inf;
end

if nargin<2
      eps = 1e-3;
end

% seek all the eigen values of A;
[V, D] = eig(A);
d=diag(D);
ind = imag(d)>eps & abs(imag(d)./(2*real(d)))>Q;
V=V(:,ind);
d = d(ind);

% sort the eigen-values in terms of the imaginary part in ascending order
[~, I] = sort(imag(d),1);

d=d(I);
V = V(:,I);% sort eigen vectors;

end

