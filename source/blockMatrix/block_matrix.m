function Z = block_matrix(rows,columns,NZMAX)
% Z = block_matrix(nx_list,ny_list)
% Allocates a block sparse matrix that stores nnz non-zero entries. 
% The matrix contains only zeros after allocation.
%
% Input:
%    rows: A vector describing the number of elements in each block row the
%          block matrix has.
% columns: A vector describing the number of elements in each block column the
%          block matrix has.
%   NZMAX: The number of elements the block matrix can hold.
%
% Output:
%       Z: Sparse matrix (all zeros)
%
% Example:  
% >> Z = block_matrix([2 3],[2 3])
% >> size(Z)
% 
% ans =
% 
%      5     5

if nargin < 3
  NZMAX = 1;
end

m = sum(rows);
n = sum(columns);

Z = spalloc(m,n,NZMAX);
