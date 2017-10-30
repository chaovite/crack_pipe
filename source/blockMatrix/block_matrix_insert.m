function B = block_matrix_insert(B,rows,columns,i,j,block)
% Z = block_matrix_insert(B,rows,columns,block)
% Inserts a submatrix into a block matrix at block position (i,j)
%
% Input:
%       B: Block matrix (see block_matrix.m for a description).
%    rows: Number of elements in each block row of the block matrix.
% columns: Number of elements in each block row of the block matrix.
%     i,j: Block indices denoting where the submatrix should be placed.
%   block: The submatrix to insert into the block matrix. 
%          The size of the submatrix must match the number given in rows(i)
%          and columns(j).
%
% Example:
% >> Z = block_matrix([2 3],[2 3])
% >> B = block_matrix_insert(Z,[2 3],[2 3],1,2,ones(2,3));
% >> B
% 
% B =
% 
%    (1,3)        1
%    (2,3)        1
%    (1,4)        1
%    (2,4)        1
%    (1,5)        1
%    (2,5)        1
%>> full(B)
%
%ans =
%
%     0     0     1     1     1
%     0     0     1     1     1
%     0     0     0     0     0
%     0     0     0     0     0
%     0     0     0     0     0



assert(i <= length(rows),          'Index out of bounds.');
assert(j <= length(columns),       'Index out of bounds.');
assert(size(block,1) == rows(i),   'Dimension mismatch');
assert(size(block,2) == columns(j),'Dimension mismatch');

r = [1 rows];
c = [1 columns];
row_offset = sum(r(1:i));
col_offset = sum(c(1:j));
B(row_offset:row_offset+r(i+1)-1,col_offset:col_offset+c(j+1)-1) = block;




