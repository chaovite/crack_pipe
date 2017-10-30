function [g, op] = grids_pipe1d(M)
%
% construct the grid for 1d pipe with multiple material interfaces
% only 1D grid in z direction is constructed.
% only first derivative.
%

% All the checks are done.
nL = length(M.L);   % number pipe sections.
n_int = nL - 1;        % number of interfaces.
dz   = M.L./M.nz;   % grid spacing, dimension [1, nL]

% need to output SAT.
dim   = M.nz + 1;
Hinv  = block_matrix(dim, dim, 0);
Dz     = block_matrix(dim, dim, 0);
z       = block_matrix(dim, 1, 0);

tmp   = 0;
for i = 1: nL
    ni   = M.nz(i);
    dzi = dz(i);
    Li   = M.L(i);
    [Hinv_i, Dz_i] = sbp_sparse(M.order, ni, dzi);
    Hinv = block_matrix_insert(Hinv, dim, dim, i, i, Hinv_i);
    Dz = block_matrix_insert(Dz, dim, dim, i, i, Dz_i);
    zi  = [0:ni]'*dzi + tmp;
    z   = block_matrix_insert(z, dim, 1, i, 1, zi);
    tmp = tmp + Li;
end

H  = inv(Hinv);
nz = sum(dim);

% construct the grids and grid properties.
g.dz = dz; 
g.z   = z; % Construct z grid
g.nz = nz; % total number grid points in the pipe.
g.dim = dim;

% grid properties, store either as a diagonal matrix.
list = {'R','rho','c', 'mu','S', 'K'};
ind_tmp = [0, cumsum(dim)];
for i = 1: length(list)
    field = list{i};
    if length(M.(field)) == 1
        g.(field) = spdiags(M.(field)*ones(nz, 1), 0, nz, nz);
    else
        field_i = ones(nz, 1);
        for j = 1:nL
            start = ind_tmp(j) + 1;
            fini   = ind_tmp( j + 1);  
             field_i(start:fini) = M.(field)(j); 
             g.(field) = spdiags(field_i, 0, nz, nz);
        end
    end
end

% construct the operators
% one dimensional operators
op.Dz     = Dz;
op.Hinv  = Hinv;
op.H       = H;

op.bc      = M.bc;
op.interfaces = M.interfaces';

% operators at the boundaries and interfaces.
% restrictions and SAT terms
e0 = spalloc(g.nz,1,1);
en = spalloc(g.nz,1,1);
e0(1)   = 1;
en(end) = 1;
op.bc.bt.e = e0;
op.bc.bt.SAT = g.c(1,1)/H(1,1);
op.bc.tp.e = en;
op.bc.tp.SAT = g.c(end,end)/H(end,end);

% restriction operators on the interface.
for i = 1: n_int
    ind_p = sum(dim(1:i)) + 1;
    ind_m = sum(dim(1:i));
    % interface minus side.
    op.interfaces{i}.ep      = spalloc(g.nz,1,1);
    % interface plus side.
    op.interfaces{i}.em     = spalloc(g.nz,1,1); 
    op.interfaces{i}.ep(ind_p)  = 1;
    op.interfaces{i}.em(ind_m) = 1;
    % interface index.
    op.interfaces{i}.indm           = ind_m;
    op.interfaces{i}.indp            = ind_p;
    % SAT at the plus and minus side.
    op.interfaces{i}.SATp  = g.c(ind_p, ind_p)/H(ind_p, ind_p);
    op.interfaces{i}.SATm = g.c(ind_m, ind_m)/H(ind_m, ind_m);
end

end

