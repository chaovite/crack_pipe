function [g, op] = grids_conduit_split(nr, R, nz, L, ind, order)
% construct the grid for the conduit with material interface splitting.
% ind: the grid point index to be 
% grid is straggered in r direction but standard in z direction.
switch order
    case 2
        min_nz = 4;
    case 4
        min_nz = 8;
    case 6
        min_nz = 12;
    case 8
        min_nz = 16;
    otherwise 
        error('>8 higher order discretization not implemented!');
end

dz = L/nz;
dr = R/nr;
truncate = 1;
[rp,rm,Pp,Pm,Qp,Qm] = sbp_staggered_strong(order,nr, dr, truncate);

nz1 = ind - 1;
nz2 = nz - nz1;

if min(nz1, nz2) < min_nz
    error(['both conduit sections need to have minimum ',...
        num2str(min_nz),' points for order ',num2str(order),'.']);
end

[Pz_inv1, Dz1] = sbp_sparse(order,nz1,dz);
[Pz_inv2, Dz2] = sbp_sparse(order,nz2,dz);
dim = [size(Pz_inv1, 1), size(Pz_inv2, 1)];

Pz_inv = block_matrix(dim, dim, 0);
Dz       = block_matrix(dim, dim, 0);
Pz_inv = block_matrix_insert(Pz_inv, dim, dim, 1, 1, Pz_inv1);
Pz_inv = block_matrix_insert(Pz_inv, dim, dim, 2, 2, Pz_inv2);
Pz        = inv(Pz_inv);

Dz = block_matrix_insert(Dz, dim, dim, 1, 1, Dz1);
Dz = block_matrix_insert(Dz, dim, dim, 2, 2, Dz2);

Rp = spdiags(rp,0,nr+1,nr+1);
Rm = spdiags(rm,0,nr,nr);

% Assemble second derivative in r direction.
D2 = inv(Rm)*inv(Pm)*Qm*Rp*inv(Pp)*Qp;

e0 = spalloc(nz+2,1,1);
en = spalloc(nz+2,1,1);

e0(1)   = 1;
en(end) = 1;

% the interface minus and plus side.
eIm = spalloc(nz+2,1,1);
eIp = spalloc(nz+2,1,1);
eIm(ind) = 1;
eIp(ind + 1) = 1;

er = ones([nr,1]);
ez = ones([nz+2,1]);

Iz = speye(nz+2);
Ir = speye(nr);

% width averaging operators.
W1 = (2/R^2)*er'*Rm*Pm; % Cylindrical width averaged operator that operators on just one row of v, dimension [1, nr]
W2 = kron(Iz,W1); % Operates on the entire v, dimension [(nz+2)*nr]

% construct the grids
g.dz = dz;
g.dr = dr;
g.z = dz*[0:ind-1, ind-1:nz]'; % Construct z grid
g.rp = rp;
g.rm = rm;
g.nz = length(g.z);
g.nr = length(rm);
[g.R, g.Z] = meshgrid(rm, g.z);

% construct the operators
% one dimensional operators
op.Dz = Dz;
op.Pz_inv = Pz_inv;
op.Pz   = Pz;
op.Rm = Rm;
op.Rp  = Rp; 
op.Pp= Pp;
op.Pm= Pm;
op.Qp= Qp;
op.Qm= Qm;
op.D2 = D2;
op.W1 = W1;
op.W2 = W2;
op.e0 = e0;
op.en = en;
op.eIm = eIm; % interface minus side.
op.eIp = eIp; % interface plus side.
op.er = er;
op.ez = ez;
op.Iz = Iz;
op.Ir = Ir;

end

