function [g, op] = grids_conduit(nr, R, nz, L, order)
% construct the grid for the conduit. 
% grid is straggered in r direction but standard in z direction.
dz = L/nz;
dr = R/nr;
truncate = 1;
[rp,rm,Pp,Pm,Qp,Qm] = sbp_staggered_strong(order,nr, dr, truncate);
[Pz_inv, Dz] = sbp_sparse(order,nz,dz);
Pz = inv(Pz_inv);

Rp = spdiags(rp,0,nr+1,nr+1);
Rm = spdiags(rm,0,nr,nr);

% Assemble second derivative in r direction.
D2 = inv(Rm)*inv(Pm)*Qm*Rp*inv(Pp)*Qp;

e0 = spalloc(nz+1,1,1);
en = spalloc(nz+1,1,1);

e0(1)   = 1;
en(end) = 1;

er = ones([nr,1]);
ez = ones([nz+1,1]);

Iz = speye(nz+1);
Ir = speye(nr);

% width averaging operators.
W1 = (2/R^2)*er'*Rm*Pm; % Cylindrical width averaged operator that operators on just one row of v, dimension [1, nr]
W2 = kron(Iz,W1); % Operates on the entire v, dimension [(nz+1)*nr]

% construct the grids
g.dz = dz;
g.dr = dr;
g.z = dz*[0:nz]'; % Construct z grid
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
op.er = er;
op.ez = ez;
op.Iz = Iz;
op.Ir = Ir;

end

