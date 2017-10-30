function [A,G,Fp,lambda,z,rm,ind, Ai, Ae] = magma_2d(nr,R,nz,L,M,order)
% [A,G,Fp,lambda,z,rm,ind] = magma_2d(nr,R,nz,L,M,order)
% magma_2d assemble the coefficient matrix and forcing array for time
% integration including boundary penalty SAT terms.
%
% Input:
% nz    : Number of grid points 
% dz    : Grid spacing
% L      : conduit lengt
% M     : Material properties structure
% order : Interior order of accuracy
%
% Output:
% A      : Spatial discretization (Sparse matrix)
% G      : Time-dependent forcing function
% Fp     : Array for pressure forcing (Sparse)
% lambda : eigenvalues of A
%
% Discretization of the magma 1d equations
% Spatial discretization:
% q_t = A*q,
% where q = [v p n h]

% Grid points and grid spacing
dz = L/nz;
dr = R/nr;

z = dz*[0:nz]'; % Construct z grid

truncate = 1;
[rp,rm,Pp,Pm,Qp,Qm] = sbp_staggered_strong(order,nr, dr, truncate);
[Pz_inv,Dz] = sbp_sparse(order,nz,dz);

rho_inv = spdiags(1./M.rho,0,nz+1,nz+1);

Rp = spdiags(rp,0,nr+1,nr+1);
Rm = spdiags(rm,0,nr,nr);

% Assemble second derivative
D2 = inv(Rm)*inv(Pm)*Qm*Rp*inv(Pp)*Qp;

%%
e0 = spalloc(nz+1,1,1);
en = spalloc(nz+1,1,1);

e0(1)   = 1;
en(end) = 1;

em = ones([nr,1]);
ez = ones([nz+1,1]);
Rm = spdiags(rm,0,nr,nr);

Iz = speye(nz+1);
Im = speye(nr);

% create index for extract solution
if strcmp(M.bottom_bc, 'crack_quasi_static')
    % Indices used to extract physical variables from q
    ind.iv = [1:nr*(nz+1)]';
    ind.ip = [nr*(nz+1)+1: nr*(nz+1)+(nz+1)]';
    ind.in = [nr*(nz+1)+(nz+1)+1: nr*(nz+1)+2*(nz+1)]';
    ind.ih = [nr*(nz+1)+2*(nz+1)+1]';
    ind.ip_c = [nr*(nz+1)+2*(nz+1)+1+1]';
else
    ind.iv = [1:nr*(nz+1)]';
    ind.ip = [nr*(nz+1)+1: nr*(nz+1)+(nz+1)]';
    ind.in = [nr*(nz+1)+(nz+1)+1: nr*(nz+1)+2*(nz+1)]';
    ind.ih = [nr*(nz+1)+2*(nz+1)+1]';
end

W1 = (2/R^2)*em'*Rm*Pm; % Cylindrical width averaged operator that operators on just one row of v, dimension [1, nr]
W2 = kron(Iz,W1); % Operates on the entire v, dimension [nz+1,(nz+1)*nr]

% size(W1),size(W2)

%rho, c, K, a, b are now diagonal matrices while M.rho, M.c, M.K, M.a, M.b remain as vectors.

rho = spdiags(M.rho,0,nz+1,nz+1);
c   = spdiags(M.c  ,0,nz+1,nz+1);
K   = spdiags(M.K  ,0,nz+1,nz+1);
a   = spdiags(M.a  ,0,nz+1,nz+1);
b   = spdiags(M.b  ,0,nz+1,nz+1);

% SAT penalties
SAT = [M.c(1)*Pz_inv(1) M.c(end)*Pz_inv(end)]; % penalty relaxation rate
%SAT = [0 0]; % no SAT boundary conditions

%% construct matrix A by parts
% contribution of A from PDE
A00 = M.mu*kron(rho_inv,D2);
A01 = kron(-inv(rho)*Dz,em)+ kron(-M.g*inv(K),em);
A02 = kron(M.g*a,em);
A03 = sparse(nr*(nz+1),1);

A10 = W2*(kron(-K*Dz,Im) + kron(M.g*rho,Im));
A11 = -K*a*b/M.tau;
A12 = -K*a/M.tau;
A13 = sparse(nz+1,1);

A20 = W2*(kron(-M.g*rho*b,Im));
A21 = -b/M.tau;
A22 = -Iz/M.tau;
A23 = sparse((nz+1),1);

A30 = W1*(kron(en',Im));
A31 = 1/(rho(end)*c(end))*en';
A32 = sparse(1,(nz+1));
A33 = -M.g/c(end);

% 
dim = [nr*(nz+1), nz+1, nz+1,1];

Ai = block_matrix(dim, dim, 0);
Ai = block_matrix_insert(Ai, dim, dim, 1, 1, A00);

% Adding SAT terms according to boundary condition
switch M.bottom_bc
    case 'p=0';
        A01 = A01 + ...
                  kron(-(SAT(1)/(rho(1)*c(1)))*(e0*e0'),em) + ...
                  kron((SAT(end)/(rho(end)*c(end)))*(en*en'),em);
              
        A03 = A03+kron(-SAT(end)*(M.g/c(end))*en,em);
        
        A11 = A11 - SAT(1)*(e0*e0') - SAT(end)*(en*en');
        A13 = A13 + SAT(end)*M.g*rho(end)*en;
    case 'v=0'
        A00 = A00 - SAT(1)*kron(en*en',Im);
        A01 = A01 + kron((SAT(end)/(rho(end)*c(end)))*(en*en'),em); 
        A03 = A03 + kron(-SAT(end)*(M.g/c(end))*en,em);
        A10 = A10 - SAT(1)*M.rho(1)*M.c(1)*W2*kron(en*en',Im); 
        A11 = A11 - SAT(end)*(en*en');
        A13 = A13 + SAT(end)*M.g*rho(end)*en;
    case 'crack_quasi_static'
        A01 = A01 + ...
                  kron(-(SAT(1)/(rho(1)*c(1)))*(e0*e0'),em) + ...
                  kron((SAT(end)/(rho(end)*c(end)))*(en*en'),em);
              
        A03 = A03+kron(-SAT(end)*(M.g/c(end))*en,em);
        
        A11 = A11 - SAT(1)*(e0*e0') - SAT(end)*(en*en');
        A13 = A13 + SAT(end)*M.g*rho(end)*en;
        
        A40 = -M.alpha*W1*kron(e0',Im);
        A41 =  M.alpha*1./(M.rho(1)*M.c(1))*e0';
        A42 = sparse(1,nz+1);
        A43 = sparse(1,1);
        A44 = -M.alpha/(M.rho(1)*M.c(1));
        
        A04 = SAT(1)/(M.rho(1)*M.c(1))*kron(e0,em);
        A14 = SAT(1)*e0;
        A24 = sparse(nz+1,1);
        A34 = sparse(1,1);
end

% Forcing Function, G(t)

G = @(t) M.pT.A*exp(-0.5*((t-M.pT.t)/M.pT.T)^2);

% assemble matrix A and 
switch M.bottom_bc
    case {'p=0','v=0'}
        A0 = horzcat(A00,A01,A02,A03);
        A1 = horzcat(A10,A11,A12,A13);
        A2 = horzcat(A20,A21,A22,A23);
        A3 = horzcat(A30,A31,A32,A33);
        % disp('Ablocks = ')
        % size(A0), size(A1), size(A2), size(A3)
        A = vertcat(A0,A1,A2,A3);
        Ae = A - Ai;
        
        % Array for pressure forcing
        Fp = vertcat(kron(-(SAT(end)/(rho(end)*c(end)))*en,em), ...
                     SAT(end)*en, ...
                     spalloc(nz+1,1,1), ...
                     -1/(rho(end)*c(end)));
    case 'crack_quasi_static'
        A0 = horzcat(A00,A01,A02,A03,A04);
        A1 = horzcat(A10,A11,A12,A13,A14);
        A2 = horzcat(A20,A21,A22,A23,A24);
        A3 = horzcat(A30,A31,A32,A33,A34);
        A4 = horzcat(A40,A41,A42,A43,A44);
        % disp('Ablocks = ')
        % size(A0), size(A1), size(A2), size(A3)
        A = vertcat(A0,A1,A2,A3,A4);

        % Array for pressure forcing
        Fp = vertcat(kron(-(SAT(end)/(rho(end)*c(end)))*en,em), ...
                     SAT(end)*en, ...
                     spalloc(nz+1,1,1), ...
                     -1/(rho(end)*c(end)),...
                     spalloc(1,1,1));
end         

% size(A0), size(A1), size(A2), size(A3)
% disp(['A = ', num2str(size(A))])
lambda=0;
%%
% Check Stability
test = true;
if test
  lambda = eig(full(A));
  disp(['Stability Test: ', num2str(max(real(lambda)))])
  %assert(max(real(lambda))<1e-10,'Unstable');
end
nq = size(A,1);
% Check energy stability
test_energy = true;
b = M.b;
b(b==0) = inf;
if test_energy
    Pz = inv(Pz_inv);
    S_diag = 0.5*[kron(Pz*M.rho,Rm*Pm*em); ...
                  (R^2/2)*Pz*(1./M.K); ...
                  (R^2/2)*Pz*(M.a./b); ...
                  (R^2/2)*M.rho(end)*M.g];
    S_diag(isnan(S_diag)) = 0;
    S = spdiags(S_diag,0,nq,nq);
%     W_diag
    E = S*A + A'*S;  
    lambda_E = eig(full(E));
    disp(['Energy Stability: ', num2str(max(real(lambda_E)))])
end