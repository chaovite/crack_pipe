function [Ai, Ae, Fp, SAT] = discretize_conduit_internal_g(geom, op, M, dim)
% discretize conduit using grid (geom), operators (op) and material properties (M).
%
% Ai and Ae contain all the contribution from [vz, pz, h, hL] of the conduit
%
% The formulation doesn't assume identity used by Karlstrom and Dunham
% 2016, thus preserving the internal gravity effect.
%
% This implementation does not consider BGR, no non-equlibrium exsolution is considered.
%

nz = geom.nz;
nr = geom.nr;
Iz = op.Iz;

% all the material properties in M should be a vector with dimension  [nz, 1]
propnames = {'Mg', 'rho','c','mu','K'};
for i = 1:length(propnames)
    name = propnames{i};
    prop = M.(name);
    if length(prop) ==  1
        M.(name)   = ones(nz, 1)*prop;
    end
end

Mg = spdiags(M.Mg, 0, nz, nz); % the coefficient that determines the gravity effect.
rho = spdiags(M.rho,0,nz,nz);
c   = spdiags(M.c  ,0,nz,nz);
K   = spdiags(M.K  ,0,nz,nz);
mu = spdiags(M.mu, 0, nz, nz);

% acoustic impedance.
Z_n = rho(end,end)*c(end,end); 
Z_0 = rho(1,1)*c(1, 1); 

% SAT penalties
SAT = [M.c(1)*op.Pz_inv(1, 1) M.c(end, end)*op.Pz_inv(end)]; % penalty relaxation rate

rho_inv = spdiags(1./M.rho,0,nz,nz);
%% construct matrix A by parts
Ai = block_matrix(dim,dim);
er = op.er;
en = op.en;
Ir  = op.Ir;

% contribution of A from PDE
A11 = kron(mu*rho_inv,op.D2);
Ai = block_matrix_insert(Ai,dim,dim,1,1, A11); % this part should be treat implicitly;

A12 = kron(-inv(rho)*op.Dz, er)+ kron(-M.g*inv(K),er);
A13 = -kron(M.g*Mg, er);
A14 = sparse(nr*(nz),1);

A21 = op.W2*(kron(-K*op.Dz,Ir) + kron(M.g*rho,Ir));

A22 = sparse(nz, nz);
A23 = sparse(nz, nz);
A24 = sparse(nz,1);

epsilon = M.epsilon; % area ratio between conduit and lake. A/A_L

A31 = Iz*op.W2;
A32 = sparse(nz,nz);
A33 = sparse(nz,nz);
A34 = sparse(nz,1);

A41 = op.W1*(kron(en',Ir));
A42 = 1/Z_n*en'; %M.S/M.SL
A43 = sparse(1,nz);
A44 = -epsilon*M.g/c(end);%M.S/M.SL*

% contribution from the SAT term of boundaries:
A12 = A12 + ...
          kron(-(SAT(1)/(rho(1)*c(1)))*(op.e0*op.e0'),er) + ...
          kron((SAT(end)/(rho(end)*c(end)))*(en*en'),er);
      
A14 = A14 + epsilon*kron(-SAT(end)*(M.g/c(end))*en,er);

A22 = A22 - SAT(1)*(op.e0*op.e0') - SAT(end)*(en*en');
A24 = A24 +epsilon * SAT(end)*M.g*rho(end)*en;


if M.interface_split
    % contribution from the SAT terms of the interface.
    % Implementation follow Karlstrom and Dunham (2016).
    %
    
    eIm = op.eIm; % minus side.
    eIp = op.eIp; % plus side.
    indm = find(eIm);
    indp = find(eIp);
    
    Zp = M.rho(indp)*M.c(indp);
    Zm = M.rho(indm)*M.c(indm);
    SATm =  M.c(indm)*op.Pz_inv(indm, indm);
    SATp =  M.c(indp)*op.Pz_inv(indp, indp);
    
    A11 = A11 - SATm*kron(eIm, kron(eIm, Ir)') -SATp* kron(eIp, kron(eIp, Ir)') + ...
                        kron(SATm*eIm + SATp*eIp, ...
                                (Zp*kron(eIp, Ir)' + Zm*kron(eIm, Ir)')/(Zp + Zm));
                            
    A12 = A12 + kron(SATm*eIm + SATp*eIp, ...
                                (er*(eIm-eIp)')/(Zp + Zm));
    
    A21 = A21 + (SATm*eIm + SATp*eIp)*Zp*Zm*(eIm-eIp)'*op.W2/(Zp+Zm);
    A22 = A22 - SATm*(eIm*eIm') - SATp*(eIp*eIp')  + ...
                        (SATm*eIm + SATp*eIp)*(Zp*eIm' + Zm*eIp')/(Zp+Zm);
end

A1 = horzcat(A11,A12,A13,A14);
A2 = horzcat(A21,A22,A23,A24);
A3 = horzcat(A31,A32,A33,A34);
A4 = horzcat(A41,A42,A43,A44);
A = vertcat(A1,A2,A3,A4);

Ae = A - Ai;

% Array for pressure forcing
Fp = vertcat(kron(-(SAT(end)/(rho(end)*c(end)))*en,er), ...
             SAT(end)*en, ...
             spalloc(nz,1,1), ...
             -1/(rho(end)*c(end)));
end

