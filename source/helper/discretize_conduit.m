function [Ai, Ae, Fp, SAT] = discretize_conduit(geom, op, M, dim)
% discretize conduit using grid (geom), operators (op) and material properties (M).
%
% Ai and Ae contain all the contribution from [vz, pz, nz, h] of the
% conduit but do not include terms from [px, vx] from the fracture.

nz = geom.nz;
nr = geom.nr;

rho = spdiags(M.rho,0,nz,nz);
c   = spdiags(M.c  ,0,nz,nz);
K   = spdiags(M.K  ,0,nz,nz);
a   = spdiags(M.a  ,0,nz,nz);
b   = spdiags(M.b  ,0,nz,nz);

% SAT penalties
SAT = [M.c(1)*op.Pz_inv(1, 1) M.c(end, end)*op.Pz_inv(end)]; % penalty relaxation rate
%SAT = [0 0]; % no SAT boundary conditions

rho_inv = spdiags(1./M.rho,0,nz,nz);
%% construct matrix A by parts
Ai = block_matrix(dim,dim);
er = op.er;
en = op.en;
Ir  = op.Ir;

% contribution of A from PDE
A11 = M.mu*kron(rho_inv,op.D2);
Ai = block_matrix_insert(Ai,dim,dim,1,1, A11); % this part should be treat implicitly;

A12 = kron(-inv(rho)*op.Dz, er)+ kron(-M.g*inv(K),er);
A13 = kron(M.g*a, er);
A14 = sparse(nr*(nz),1);

A21 = op.W2*(kron(-K*op.Dz,Ir) + kron(M.g*rho,Ir));
A22 = -K*a*b/M.tau;
A23 = -K*a/M.tau;
A24 = sparse(nz,1);

A31 = op.W2*(kron(-M.g*rho*b,Ir));
A32 = -b/M.tau;
A33 = -op.Iz/M.tau;
A34 = sparse((nz),1);

A41 = op.W1*(kron(en',Ir)); %M.S/M.SL
A42 = 1/(rho(end)*c(end))*en'; %M.S/M.SL
A43 = sparse(1,(nz));
A44 = -M.g/c(end);%M.S/M.SL*

% contribution from the SAT term of boundaries:
A12 = A12 + ...
          kron(-(SAT(1)/(rho(1)*c(1)))*(op.e0*op.e0'),er) + ...
          kron((SAT(end)/(rho(end)*c(end)))*(en*en'),er);
      
A14 = A14+kron(-SAT(end)*(M.g/c(end))*en,er);

A22 = A22 - SAT(1)*(op.e0*op.e0') - SAT(end)*(en*en');
A24 = A24 + SAT(end)*M.g*rho(end)*en;


if M.interface_split
    % contribution from the SAT terms of the interface.
    % Implementation follow Karlstrom and Dunham (2016).
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

