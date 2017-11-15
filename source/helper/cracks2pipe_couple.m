function [A, Ap] = cracks2pipe_couple(pipe, cracks, dim)
% A = cracks2pipe_couple(pipe, cracks)
% calculate the contribution from cracks to the pipe.
%
% cracks, a structure that contains all crack objects.
% pipe: pipe object
% dim: dimension of the entire model including pipe and crack.
%
% A structure, each field is a name of crack interface.
% A is crack's contribution to the pipe.
%
% Ap structure, each field is a name of crack interface.
% Ap the matrix used to update pressure of each crack before multiplied by
% the dense matrix Kt. 
%
% Note that, a crack in the middle and a crack at the bottom has different
% contribution to the conduit.
%

% crack names.
names = fields(cracks);
nc = length(names);
dim_p = pipe.dimensions();

% calculate A
for i = 1: nc
    name = names{i};
    dim_f  = cracks.(name).dimensions();
    A12 = block_matrix(dim_p, dim_f, 0);
    dc     = cracks.(name).dc;
    switch name
        case 'bt'
            % a crack intersect pipe at the bottom.
             % crack to pipe
             SAT = pipe.op.bc.bt.SAT;
             c      = pipe.grd.c(1,1);
             rho   = pipe.grd.rho(1,1);
             e0    = pipe.op.bc.bt.e;
             
             A12_11 = SAT/(rho * c) * e0 * dc'; %[vz, px]
             A12_21 = SAT*e0 * dc'; %[pz, px]
        otherwise
            strs   =strsplit(name,'_'); % interface_1 is 1.
            nface = str2double(strs{2});
            % a crack intersect pipe in the middle.
            % get the index of the interface
            face = pipe.op.interfaces{nface};
            indm = face.indm;
            indp = face.indp;
            % plus:1, minus: 2
            SAT1 = face.SATp;
            SAT2 = face.SATm;
            Z1      = pipe.grd.c(indp,indp) * pipe.grd.rho(indp,indp);
            Z2      = pipe.grd.c(indm,indm) * pipe.grd.rho(indm,indm);
            e1       = face.ep;
            e2       = face.em;
            
            A12_21 = (SAT1*e1 + SAT2*e2) * dc'; %[pz, px]
            A12_11 = 1/(Z1*Z2)*(Z2*SAT1*e1 - Z1*SAT2*e2 )* dc'; %[vz, px]
    end
     A12   = block_matrix_insert(A12, dim_p, dim_f, 1, 1, A12_11);
     A12   = block_matrix_insert(A12, dim_p, dim_f, 2, 1, A12_21);
     
     % put A12 into A.
     A.(name) = A12;
end

% Calculate Ap, used to update crack pressure.

% calculate Ap, which is used to calculate pressure, before * Kt
% this is done to avoid extract individual unknowns during time
% integration.
% p = Ap * u.
% 

for i = 1: nc
    name = names{i};
    dim_f = cracks.(name).dimensions();
    
    idx_crack   = i + 1; % index of block in entire unknown.
    npxy   = dim_f(1); % pxy unknowns 
    Api      = block_matrix(npxy, dim, 0);
    Hinv    = inv(cracks.(name).op.p.Pxy2);
    
    Api_pipe        = block_matrix(npxy, dim_p, 0);
    Api_crack      = block_matrix(npxy, dim_f, 0);

    % we need to insert [1, 1], [1, idx_f]
    dc     = cracks.(name).dc;
    w0    = cracks.(name).M.w0;
    Dx2  =cracks.(name).op.vx.Dx2;
    Dy2  =cracks.(name).op.vy.Dy2;
    
    switch name
        case 'bt'
            % a crack intersect pipe at the bottom.
             % crack to pipe
             c      = pipe.grd.c(1,1);
             rho   = pipe.grd.rho(1,1);
             S      = pipe.grd.S(1, 1);
             Z      = rho*c; 
             e0    = pipe.op.bc.bt.e;
             % vz -> pxy
             Api_pipe11 =    -S/w0*Hinv*dc*e0';
             % pz -> pxy
             Api_pipe12 =     S/w0/Z*Hinv*dc*e0';
             % pxy -> pxy, 
             % This term might be unstable and should be moved to the
             % implicit side.
%              Api_crack11 =  -S/w0/Z*Hinv*(dc*dc');
             % vx -> pxy
             Api_crack12 =  -Dx2;
             % vy -> pxy
             Api_crack13 =  -Dy2 ;
        otherwise
            strs   =strsplit(name,'_'); % interface_1 is 1.
            nface = str2double(strs{2});
            % a crack intersect pipe in the middle.
            % get the index of the interface
            % 
            face = pipe.op.interfaces{nface};
            indm = face.indm;
            indp = face.indp;
            % plus:1, minus: 2
            Z1      = pipe.grd.c(indp,indp) * pipe.grd.rho(indp,indp);
            Z2      = pipe.grd.c(indm,indm) * pipe.grd.rho(indm,indm);
            S        = pipe.grd.S(indp, indp); % Note that this assumes the area doesn't jump. 
            e1      = face.ep;
            e2      = face.em;
            
             % vz -> pxy
             Api_pipe11 =    -S/w0*Hinv*dc*(e1' - e2');
             % pz -> pxy
             Api_pipe12 =     S/w0*Hinv*dc*(e1'/Z1+e2'/Z2);
             % pxy -> pxy
%              Api_crack11 =  -S/w0*(1/Z1+1/Z2)*Hinv*(dc*dc');
             % vx -> pxy
             Api_crack12 =  -Dx2;
             % vy -> pxy
             Api_crack13 =  -Dy2 ;
    end
    
    % insert Api_pipe11 into Api_pipe, vz->pxy
    Api_pipe = block_matrix_insert(Api_pipe, npxy, dim_p, 1, 1, Api_pipe11);
    % insert Api_pipe12 into Api_pipe, pz ->pxy.
    Api_pipe = block_matrix_insert(Api_pipe, npxy, dim_p, 1, 2, Api_pipe12);
    % insert Api_crack11 into Api_crack, pxy->pxy
%     Api_crack = block_matrix_insert(Api_crack, npxy, dim_f, 1, 1, Api_crack11);
    % insert Api_crack12 into Api_crack, vx->pxy
    Api_crack = block_matrix_insert(Api_crack, npxy, dim_f, 1, 2, Api_crack12);
    % insert Api_crack13 into Api_crack, vy->pxy
    Api_crack = block_matrix_insert(Api_crack, npxy, dim_f, 1, 3, Api_crack13);
    % insert Api_crack and Api_pipe into Api.
    Api =  block_matrix_insert(Api, npxy, dim, 1, 1, Api_pipe);
    Api =  block_matrix_insert(Api, npxy, dim, 1, idx_crack, Api_crack);
    % put Api into A.
    Ap.(name) = Api;
end

end

