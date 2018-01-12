classdef pipe1d
    % Implement 1d pipe model, acoustic-gravity wave, velocity is cross-sectional averaged.
    % assume Poiseuille flow, du/dt + dp/dz = -alpha*u, alpha = 8*mu/R^2.
    % Don't consider non-equilibrium gas exsolution.
    % The pipe can have multiple interfaces, either 'm' (material
    % properties) or 'c' (coupled to a crack)
   
    properties
    M;              % material properties, including conduit geometry and grid configuration.
    op;             % discretization and width averaging operators
    grd;            % grid in z.
    u;               % unknowns in the conduit, uz, pz, h.
    indu;           % field index
    Ae;             % the part of A that treated explicitly PDE part.
    Ai;              % the part of A that treated implicitly.
    Fp;
    E;               % energy norm
    
    % material properties in M:
    % nr, R, nz, L, order. S (conduit cross-section area), 
    % rho, c, K, a, b.
    end
    
    methods
        function obj = pipe1d(M)
            [check, M] = pipe1d_precheck(M);
            if ~check
                return
            else
                obj.M = M;
            end
            % discretize:
            [obj.grd, obj.op] = grids_pipe1d(M);
            obj = obj.field_index();
            [obj.Ai, obj.Ae, obj.Fp] = obj.build();
            obj = obj.init();
        end
        
    function obj = init(obj)
        dim = obj.dimensions();
        obj.u = zeros(sum(dim), 1);
    end

   function obj = update(obj, u1)
       obj.u = u1;
   end
    
    function dim = dimensions(obj)
        % dimension for vz, pz, h, eta's for each material interface.
        % eta is the moving material interface position.
        if strcmp(obj.M.bc.tp.type, 'h')
          % moving top surface equation for h.
          % Only this one is implemented.
          dim = [obj.grd.nz, obj.grd.nz, 1];
          else
              dim = [obj.grd.nz, obj.grd.nz];
        end
      % loop over all the material interfaces.
      faces = obj.M.material_interfaces;
      nfaces = length(faces);
      for i = 1: nfaces
        % add addition eta's for each material interface.
          dim = [dim, 1];
      end
    end

    function  [vz, pz, h, etas] = fields(obj, u)
%     dimension for vz, pz, h, etas.
        indv  = obj.indu.indv;
        indp  = obj.indu.indp;
        indh  = obj.indu.indh;
        indetas = obj.indu.indetas;
        vz     = u(indv);
        pz     = u(indp);
        h       = u(indh);
        etas  = u(indetas); % moving material interface position.
    end
     
    function  obj = field_index(obj)
%     dimension for vz, pz, h, eta's
        dim   = obj.dimensions();
        obj.indu.indv  = [1:dim(1)]';
        obj.indu.indp  = [(dim(1)+1): sum(dim(1:2))]';
        obj.indu.indh  = sum(dim(1:2)) + 1;
        obj.indu.indetas  = [sum(dim(1:2)) + 2: sum(dim)]';
    end
     
    function [Ai, Ae, Fp] = build(obj)
        
        dim  = obj.dimensions();
        nz    = obj.grd.nz; % grid points for 
        rho    = obj.grd.rho;
        rhoi   = inv(rho);
        K       = obj.grd.K; 
        Ki      = inv(K);
        Dz     = obj.op.Dz;
        g        = obj.M.g;
        S       = obj.grd.S;
        Si      = inv(S);
        alpha = 8*obj.grd.mu*inv(obj.grd.R.^2);
        % build part of the matrix that needs to be integrated implicitly Ai.
        Ai11 = -alpha*rhoi;
        Ai = block_matrix(dim, dim, 0);
        Ai = block_matrix_insert(Ai, dim, dim, 1, 1, Ai11);
        
        % Ae is the PDE part of the operators.
        
        % save each block of Ae into a cell array. Ae_cell.
        % this treatment is to simplify the code when inserting block
        % matrices.
        Ae = block_matrix(dim, dim, 0);
        Ae_cell = cell(length(dim), length(dim));
        
        for i = 1:length(dim)
            for j = 1: length(dim)
                Ae_cell{i, j} = sparse(dim(i), dim(j));
            end
        end
        
        Ae_cell{1,2} = - rhoi*Dz - g*Ki;
        Ae_cell{2,1} = - K*Dz + rho*g;
        
        idx_block = 3;
        
        % build interface sconditions into Ae.
         for i = 1: length(obj.op.interfaces)
                face = obj.op.interfaces{i};
                type = face.type;
                % plus (1) and minus (2) index.
                e1    =  face.ep; %plus
                e2   =  face.em; % minus
                indp = face.indp;
                indm = face.indm;
                rho1 = obj.grd.rho(indp, indp);
                c1    = obj.grd.c(indp, indp);
                rho2 = obj.grd.rho(indm, indm);
                c2    = obj.grd.c(indm, indm);
                Z1    = rho1*c1;% plus impedance
                Z2    = rho2*c2;% minus impedance
                
                SAT1= face.SATp;
                SAT2= face.SATm;
                switch type 
                    case 'm'
                        % add terms related to material interfaces.
                        % block index of this material interface.
                        % vz, pz, h, eta1, eta2, eta3 ...
                        % idx_block start from 4 and increment by one after
                        % each material interface is processed.
                        
                        idx_block = idx_block + 1; 
                        % vz <-- vz
                        Ae_cell{1,1} = Ae_cell{1,1} - SAT1*Z2/(Z1+Z2)*e1*(e1' - e2')  ...
                                    + SAT2* Z1/(Z1+Z2)*e2*(e1' - e2');
                        % vz <-- pz
                        Ae_cell{1,2} = Ae_cell{1,2}+ 1/(Z1+Z2)*(SAT1*e1*(e2' - e1') + SAT2*e2*(e2' - e1'));
                        % pz <-- pz
                        Ae_cell{2,2} = Ae_cell{2,2}  - SAT1*Z1/(Z1+Z2)*(e1*(e1'-e2')) ...
                                    + SAT2*Z2/(Z1+Z2)*(e2*(e1'-e2'));
                        % pz <-- vz
                        Ae_cell{2,1} = Ae_cell{2,1} - Z1*Z2/(Z1+Z2)*(SAT1*e1*(e1' - e2') + SAT2*e2*(e1' - e2'));
                        
                        % vz  <-- eta
                        Ae_cell{1, idx_block} =  Ae_cell{1, idx_block} + (rho1-rho2)*g/(Z1+Z2)*(SAT1*e1 + SAT2*e2);
                        % pz  <-- eta
                        Ae_cell{2, idx_block} =  Ae_cell{2, idx_block} + (rho1-rho2)*g/(Z1+Z2)*(SAT1*Z1*e1 - SAT2*Z2*e2);
                        
                        % eta <-- eta
                        Ae_cell{idx_block, idx_block} =  Ae_cell{idx_block, idx_block} + (rho1-rho2)*g/(Z1+Z2);
                        % kinematic condition for eta.
                        % eta <-- vz
                        Ae_cell{idx_block, 1} =  Ae_cell{idx_block, 1} + 1/(Z1+Z2)*(Z1*e1' + Z2*e2');
                        % eta <-- pz
                        Ae_cell{idx_block, 2} =  Ae_cell{idx_block, 2} + 1/(Z1+Z2)*(-e1' + e2');
                    case 'c'
                        % add terms contributed by from pipe unknowns.
                        % coupling with crack doesn't introduce additional
                        % unknown to the pipe.
                        Ae_cell{1,2} = Ae_cell{1,2} + 1/(Z1*Z2)*(Z1*SAT2*(e2*e2') - Z2*SAT1*(e1*e1'));
                        Ae_cell{2,2} = Ae_cell{2,2} - SAT1*(e1*e1') - SAT2*(e2*e2');
                    otherwise
                end
         end
         
         % Note that, if the interface is coupled to the crack, additional
         % terms need to added, but this requires knowing the information
         % of the crack, which will be added when building a coupled model.

         % check if the bottom is p=0 or v=0 boundary condition. if so, add
         % these terms into the Ae.
         bt = obj.op.bc.bt;
         type = bt.type;
         SAT = bt.SAT;
         e = bt.e;
         Z    = obj.grd.rho(1, 1) * obj.grd.c(1, 1);
         switch type
             case 'v'
                 % v0 = 0
                 Ae_cell{1,1} = Ae_cell{1,1} - SAT*(e*e');
                 Ae_cell{2,1} = Ae_cell{2,1} - SAT*Z*(e*e');
             case 'p'
                 % p0 = 0
                 Ae_cell{2,2} = Ae_cell{2,2} - SAT*(e*e');
                 Ae_cell{1,2} = Ae_cell{1,2} - SAT/Z*(e*e');
             case 'c'
                 % for crack bottom, build in the contribution from pipe
                 % pressures.
                 Ae_cell{2,2} = Ae_cell{2,2} - SAT*(e*e');
                 Ae_cell{1,2} = Ae_cell{1,2} - SAT/Z*(e*e');
         end
         
         % compute the source term vector from surface excitation.
         tp = obj.op.bc.tp;
         SAT = tp.SAT;
         e      = tp.e;
         Z      = obj.grd.rho(end, end) * obj.grd.c(end, end);
         rho_g = obj.grd.rho(end, end)*obj.M.g;
         
         % top moving surface.
         Fp = block_matrix(dim, 1, 0);
         Ae_cell{1,2} =  Ae_cell{1,2} + SAT/Z*(e*e');
         Ae_cell{1,3} =  Ae_cell{1,3} - SAT/Z*rho_g*e;
         Fp1    = -SAT/Z*e;
         Ae_cell{2,2} = Ae_cell{2,2} - SAT*(e*e');
         Ae_cell{2,3} = Ae_cell{2,3} + SAT*rho_g*e;
         Fp2    = SAT*e;
         Ae_cell{3,1}  = Ae_cell{3,1} + e';
         Ae_cell{3,2}  = Ae_cell{3,2} + e'/Z;
         Ae_cell{3,3}  = Ae_cell{3,3} - rho_g/Z;
         Fp3    = 1/Z;
         
         % insert Ae_cell into Ae, not good practice but okay for 1D problem.
         for i = 1: length(dim)
             for j=1: length(dim)
                 Ae = block_matrix_insert(Ae, dim, dim, i, j, Ae_cell{i, j});
             end
         end
         % source term only for top bc.
         Fp = block_matrix_insert(Fp, dim, 1, 1, 1, Fp1);
         Fp = block_matrix_insert(Fp, dim, 1, 2, 1, Fp2);
         Fp = block_matrix_insert(Fp, dim, 1, 3, 1, Fp3);
    end
    
        function [cmax, hmin] = getCFL(obj)
            % get min(c) and max(dz) for CFL condition.
            cmax = max(obj.M.c);
            hmin  = min([obj.grd.dz]);
        end
    
    function [du]= fun_integrate(obj, u, t)
        % a function for integrating the explicit
        du    =  obj.Ae*u + obj.Fp*obj.op.bc.tp.f(t);
    end
    
    function E = energy_norm(obj)
        % calculate the energy norm of the pipe
        % vz, pz, h, etas.
          dim = obj.dimensions();
          E    = block_matrix(dim,dim,1);
          rho = obj.grd.rho;
          K    = obj.grd.K;
          S    = obj.grd.S;
          nz   = obj.grd.nz;
          E11  = rho*obj.op.H;
          E22  = S*inv(K)*obj.op.H;
          E33  = S(end, end)*rho(end, end)*obj.M.g;
          
          E    = block_matrix_insert(E, dim, dim, 1, 1, E11);% kinetic energy
          E    = block_matrix_insert(E, dim, dim, 2, 2, E22);% compressibility
          E    = block_matrix_insert(E, dim, dim, 3, 3, E33);% gravity of top surface.
          % gravity of material interfaces, etas.
          
          %get all the material interfaces.
          nint  = length(obj.M.interfaces);
          
          inds = [];
          for i = 1: nint
              if obj.M.interfaces{i}.type=='m'
                  inds = [inds; i];
              end
          end
          
          % loop over all the material interfaces.
          nfaces = length(inds);
          j = 3; % vz,p,h
          for i = 1: nfaces
            % add addition eta's for each material interface.
              j = j + 1;
              ind_m = obj.op.interfaces{inds(i)}.indm;
              ind_p  = obj.op.interfaces{inds(i)}.indp;
              rho_m = obj.grd.rho(ind_m, ind_m);
              rho_p  = obj.grd.rho(ind_p, ind_p);
              drho    = rho_m - rho_p;
              S_m = obj.grd.S(ind_m, ind_m);
              S_p  = obj.grd.S(ind_p, ind_p);
              S_j   = (S_m + S_p)/2; 
              % get the jump in density across the interfaces
              E    = block_matrix_insert(E, dim, dim, j,  j, drho*S_j*obj.M.g);% gravity of top surface.
          end

    end
    end
    
end

