classdef conduit
  % multiphase magma column, acoustic-gravity wave.
  properties
    M;              % material properties, including conduit geometry and grid configuration.
    op;             % discretization and width averaging operators
    geom;        % grid in z, x direction.
    u;               % unknowns in the conduit, vz, pz, nz, h.
    SAT;           % the SAT term at the location coupled to the crack.
    A;               % the du/dt = A*u+Fp*G.
    Fp;             % operators for external forcing or source.
    Ae;             % the part of A that treated explicitly.
    Ai;              % the part of A that treated implicitly.
    E;               % energy norm
    
    % material properties in M:
    % nr, R, nz, L, order. S (conduit cross-section area), 
    % rho, c, K, a, b.
  end

  methods

    function obj = conduit(M)
      obj.M        = M;
      
      %construct grid and operators
      if ~ isfield(M,'interface_split') || ~M.interface_split
          [obj.geom, obj.op] = grids_conduit(M.nr, M.R, ...
                                             M.nz, M.L, M.order);
      else
          if isfield(M, 'split_index')
              [obj.geom, obj.op] = grids_conduit_split(M.nr, M.R, ...
                                                 M.nz, M.L, M.split_index, M.order);
          else
              error('Please specify split index');
          end
      end
      
      % construct the discretization. A, Ae, Ai, Fp.
      dim = obj.dimensions();
      if  isfield(obj.M, 'with_exsolution') &&  ~obj.M.with_exsolution
           [obj.Ai, obj.Ae, obj.Fp, obj.SAT] = discretize_conduit_no_ex(obj.geom, obj.op, M, dim);
      else
          [obj.Ai, obj.Ae, obj.Fp, obj.SAT] = discretize_conduit(obj.geom, obj.op, M, dim);
      end

      obj.A = obj.Ai + obj.Ae;
      obj.u = zeros(sum(obj.dimensions()),1);
      obj.E = obj.energy_norm();
    end
    
    function obj = init(obj, u0)
        dim = obj.dimensions();
        check_dim = length(u0) == sum(dim);
        if check_dim
            obj.u = u0;
        else
            error('Dimension of u0 is doesn''t match with the grid')
        end
    end

   function obj = update(obj,U)
       obj.u = U;
   end
    
    function dim = dimensions(obj)
        % dimension for vz, pz, nz, h.
      if isfield(obj.M, 'with_exsolution') &&  ~obj.M.with_exsolution
          dim = [obj.geom.nz*obj.geom.nr,...
                    obj.geom.nz,...
                    1];
      else
                    dim = [obj.geom.nz*obj.geom.nr,...
                    obj.geom.nz,...
                    obj.geom.nz,...
                    1];
      end
          
    end

    function u = field(obj,U,num)
        % return the solution, in the order of [1, 2, 3, 4] = [vz, pz, nz, h]
      indices = field_indices(obj.dimensions(),num);
      u = U(indices);
    end
    
    function [Epc, Epg, Ek, Evis]= get_energetics(obj, u)
        % [Epc, Epg, Ek, Evis]= get_energetics(obj, u)
        %
        % Epc : potential energy from fluid compressibility (>0)
        % Epg : potential energy from gravity (>0)
        % Ek   : fluid kinetic energy. (>0)
        % Evis : energy dissipation rate from fluid viscosity. (<0)
        %
        % TODO: exsolution is not implemented.
        %
        
        vz = obj.field(u, 1);
        pz = obj.field(u, 2);
        h = obj.field(u, 3);
        
        Hz = obj.op.Pz;
        
        % rho, K, S must be a vector of dimension [nz, 1]
        rho = obj.M.rho;
        K    = obj.M.K;
        S    = obj.M.S;
        g    = obj.M.g;
        
        nz  = obj.geom.nz;
        nr   = obj.geom.nr;
        
        if ~isfield(obj.M, 'SL')
            SL = S(end);
        else
            SL = obj.M.SL;
        end
        
        % potential energy from gravity
        Epg = 0.5 * rho(end)*g*(h'*h)*SL;
        
        % potential energy from fluid compressibility.
        Epc = 0.5 * pz'*spdiags(S./K, 0, nz, nz)*Hz*pz;
        
        % fluid kinetic energy.
        RHO = spdiags(rho, 0, nz, nz);
        Rm   = obj.op.Rm;
        Pm   = obj.op.Pm;
        Rp    = obj.op.Rp;
        Pp    = obj.op.Pp;
        
        Hv = 2*pi*kron(RHO.*Hz, Rm*Pm);
        Ek   = 0.5*vz'*Hv*vz;
        
        % fluid viscous energy dissipation rate.
        D1r2 = obj.op.D1r2;
        dvz   =  D1r2*vz;
        
        if length(obj.M.mu)==1
            mu = ones(nz, 1)*obj.M.mu;
        end
        
        MU    =  spdiags(mu, 0, nz, nz);
        Hdvz = 2*pi*kron(MU.*Hz, Rp*Pp);
        Evis  = - dvz'*Hdvz*dvz;
    end

    function E = energy_norm(obj)
      dim = obj.dimensions();
      E    = block_matrix(dim,dim,1);
      
      if  isfield(obj.M, 'with_exsolution') &&  ~obj.M.with_exsolution
           % TODO: Add each term in the mechanical energy balance norm for
           % the case without exsolution.
          rho = obj.M.rho;
          K    = obj.M.K;
          S    = obj.M.S;
          nz   = obj.geom.nz;
          E11  = kron(spdiags(rho,0, nz, nz)*obj.op.Pz, 2*pi*obj.op.Rm*obj.op.Pm);
          E22  = spdiags(S./K,0, nz, nz) * obj.op.Pz;
          E33  = S(end)*rho(end)*obj.M.g;
          E    = block_matrix_insert(E, dim, dim, 1, 1, E11);% kinetic energy
          E    = block_matrix_insert(E, dim, dim, 2, 2, E22);% kinetic energy
          E    = block_matrix_insert(E, dim, dim, 3, 3, E33);% kinetic energy
      else
          rho = obj.M.rho;
          K    = obj.M.K;
          a    = obj.M.a;
          b    = obj.M.b;
          b(b==0) = inf;% avoid infinite a./b. b=0 is below the exsolution depth, where n is always equal to zero.
          S    = obj.M.S;
          nz   = obj.geom.nz;
          E11  = kron(spdiags(rho,0, nz, nz)*obj.op.Pz, 2*pi*obj.op.Rm*obj.op.Pm);
          E22  = spdiags(S./K,0, nz, nz) * obj.op.Pz;
          E33  = spdiags(S*a./b,0, nz, nz) * obj.op.Pz;
          E44  = S(end)*rho(end)*obj.M.g;
          E    = block_matrix_insert(E, dim, dim, 1, 1, E11);% kinetic energy
          E    = block_matrix_insert(E, dim, dim, 2, 2, E22);% kinetic energy
          E    = block_matrix_insert(E, dim, dim, 3, 3, E33);% kinetic energy
          E    = block_matrix_insert(E, dim, dim, 4, 4, E44);% kinetic energy.
      end
      
    end

  end

end
