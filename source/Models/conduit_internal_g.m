classdef conduit_internal_g
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

    function obj = conduit_internal_g(M)
      obj.M        = M;
      
      %construct grid and operators
      if ~ isfield(M,'interface_split') || ~M.interface_split
          [obj.geom, obj.op] = grids_conduit(M.nr, M.R, ...
                                             M.nz, M.L, M.order, M.order_r);
      else
          if isfield(M, 'split_index')
              [obj.geom, obj.op] = grids_conduit_split(M.nr, M.R, ...
                                                 M.nz, M.L, M.split_index, M.order, M.order_r);
          else
              error('Please specify split index');
          end
      end
      
      % construct the discretization. A, Ae, Ai, Fp.
      dim = obj.dimensions();
      
     [obj.Ai, obj.Ae, obj.Fp, obj.SAT]  = discretize_conduit_internal_g(obj.geom, obj.op, M, dim);

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
   
   function [cmax, hmin] = getCFL(obj)
        % get min(c) and max(dz) for CFL condition.
        cmax = max(obj.M.c);
        hmin  = min([obj.geom.dz]);
    end
    
    function dim = dimensions(obj)
        % dimension for vz, pz, h, hL.
        dim = [obj.geom.nz*obj.geom.nr,...
        obj.geom.nz,...
        obj.geom.nz,...
        1];
          
    end

    function u = field(obj,U,num)
        % return the solution, in the order of [1, 2, 3, 4] = [vz, pz, h, hL]
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
        hL = obj.field(u, 4); 
        
        Hz = obj.op.Pz;
        
        % rho, K, S must be a vector of dimension [nz, 1]
        rho = obj.M.rho;
        K    = obj.M.K;
        S    = obj.M.S;
        Mg = obj.M.Mg;
        g    = obj.M.g;
        
        nz  = obj.geom.nz;
        nr   = obj.geom.nr;
        epsilon = obj.M.epsilon; % A/AL;
        
        % potential energy from gravity (two contribution, internal and surface)
        Epg = 0.5 * epsilon*rho(end)*g*(hL'*hL)*S(end, end) + ...
                  0.5 * h'*spdiags(g*rho.*Mg.*S, 0, nz, nz)*Hz*h; 
        
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
        
        mu = obj.M.mu;
        
        if length(mu)==1
            mu = ones(nz, 1)*mu;
        else
            
        end
        
        MU    =  spdiags(mu, 0, nz, nz);
        Hdvz = 2*pi*kron(MU.*Hz, Rp*Pp);
        Evis  = - dvz'*Hdvz*dvz;
    end

    function E = energy_norm(obj)
      dim = obj.dimensions();
      E    = block_matrix(dim,dim,1);
      rho = obj.M.rho;
      g    = obj.M.g;
      K    = obj.M.K;
      Mg = obj.M.Mg;
      S    = obj.M.S;
      nz   = obj.geom.nz;
      
      E11  = kron(spdiags(rho,0, nz, nz)*obj.op.Pz, 2*pi*obj.op.Rm*obj.op.Pm);
      E22  = spdiags(S./K,0, nz, nz) * obj.op.Pz;
      E33  = spdiags(rho.*Mg.*g.*S,0, nz, nz) * obj.op.Pz;
      E44  = S(end)*rho(end)*g;
      
      E    = block_matrix_insert(E, dim, dim, 1, 1, E11);% kinetic energy
      E    = block_matrix_insert(E, dim, dim, 2, 2, E22);% potential energy from pressure
      E    = block_matrix_insert(E, dim, dim, 3, 3, E33);% potential energy from internal gravity.
      E    = block_matrix_insert(E, dim, dim, 4, 4, E44);% potential energy from surface gravity.
      
    end

  end

end
