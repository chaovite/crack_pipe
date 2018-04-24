classdef coupledModel
  % A coupled conduit-crack model
  properties
    conduit;              % material properties, including conduit geometry
    frac;        % configurations for grid.
    A;               % the du/dt = A*u+Fp*G.
    Fp;             % operators for external source from both the conduit and fracture.
    Ae;             % the part of A that treated explicitly.
    Ai;              % the part of A that treated implicitly.
    E;               % energy norm
    u;               % unknowns vz, pz, nz, h, px, vx.
  end

  methods
    function obj = coupledModel(conduit, frac)
        %% construct a coupled model with a conduit and frac
      obj.conduit  = conduit;
      obj.frac        = frac;
      [obj.A, obj.Ai, obj.Ae, obj.Fp] = obj.build();% build the operators.
      obj.u = zeros(sum(obj.dimensions()),1);
      obj.E = obj.energy_norm();
      
      % empty the memory of A, Ai, Ae matrices from conduit and crack.
      obj.conduit.A = [];
      obj.conduit.Ai = [];
      obj.conduit.Ae = [];
      obj.frac.A = [];
      obj.frac.Ai = [];
      obj.frac.Ae = [];
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
    
    function dim = dimensions(obj)
        %[conduit, frac].
        % index
        dim_c = obj.conduit.dimensions();
        dim_f =  obj.frac.dimensions();
        dim = [sum(dim_c),sum(dim_f)];
    end

    function u = field(obj,U,num)
        % if num = [1,1], return the vz in the conduit.
        % if num = [1], return all the unknowns in the conduit.
        if length(num)==1
          indices = field_indices(obj.dimensions(),num);
          u = U(indices);
        else if length(num) ==2
                indices1 = field_indices(obj.dimensions(),num(1));
                switch num(1)
                    case 1
                        u = obj.conduit.field(U(indices1), num(2));
                    case 2 
                        u = obj.frac.field(U(indices1), num(2));
                end
            end
        end
    end
    
        function varargout = fields(obj, U)
            vars_conduit = length(obj.conduit.dimensions());
            vars_crack    = length(obj.frac.dimensions());
            nout = vars_conduit + vars_crack;
            nz = obj.conduit.geom.nz;
            nr = obj.conduit.geom.nr;
            varargout = cell(nout, 1);
            
            for k = 1: vars_conduit
                var = obj.field(U, [1, k]);
                if k ==1
                    var = reshape(var, nr, nz);
                end
                varargout{k} = var; 
            end
       
           % solution in the fluid.
           frac_model = '3D';
           if ~isfield(obj.frac.M,'nz')
               frac_model = '2D';
           end
           switch frac_model
               case '3D'
                   varargout{vars_conduit + 1} = obj.frac.geom.p.grd(obj.field(U,[2,1]));
                   varargout{vars_conduit + 2}  = obj.frac.geom.vx.grd(obj.field(U,[2, 2]));
                   varargout{vars_conduit + 3} = obj.frac.geom.vy.grd(obj.field(U,[2, 3]));
               case '2D'
                   varargout{vars_conduit + 1}  = obj.field(U,[2 1]);
                   varargout{vars_conduit + 2} = obj.frac.geom.mm.grd(obj.field(U,[2 2]));
           end
           
        end
    
    function [A, Ai, Ae, Fp] = build(obj)
        % this method build the A, Ai, Ae, Fp from both conduit and
        % fluid-filled fractures and aggregate into large A, Ai, Ae, Fp.
          dim = obj.dimensions();
          Ae = block_matrix(dim, dim);
          Ai = block_matrix(dim, dim);
          % construct the block diagonal terms
          Ae = block_matrix_insert(Ae,dim,dim,1,1, obj.conduit.Ae);
          Ae = block_matrix_insert(Ae,dim,dim,2,2, obj.frac.Ae);
          Ai = block_matrix_insert(Ai,dim,dim,1,1, obj.conduit.Ai);
          Ai = block_matrix_insert(Ai,dim,dim,2,2, obj.frac.Ai);
         
          % construct the forcing terms:
          nc = 2; % number of components. one for conduit and another one for frac.
          Fp = block_matrix(dim, ones(1, nc));
          Fp = block_matrix_insert(Fp,dim, [1 1], 1,1, obj.conduit.Fp);
          Fp = block_matrix_insert(Fp,dim, [1 1], 2,2, obj.frac.Fp);
          
          % adding the coupling terms, all the coupling terms are explicit
         % in this case.
         Ae_c =  block_matrix(dim, dim, 0);
         
         dim_c = obj.conduit.dimensions();
         dim_f  = obj.frac.dimensions(); 
         
         A12 = block_matrix(dim_c, dim_f, 0);
         A21 = block_matrix(dim_f, dim_c, 0);
         A22 = block_matrix(dim_f, dim_f, 0);
         
         % fracture to conduit.
         SAT = obj.conduit.SAT;
         M_c   = obj.conduit.M;
         op_c = obj.conduit.op;
         A12_11 = SAT(1)/(M_c.rho(1)*M_c.c(1)) * kron(op_c.e0, op_c.er) * obj.frac.dc'; %[vz, px]
         A12_21 = SAT(1)*op_c.e0*obj.frac.dc'; %[pz, px]
         A12   = block_matrix_insert(A12, dim_c, dim_f, 1, 1, A12_11);
         A12   = block_matrix_insert(A12, dim_c, dim_f, 2, 1, A12_21);
         
         % conduit to fracture.
         if ~isfield(obj.frac.M,'nz')
             alpha = obj.conduit.M.S(1)/obj.frac.M.S;
         else
             alpha = obj.conduit.M.S(1)/obj.frac.M.w0;
         end
         Kt = obj.frac.M.K_t;
          if ~isfield(obj.frac.M,'nz')
             H = obj.frac.op.p.Px;
          else
             H = obj.frac.op.p.Pxy2;
          end
         d = obj.frac.dc;
         
         A21_11 = - alpha* Kt*inv(H)*d*op_c.e0'*op_c.W2;
         A21_12 = - alpha* Kt*inv(H)*d*(-op_c.e0'/(M_c.rho(1)*M_c.c(1)));
         A21  = block_matrix_insert(A21, dim_f, dim_c, 1 , 1, A21_11);
         A21  = block_matrix_insert(A21, dim_f, dim_c, 1 , 2, A21_12);
         
         Ae_c = block_matrix_insert(Ae_c,dim, dim, 1, 2, A12);
         Ae_c = block_matrix_insert(Ae_c,dim, dim, 2, 1, A21);
         
         % fracture to fracture.
         A22_11 = - alpha* Kt*inv(H)*d*(d'/(M_c.rho(1)*M_c.c(1)));
         A22 = block_matrix_insert(A22, dim_f, dim_f, 1,1, A22_11);
         Ae_c = block_matrix_insert(Ae_c,dim, dim, 2, 2, A22);
         
         Ae = Ae + Ae_c;
         A   = Ae + Ai;
    end
    
    function Hp = get_Hp(obj, xq, yq)
        % get the observation matrix.
            Nq  = length(xq);
            dim = obj.dimensions();
            dim_c = obj.conduit.dimensions();
            dim_f  = obj.frac.dimensions(); 
            N    = sum(dim);
            Hp = zeros(3*Nq, N);
            Hp_i = obj.frac.disloc3d.get_Hp(xq, yq);
            indp  = sum(dim_c) + [1: dim_f(1)]';
            Hp(:, indp) = Hp_i;
            Hp = sparse(Hp);
    end
    
    function Hp_tilt = get_Hp_tilt(obj, xq, yq, omega)
        % get the observation matrix for surface tilt at angular frequency
        % omega.
            obj.conduit.M.g     = 9.8; % gravitational acceleration.
            Nq  = length(xq);
            dim = obj.dimensions();
            dim_c = obj.conduit.dimensions();
            dim_f  = obj.frac.dimensions(); 
            N    = sum(dim);
            Hp_tilt = zeros(3*Nq, N);
            Hp_i = obj.frac.disloc3d.get_Hp_tilt(xq, yq);
            indp  = sum(dim_c) + [1: dim_f(1)]';
            Hp_tilt(:, indp) = Hp_i*g/(1i*omega)^2;
            Hp_tilt = sparse(Hp_tilt);
    end
    
    function obj = update(obj,U)
       obj.u = U;
       obj.conduit = obj.conduit.update(obj.field(U,1));
       obj.frac = obj.frac.update(obj.field(U,2));
    end
    
    function [Epf_p, Epg_p, Ek_p, Evis_p, Epf_c, Epw_c, Ek_c, Evis_c] = get_energetics(obj, u)
        % do the energetics analysis given a solution vector u.
        upipe = obj.field(u, 1);
        ufrac = obj.field(u,  2);
        [Epf_c, Epw_c, Ek_c, Evis_c]= obj.frac.get_energetics(ufrac);
        [Epf_p, Epg_p, Ek_p, Evis_p]= obj.conduit.get_energetics(upipe);
    end
    
    function [D, V, D_full,V_full] = eigs(obj)
        %   returns eigenvalues in the order of descending period. 
        %       D: a vector contain the eigen value
        %       V: a matrix, each column contains the corresponding eigenvector
        %           given each eigenvalue in D
        % seek all the eigen values of A;
        [V_full, D_full] = eig(full(obj.A));

        v=V_full;
        lambda=D_full;

        lambda=diag(lambda);

        %select modes with non-zero imaginary part.

        % sort the eigen-values in terms of the imaginary part in descend order
        [~, I] = sort(imag(lambda),1);
        lambda=lambda(I);

        v = v(:,I);% sort eigen vectors;

        % select eigen values with non-zero imaginary parts
        eps  = 1e-6;
        D     = lambda(imag(lambda)>eps);
        V     = v(:,imag(lambda)>eps); 
    end
    
    function E = energy_norm(obj)

      dim = obj.dimensions();
      E    = block_matrix(dim,dim,1);
      % TODO: Add each term in the mechanical energy balance for the coupled conduit and crack system.
      E    = block_matrix_insert(E, dim,dim, 1,1, obj.conduit.E);
      E    = block_matrix_insert(E, dim,dim, 2,2, obj.frac.E);
    end


  end

end
