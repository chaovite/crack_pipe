classdef coupledModel_NoMatrix
  % A coupled conduit-crack model without assembling a complete big A matrix.
  % This model can be messy and is highly optimized for the code to run faster.
  % This model works only for 3D crack.
  
  properties
    cond;           %  material properties, including conduit geometry
    frac;             %  configurations for grid.
    u;                 %  unknown fields: [vz; pz; nz; h; pxy; vx; vy]
    ind_u;          %  the index of u.
    
    %                 du = ft(t) + Ai*u.
    Ai;                %  implicit part of matrix A.
    Ae;               %  explict part of matrix A less the terms involving dense matrix multiplification.
    ft;                 %  a function handle ft(t) = Ae*u + Fp(t), 
                        %       with both explicit terms and source time function.
                        
    Fp;               %  operators for external source from both the conduit and fracture. 
    E;                 %  energy norm 
  end
  
  methods
      %%
    function obj = coupledModel_NoMatrix(Mc, Mf)
      % construct a coupled model with a conduit and frac
      % build the conduit discretization
      obj.cond  = conduit(Mc); 
      
      % build the fracture discrization, terms that relates vx, vy to pxy 
      % are not included to avoid assemble a huge dense matrix.
      obj.frac    = frac3d(Mf);
      
      obj.ind_u = obj.fields_index();
      
      % build the discretization of the coupled model.
      [obj.Ai, obj.Ae, obj.Fp] = obj.build();
      
      % build the time integration function handle
      obj.ft = @(u, t) fun_integrate(obj, u, t);
      
      obj.u = zeros(sum(obj.dimensions()),1);
      
      obj.E = obj.energy_norm();
    end
    
    %%
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
        dim_c = obj.cond.dimensions();
        dim_f =  obj.frac.dimensions();
        dim = [sum(dim_c), sum(dim_f)];
    end

    function u = field(obj,U,num)
        %   [vz, pz, nz, h;
        %     pxy, vx, vy].
        %
        % if num = [1,1], return the vz in the conduit.
        % if num = [1], return all the unknowns in the conduit.
        if length(num)==1
          indices = field_indices(obj.dimensions(),num);
          u = U(indices);
        else if length(num) ==2
                indices1 = field_indices(obj.dimensions(),num(1));
                switch num(1)
                    case 1
                        u = obj.cond.field(U(indices1), num(2));
                    case 2 
                        u = obj.frac.field(U(indices1), num(2));
                end
            end
        end
    end
    
        function varargout = fields(obj, U)
            nvars_conduit = length(obj.cond.dimensions());
            vars_crack    = length(obj.frac.dimensions());
            nout = nvars_conduit + vars_crack;
            nz = obj.cond.geom.nz;
            nr = obj.cond.geom.nr;
            varargout = cell(nout, 1);
            
            for k = 1: nvars_conduit
                var = obj.field(U, [1, k]);
                if k ==1
                    var = reshape(var, nr, nz);
                end
                varargout{k} = var; 
            end
       
            varargout{nvars_conduit + 1} = obj.frac.geom.p.grd(obj.field(U,[2,1]));
            varargout{nvars_conduit + 2} = obj.frac.geom.vx.grd(obj.field(U,[2, 2]));
            varargout{nvars_conduit + 3} = obj.frac.geom.vy.grd(obj.field(U,[2, 3]));
        end
        
    function fields_index = fields_index(obj)
        % return the index of each field.
         
         dim_c = obj.cond.dimensions();
         dim_f  = obj.frac.dimensions(); 
         dims = [dim_c, dim_f];
         
         % fields in the conduit.
         fields_index.vz = [1: dims(1)];
         fields_index.pz = [(dims(1))+1: sum(dims(1:2))];
         fields_index.nz = [(sum(dims(1:2))+1): sum(dims(1:3))];
         fields_index.h = [(sum(dims(1:3))+1): sum(dims(1:4))];
         
         % fields in the crack
         fields_index.pxy = [(sum(dims(1:4))+1): sum(dims(1:5))];
         fields_index.vx = [(sum(dims(1:5))+1): sum(dims(1:6))];
         fields_index.vy = [(sum(dims(1:6))+1): sum(dims(1:7))];
           
    end
    
    function [Ai, Ae, Fp] = build(obj)
        % this method build the Ai, Ae, Fp from both conduit and
        % fluid-filled fractures.
        
        % Note that all the terms involving multiplication of dense matrix
        % K_t are not included in Ae and are handled separately.
        
          dim = obj.dimensions();
          
          Ae = block_matrix(dim, dim); 
          Ai = block_matrix(dim, dim);
          % construct the block diagonal terms
          Ae = block_matrix_insert(Ae,dim,dim,1,1, obj.cond.Ae);
          Ae = block_matrix_insert(Ae,dim,dim,2,2, obj.frac.Ae);
          Ai = block_matrix_insert(Ai,dim,dim,1,1, obj.cond.Ai);
          Ai = block_matrix_insert(Ai,dim,dim,2,2, obj.frac.Ai);
         
          % construct the forcing terms:
          
          % number of components. one for conduit and another one for frac.
          % This is hard coded for now. 
          nc = 2; 
           
          Fp = block_matrix(dim, ones(1, nc));
          Fp = block_matrix_insert(Fp,dim, [1 1], 1,1, obj.cond.Fp);
          Fp = block_matrix_insert(Fp,dim, [1 1], 2,2, obj.frac.Fp);
          
          % adding the coupling terms, all the coupling terms are explicit
         % in this case.
         Ae_c =  block_matrix(dim, dim, 0);
         
         dim_c = obj.cond.dimensions();
         dim_f  = obj.frac.dimensions(); 
         
         A12 = block_matrix(dim_c, dim_f, 0);
%          A21 = block_matrix(dim_f, dim_c, 0);
%          A22 = block_matrix(dim_f, dim_f, 0);
         
         % fracture to conduit.
         SAT = obj.cond.SAT;
         M_c   = obj.cond.M;
         op_c = obj.cond.op;
         A12_11 = SAT(1)/(M_c.rho(1)*M_c.c(1)) * kron(op_c.e0, op_c.er) * obj.frac.dc'; %[vz, px]
         A12_21 = SAT(1)*op_c.e0*obj.frac.dc'; %[pz, px]
         A12   = block_matrix_insert(A12, dim_c, dim_f, 1, 1, A12_11);
         A12   = block_matrix_insert(A12, dim_c, dim_f, 2, 1, A12_21);
         
         % conduit to fracture.
%          if ~isfield(obj.frac.M,'nz')
%              alpha = obj.cond.M.S/obj.frac.M.S;
%          else
%              alpha = obj.cond.M.S;
%          end
%          Kt = obj.frac.M.K_t;
%           if ~isfield(obj.frac.M,'nz')
%              H = obj.frac.op.p.Px;
%           else
%              H = obj.frac.op.p.Pxy2;
%           end
%          d = obj.frac.dc;
         
         % all the terms involving multification of dense Kt need to
         % handled separately.
         
%          A21_11 = - alpha* Kt*inv(H)*d*op_c.e0'*op_c.W2;
%          A21_12 = - alpha* Kt*inv(H)*d*(-op_c.e0'/(M_c.rho(1)*M_c.c(1)));
%          A21  = block_matrix_insert(A21, dim_f, dim_c, 1 , 1, A21_11);
%          A21  = block_matrix_insert(A21, dim_f, dim_c, 1 , 2, A21_12);
         
         Ae_c = block_matrix_insert(Ae_c,dim, dim, 1, 2, A12);
%          Ae_c = block_matrix_insert(Ae_c,dim, dim, 2, 1, A21);
         
         % fracture to fracture.
%          A22_11 = - alpha* Kt*inv(H)*d*(d'/(M_c.rho(1)*M_c.c(1)));
%          A22 = block_matrix_insert(A22, dim_f, dim_f, 1,1, A22_11);
%          Ae_c = block_matrix_insert(Ae_c,dim, dim, 2, 2, A22);
         
         Ae = Ae + Ae_c;
    end
    
    function obj = update(obj, U)
       obj.u = U;
       obj.cond = obj.cond.update(obj.field(U,1));
       obj.frac = obj.frac.update(obj.field(U,2));
    end
    
    function ft = fun_integrate(obj, u, t)
        % defines the integration function ft(u, t), additional parameters are
        % parsed in through obj.

        % contribution less the terms involving dense matrix.
        ft = obj.Ae*u + obj.Fp(:,1)*obj.cond.M.G(t) + obj.Fp(:,2)*obj.frac.M.G(t);

        % add contribution to pxy, from vx and vy.
        ind_vz =   obj.ind_u.vz;
        ind_pz =   obj.ind_u.pz;
        ind_pxy =   obj.ind_u.pxy;
        ind_vx =   obj.ind_u.vx;
        ind_vy =   obj.ind_u.vy;
        
        Kt = obj.frac.M.K_t;
        alpha = obj.cond.M.S;
        H_inv = inv(obj.frac.op.p.Pxy2);
        d = obj.frac.dc;
        rho1 = obj.cond.M.rho(1);
        c1    = obj.cond.M.c(1);
        e0      = obj.cond.op.e0;
         
%         ft(ind_pxy) = ft(ind_pxy) - Kt * (...
%                                 obj.frac.op.ux.Dx2*obj.frac.op.vx.Wz3*u(ind_vx) ...
%                              + obj.frac.op.uy.Dy2*obj.frac.op.vy.Wz3*u(ind_vy) ...
%                              + alpha * inv(H)*d*e0'*obj.cond.op.W2*u(ind_vz)...
%                              + alpha * inv(H)*d*(-e0'/(rho1 * c1)) * u(ind_pz)...
%                              + alpha * inv(H)*d*(d'/(rho1 * c1)) * u(ind_pxy));

        ft(ind_pxy) = ft(ind_pxy) - Kt * (...
                                obj.frac.op.ux.Dx2*(obj.frac.op.vx.Wz3*u(ind_vx)) ...
                             + obj.frac.op.uy.Dy2*(obj.frac.op.vy.Wz3*u(ind_vy)) ...
                             + alpha * H_inv*(d*(e0'*(obj.cond.op.W2*u(ind_vz))))...
                             + alpha * H_inv*(d*((-e0' * u(ind_pz))/(rho1 * c1)))...
                             + alpha * H_inv*(d*(d' * u(ind_pxy)/(rho1 * c1))));
    end

    function E = energy_norm(obj)

      dim = obj.dimensions();
      E    = block_matrix(dim,dim,1);
      % TODO: Add each term in the mechanical energy balance for the coupled conduit and crack system.
      E    = block_matrix_insert(E, dim,dim, 1,1, obj.cond.E);
      E    = block_matrix_insert(E, dim,dim, 2,2, obj.frac.E);
    end

  end

end
