classdef fluid
  % acoustic
  properties
    metrics
    op
    geom
    A
    E
    B
    Ae;
    Ai;
    u
    w0;
    % fluid properties
    rho;
    K;
    mu;
    % solid properties
    G; % shear modulus
    nu; %Poisson ratio
    Ks; 
    Ks_inv;
    K_t;
  end

  methods
    function obj = fluid(geom,metrics,op,w0,rho,K,mu, G, nu, isrigid)
        if nargin<10
            isrigid = false;
        end
      obj.metrics  = metrics;
      obj.op       = op;
      obj.geom     = geom;
      obj.w0       = w0;
      obj.rho      = rho;
      obj.K        = K;
      obj.mu       = mu;
      obj.G         = G;
      obj.nu        = nu;
      if isrigid
          obj.K_t      = K;
      else
          [obj.Ks, obj.Ks_inv] = KernelBEM2D(geom.pp.x, obj.G, obj.nu);
          obj.K_t      = inv(1/K * speye(op.pp.n1, op.pp.n1) + obj.Ks_inv/obj.w0);
      end
      obj.u = zeros(sum(obj.dimensions()),1);
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
      dim = [obj.op.p.n1,...
             obj.op.m.n1*obj.op.mm.n2,...
        ];

    end

    function u = field(obj,U,num)

      indices = field_indices(obj.dimensions(),num);
      u = U(indices);
   
    end

    function e = source(obj,xs,order)

      % Source term
      dim  = obj.dimensions();
      es = spalloc(dim(1),1,1);
      x = obj.geom.p.x;
      [val index] = min(abs(x-xs));
      range = index  + (-order/2:order/2);
      range = max(min(range), 1): max(range);
      l = interp(x(range),xs);
      es(index) = 1;
      es(range) = double(l);
      % Scale by material properties and grid spacing
      es = obj.K_t*inv(obj.op.p.Px)*es;
      e = block_matrix(dim,1,0);
      e = block_matrix_insert(e,dim,1,1,1,es);
    end


    function [Ae, Ai] = interior(obj)
             
      dim = obj.dimensions();
      op = obj.op;

      Ae = block_matrix(dim,dim,0);
      Ai = block_matrix(dim,dim,0);

      e   = ones(op.mm.n2,1);

%       K    = obj.K;
      rhoi = 1/obj.rho;

     % TODO: load variable coefficients, stretching the grid in y direction.
     
     Ix = speye(op.m.n1, op.m.n1);
     
     J_eta_ypi  = obj.geom.p.Jyi; % Jacobian matrix yp_eta
     J_eta_ymi = obj.geom.m.Jyi; % Jacobian matrix ym_eta
     J_eta_ym = obj.geom.m.Jy;% Jacobian matrix eta_ym
      
     % Wave propagation
     Ae12 = -obj.K_t*kron(op.pm.Dx1, 1/obj.w0*e'*J_eta_ym*op.mm.Py);
     Ae = block_matrix_insert(Ae,dim,dim,1,2, Ae12);
     Ae = block_matrix_insert(Ae,dim,dim,2,1,-rhoi*kron(op.mm.Dx1,e));

     % Diffusion
     D2 = spalloc(dim(2),dim(2),dim(2));
     D2 = rhoi*obj.mu *kron(Ix, J_eta_ymi) * op.mm.Dy * kron(Ix, J_eta_ypi)* op.mp.Dy ;
     Ai = block_matrix_insert(Ai, dim, dim, 2, 2, D2);
      
    end

    function E = energy_norm(obj)

      dim = obj.dimensions();
      E    = block_matrix(dim,dim,1);
      % TODO: Add each term in the mechanical energy balance
      E    = block_matrix_insert(E,dim,dim,1,1,obj.w0*inv(obj.K_t)*obj.op.pp.Px);
      E    = block_matrix_insert(E,dim,dim,2,2,obj.rho*kron(obj.op.mm.Px, obj.geom.m.Jy * obj.op.mm.Py));
    end


  end

end
