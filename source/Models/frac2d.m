classdef frac2d
    % 2d fluid filled fracture
    properties
        M;              % material properties, including crack geometry, grid configuration.
        op;             % difference operators
        geom;        % grids. staggered.
        u;               % unknowns in the crack, px, vx.
        A;               % the du/dt = A*u+Fp*G.
        Fp;             % operators for external forcing or source. (all zero in this case)
        dc;             % the discretized delta function of at the coupling location.
        Ae;             % the part of A that treated explicitly.
        Ai;              % the part of A that treated implicitly.
        disloc2d;    % the 2D dislocation model.
        E;               % energy norm

        
        % M:
        % geometry and grid configurations:
        % w0, L, nx, ny,  order, interp_order, 
        % xs (source location in crack),
        % xc (coupling location between conduit and crack.)
        % lc (the characteristic length coupled to conduit)
        % S: fracture cross section area (w0*lc)
        % isrigid, r_g, r_bl (grid stretch parameters)
        % fluid properties: K, mu, rho
        % solid properties: G; nu
    end
    
    methods
        
        function obj = frac2d(M)
            if ~isfield(M, 'isrigid')
                M.isrigid = false;
            end
            % construct the grid.
            grid_types = {'pp','mm','pm','mp','p','m'};
            truncate = true;
            for i=1:length(grid_types)
                gt = grid_types{i};
                geom.(gt)    = stretched_grid(gt,M.order,'strong',M.nx,M.ny,M.L,M.w0,M.r_g,M.r_bl,truncate);
                op.(gt)         = operators(gt,M.order,geom.(gt),'strong');
            end
            
            obj.op       = op;
            obj.geom     = geom;
            disp('creating disloc2d model....');
            obj.disloc2d = disloc2dGrid(M.Xc, M.Yc, M.Zc, M.strike, M.dip, obj.geom.p.x, M.Lz, M.Gs, M.nu);
            
            if M.isrigid
                M.K_t      = M.K;
            else
                M.Ks   = obj.disloc2d.K;
                M.Ks_inv = obj.disloc2d.K_inv;
                M.K_t      = inv(1/M.K * speye(size(M.Ks_inv)) + M.Ks_inv/M.w0);
            end
            
            obj.M = M;
            
            % construct the discretization. A, Ae, Ai, Fp.
            % Note that the coupling terms that require the knowledge of
            % conduit is not included here. These  terms will be added when
            % constructing the coupled model.
            dim = obj.dimensions();
            [obj.Ae, obj.Ai] = obj.discretization();
            obj.A = obj.Ae + obj.Ai;
            obj.Fp  = block_matrix(dim,1, 0);
            
            if isfield(M, 'xs')
                [obj.Fp, ~ ] = obj.delta(M.xs, M.interp_order);% delta function source in the crack.
            end
            
            [~, obj.dc ] = obj.delta(M.xc, M.interp_order); % delta functon, coupled to the conduit.
            
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
        
        
        function dim = dimensions(obj)
            dim = [obj.op.p.n1,...
                obj.op.m.n1*obj.op.mm.n2,...
                ];
        end
        
       function obj = update(obj,U)
           obj.u = U;
        end
        
        function u = field(obj,U,num)
            
            indices = field_indices(obj.dimensions(),num);
            u = U(indices);
            
        end
        
        function w = get_w(obj)
            % given the crack opening.
            p  = obj.field(obj.u,1);
            w  = obj.M.Ks_inv*p;
        end
        
        function [e, d]= delta(obj, xs, interp_order)
            
            % Source term
            dim  = obj.dimensions();
            d = delta1d(obj.geom.p.x, xs, interp_order);% discretization of the delta function.
            
            % Scale by material properties and grid spacing
            es = obj.M.K_t*(obj.op.p.Px\d);
            e = block_matrix(dim,1,0);
            e = block_matrix_insert(e,dim,1,1,1,es);
        end
        
        
        function [Ae, Ai] = discretization(obj)
            
            dim = obj.dimensions();
            
            Ae = block_matrix(dim,dim,0);
            Ai = block_matrix(dim,dim,0);
            
            e   = ones(obj.op.mm.n2,1);
            
            rhoi = 1/obj.M.rho;
            
            Ix = speye(obj.op.m.n1, obj.op.m.n1);
            
            J_eta_ypi  = obj.geom.p.Jyi; % Jacobian matrix yp_eta
            J_eta_ymi = obj.geom.m.Jyi; % Jacobian matrix ym_eta
            J_eta_ym = obj.geom.m.Jy;% Jacobian matrix eta_ym
            
            % Wave propagation
            Ae12 = -obj.M.K_t*kron(obj.op.pm.Dx1, 1/obj.M.w0*e'*J_eta_ym*obj.op.mm.Py);
            Ae21 = -rhoi*kron(obj.op.mm.Dx1,e);
            Ae = block_matrix_insert(Ae,dim,dim,1,2, Ae12);
            Ae = block_matrix_insert(Ae,dim,dim,2,1, Ae21);
            
            % Diffusion
            D2 = spalloc(dim(2),dim(2),dim(2));
            D2 = rhoi*obj.M.mu *kron(Ix, J_eta_ymi) * obj.op.mm.Dy * kron(Ix, J_eta_ypi)* obj.op.mp.Dy ;
            Ai = block_matrix_insert(Ai, dim, dim, 2, 2, D2);
            
        end
        
        function E = energy_norm(obj)
            % Energy norm, E*A+A'*E.
            dim = obj.dimensions();
            E    = block_matrix(dim,dim,1);
            E    = block_matrix_insert(E,dim,dim,1,1,obj.M.w0*inv(obj.M.K_t)*obj.op.pp.Px);
            E    = block_matrix_insert(E,dim,dim,2,2,obj.M.rho*kron(obj.op.mm.Px, obj.geom.m.Jy * obj.op.mm.Py));
        end
        
    end
    
end
