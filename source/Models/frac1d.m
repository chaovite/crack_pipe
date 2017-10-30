classdef frac1d
    % 1d fluid filled fracture 
    % 2D plane strain for elasticity, 1D width averaged formulation for fluid dynamics
    % 1D linear acoustics + static elasticity.
    % Weak enforcement of b.c. at both ends.
    
    properties
        M;              % material properties, including crack geometry, grid configuration.
        op;             % operators
        geom;        % grids. staggered.
        u;               % unknowns in the crack, px, vx.
        A;               % the du/dt = A*u+Fp*G.
        Fp;             % operators for external forcing or source. (all zero in this case)
        dc;             % the discretized delta function of at the coupling location.
        Ae;             % the part of A that treated explicitly.
        Ai;              % the part of A that treated implicitly.
        E;               % energy norm
        
        % M:
        % geometry and grid configurations:
        % w0, Lx, nx, order, interp_order, 
        % xs, ys (source location in crack),
        % xc, yc (coupling location between conduit and crack.)
        % bc: a structure defines the boundary conditions, v0, p0, vn, pn.
        % lc (the characteristic length coupled to conduit)
        % isrigid
        % fluid properties: K, mu, rho
        % solid properties: G; nu
    end
    
    methods
        
        function obj = frac1d(M)
            if ~isfield(M, 'isrigid')
                M.isrigid = true;
            end
            % construct the grid.
            [obj.geom, obj.op]  = grids_frac1d(M.nx, M.Lx, M.order, 'weak');
            
            if M.isrigid
                M.K_t      = M.K;
            else
                [M.Ks, M.Ks_inv] = KernelBEM2D(obj.geom.p.x, M.Gs, M.nu);
                M.K_t      = inv(1/M.K * speye(obj.geom.p.nx, obj.geom.p.nx) + M.Ks_inv/M.w0);
            end
            
            obj.M = M;
            
            % construct the discretization. A, Ae, Ai, Fp.
            % Note that the coupling terms that require the knowledge of
            % conduit is not included here. These  terms will be added when
            % constructing the coupled model.
            dim = obj.dimensions();
            obj.Fp = struct();
            
            [obj.Ae, obj.Ai, obj.Fp.bc] = obj.discretization();
            
            obj.A = obj.Ae + obj.Ai;
            obj.Fp.ps  = block_matrix(dim,1, 0);% point source.
            obj.dc  = block_matrix(dim,1, 0);
            
            % delta function source in the crack.
            if isfield(M, 'xs')
                [obj.Fp.ps, ~ ] = obj.delta(M.xs, M.interp_order);
            end
            
            % delta functon, coupled to the conduit
            if isfield(M, 'xc')
                [~, obj.dc] = obj.delta(M.xc, M.interp_order);
            end
            
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
            dim = [obj.geom.p.nx, obj.geom.v.nx];
        end
        
       function obj = update(obj,U)
           obj.u = U;
        end
        
        function u = field(obj,U,num)
            indices = field_indices(obj.dimensions(), num);
            u = U(indices);
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
        
        
        function [Ae, Ai, Fp_bc] = discretization(obj)
            % order of unknown: p and v.
            dim = obj.dimensions();
            
            Ae = block_matrix(dim,dim,0);
            Ai = block_matrix(dim,dim,0);
            Fp_bc = struct();
            rhoi = 1/obj.M.rho;
            mu  = obj.M.mu;
            k     = mu/(obj.M.w0^2/12);% friction parameter. mu/(w^2/12)

            % PDE contribution
            Ae12 = -obj.M.K_t * obj.op.p.Dx;
            Ae21 = -rhoi*obj.op.v.Dx;
            Ai22 = -k*rhoi*speye(dim(2));% the dissipation term.
            
            e0_p = obj.op.p.restrictions.e0;
            en_p = obj.op.p.restrictions.en;
            e0_v = obj.op.v.restrictions.e0;
            en_v = obj.op.v.restrictions.en;
            H_p  = obj.op.p.Px;
            H_v  = obj.op.v.Px; 

            % Boundary conditions:
            bc_types = fields(obj.M.bc);
            
            for i = 1: length(bc_types)
                type = bc_types{i};
                switch type
                    case 'p0'
                        Ae21 = Ae21 - rhoi*inv(H_v)*e0_v*e0_p';
                        Fp_bc.p0 =  block_matrix(dim,1, 0);
                        Fp_bc.p0 = block_matrix_insert(Fp_bc.p0,dim,1, 2, 1, rhoi*inv(H_v)*e0_v);
                    case 'v0'
                        Ae12 = Ae12 - obj.M.K_t*inv(H_p)*e0_p*e0_v';
                        Fp_bc.v0 =  block_matrix(dim,1, 0);
                        Fp_bc.v0 = block_matrix_insert(Fp_bc.v0,dim,1, 1, 1,  obj.M.K_t*inv(H_p)*e0_p);
                    case 'pn'
                        Ae21 = Ae21 + rhoi*inv(H_v)*en_v*en_p';
                        Fp_bc.pn =  block_matrix(dim,1, 0);
                        Fp_bc.pn = block_matrix_insert(Fp_bc.pn,dim,1, 2, 1, - rhoi*inv(H_v)*en_v);
                    case 'vn'
                        Ae12 = Ae12 + obj.M.K_t*inv(H_p)*en_p*en_v';
                        Fp_bc.vn =  block_matrix(dim,1, 0);
                        Fp_bc.vn = block_matrix_insert(Fp_bc.vn,dim,1, 1, 1, - obj.M.K_t*inv(H_p)*en_p);
                end
            end
            
            Ae = block_matrix_insert(Ae,dim,dim,1,2,  Ae12);
            Ae = block_matrix_insert(Ae,dim,dim,2,1,  Ae21);
            Ai = block_matrix_insert(Ai, dim, dim, 2, 2, Ai22);
            
        end
        
        function E = energy_norm(obj)
            % Energy norm, E*A+A'*E
            dim = obj.dimensions();
            E    = block_matrix(dim,dim,1);
            E    = block_matrix_insert(E,dim,dim,1,1, obj.M.w0*inv(obj.M.K_t)*obj.op.p.Px);
            E    = block_matrix_insert(E,dim,dim,2,2, obj.M.w0*obj.M.rho*obj.op.v.Px);
        end
        
    end
    
end
