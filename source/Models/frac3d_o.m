classdef frac3d_o
    % 3d fluid filled fracture, orginal version
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
        disloc3d;    % the 3d dislocation model for the crack.
        E;               % energy norm
        
        % M:
        % geometry and grid configurations:
        % w0, Lx, Ly, nx, ny,nz, order, interp_order, 
        % xs, ys (source location in crack),
        % xc, yc (coupling location between conduit and crack.)
        % lc (the characteristic length coupled to conduit)
        % isrigid, r_g, r_bl (grid stretch parameters)
        % fluid properties: K, mu, rho
        % solid properties: G; nu
    end
    
    methods
        
        function obj = frac3d_o(M)
            if ~isfield(M, 'isrigid')
                M.isrigid = false;
            end
            % construct the grid.
            [obj.geom, obj.op] = grids_frac3d(M.nx,M.ny,M.nz, M.Lx, M.Ly, M.w0, M.order, 'strong',M.r_g, M.r_bl);
            
            if M.isrigid
                M.K_t      = M.K;
            else
                % TO DO: implement 3D DDM to calculate Ks.
                obj.disloc3d = disloc3dGrid(M.Xc, M.Yc, M.Zc, M.strike, M.dip, obj.geom.p.X, obj.geom.p.Y, M.Gs, M.nu);
                M.Ks   = obj.disloc3d.K;
                M.Ks_inv = obj.disloc3d.K_inv;
                M.K_t      = inv(1/M.K * speye(size(M.Ks_inv)) + M.Ks_inv/M.w0);
            end

            obj.M = M;
            
            % construct the discretization. A, Ae, Ai, Fp.
            % Note that the coupling terms that require the knowledge of
            % conduit is not included here. These  terms will be added when
            % constructing the coupled model.
            dim = obj.dimensions();
            obj.E = obj.energy_norm();
            
            [obj.Ae, obj.Ai] = obj.discretization();
            
            obj.A = obj.Ae + obj.Ai;
            obj.Fp  = block_matrix(dim,1, 0);
            obj.dc  = block_matrix(dim,1, 0);
            
            % delta function source in the crack.
            if isfield(M, 'xs') && isfield(M, 'ys')
                [obj.Fp, ~ ] = obj.delta(M.xs, M.ys, M.interp_order);
            end
            
            % delta functon, coupled to the conduit
            if isfield(M, 'xc') && isfield(M, 'yc')
                [~, obj.dc] = obj.delta(M.xc, M.yc, M.interp_order);
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
            dim = [obj.geom.p.nx*obj.geom.p.ny,...
                      obj.geom.vx.nx*obj.geom.vx.ny*obj.geom.vx.nz,...
                      obj.geom.vy.nx*obj.geom.vy.ny*obj.geom.vy.nz];
        end
        
       function obj = update(obj,U)
           obj.u = U;
        end
        
        function u = field(obj,U,num)
            indices = field_indices(obj.dimensions(),num);
            u = U(indices);
        end
        
        function [e, d]= delta(obj, xs, ys, interp_order)
            % Source term
            dim  = obj.dimensions();
            d = delta2d(obj.geom.p.x,xs, obj.geom.p.y,ys, interp_order); % weights of discretization of the delta function.

            % Scale by material properties and grid spacing
            es = obj.M.K_t*(obj.op.p.Pxy2\d);
            e = block_matrix(dim,1,0);
            e = block_matrix_insert(e,dim,1,1,1,es);
        end
        
        function [Ae, Ai] = discretization(obj)
            
            dim = obj.dimensions();
            
            Ae = block_matrix(dim,dim,0);
            Ai = block_matrix(dim,dim,0);
            
            rhoi = 1/obj.M.rho;
            mu  = obj.M.mu;

            % Wave propagation
            Ae12 = -obj.M.K_t * obj.op.ux.Dx2*obj.op.vx.Wz3;
            Ae13 = -obj.M.K_t * obj.op.uy.Dy2*obj.op.vy.Wz3;
            Ae21 = -rhoi*obj.op.p.Dx3*obj.op.p.Ez;
            Ai22 = mu*rhoi*obj.op.vx.D2z3;
            Ae31 = -rhoi*obj.op.p.Dy3*obj.op.p.Ez;
            Ai33 = mu*rhoi*obj.op.vy.D2z3;

            Ae = block_matrix_insert(Ae,dim,dim,1,2, Ae12);
            Ae = block_matrix_insert(Ae,dim,dim,1,3, Ae13);
            Ae = block_matrix_insert(Ae,dim,dim,2,1, Ae21);
            Ae = block_matrix_insert(Ae,dim,dim,3,1, Ae31);
            Ai = block_matrix_insert(Ai, dim, dim, 2, 2, Ai22);
            Ai = block_matrix_insert(Ai, dim, dim, 3, 3, Ai33);
        end
        
        function E = energy_norm(obj)
            % Energy norm, E*A+A'*E
            dim = obj.dimensions();
            E    = block_matrix(dim,dim,1);
            E    = block_matrix_insert(E,dim,dim,1,1,obj.M.w0*inv(obj.M.K_t)*obj.op.p.Pxy2);
            E    = block_matrix_insert(E,dim,dim,2,2,obj.M.rho*obj.op.vx.Pxyz3);
            E    = block_matrix_insert(E,dim,dim,3,3,obj.M.rho*obj.op.vy.Pxyz3);
        end
        
    end
    
end
