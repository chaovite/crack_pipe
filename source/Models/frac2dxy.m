classdef frac2dxy
    % fluid filled fracture has two length scale and one width.
    % assume Poiseuille flow for viscosity.
    % right now, the material properties, rho, c, mu are treated as
    % constant but requires trivial work to make them vary in space.
    %
    
    properties
        M;              % material properties, including crack geometry, grid configuration.
        op;             % operators
        grd;            % grids. staggered.
        u;               % unknowns in the crack, p, vx, vy
        Fp;             % operators for external forcing or source. (all zero in this case)
        dc;             % the discretized delta function of at the coupling location.
        Ae;             % the part of A that treated explicitly, 
                           % do not include A12, A13, which includes dense
                          % matrix multiplification.
        indu;           % field index, a structure contains indp, indvx, indvy
        Ai;              % the part of A that treated implicitly.
        disloc3d;    % the 3d dislocation model for the crack.
        kernel_fft;   % the elastic kernel for fft.
        mask_fft;    % the mask for pressure for fft.
        E;               % energy norm
        
        % M:
        % name: name of the interface in the pipe this crack is coupled to.
        % geometry and grid configurations:
        % w0, Lx, Ly, nx, ny, order, interp_order, 
        % xs, ys (source location in crack),
        % xc, yc (coupling location between conduit and crack.)
        % lc (the characteristic length coupled to conduit)
        % isrigid, r_g, r_bl (grid stretch parameters)
        % fluid properties: K, mu, rho
        % solid properties: G; nu
    end
    
    methods
        
        function obj = frac2dxy(M)
            if ~isfield(M, 'isrigid')
                M.isrigid = true;
            end
            % construct the grid.
            [obj.grd, obj.op] = grids_frac2dxy(M.nx, M.ny, M.Lx, M.Ly, M.order);
            obj.indu = fields_index(obj);
%             tic;
%             disp('creating dislocation model....');
            if M.isrigid
                M.K_t      = M.K;
            else
                obj.disloc3d = disloc3dGrid(M.Xc, M.Yc, M.Zc, M.strike, M.dip, obj.grd.p.X, obj.grd.p.Y, M.Gs, M.nu);
                M.Ks   = obj.disloc3d.K;
                M.Ks_inv = obj.disloc3d.K_inv;
                M.K_t      = inv(1/M.K * speye(size(M.Ks_inv)) + M.Ks_inv/M.w0);
            end
%             toc;
            obj.M = M;
            
            if isfield(M, 'use_fft') && M.use_fft
                [obj.kernel_fft, obj.mask_fft]= obj.fft_kernel(M.npad_fft);
            end
            
            dim = obj.dimensions();
            [obj.Ae, obj.Ai] = obj.discretization();
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
            obj = obj.init();
            obj.E = obj.energy_norm();
        end
        
        function obj = init(obj)
            dim = obj.dimensions();
            obj.u = zeros(sum(dim), 1);
        end
        
        function dim = dimensions(obj)
            dim = [obj.grd.p.nx*obj.grd.p.ny,...
                      obj.grd.vx.nx*obj.grd.vx.ny,...
                      obj.grd.vy.nx*obj.grd.vy.ny];
        end
        
       function obj = update(obj, u1)
           obj.u = u1;
       end
        
       function ind = fields_index(obj)
%     dimension for p, vx, vy.
        dim   = obj.dimensions();
        ind.indp  = [1:dim(1)]';
        ind.indvx  = [(dim(1)+1): sum(dim(1:2))]';
        ind.indvy  = [sum(dim(1:2)) + 1: sum(dim(1:3))]';
       end
        
        function u = field(obj, U, num)
            indices = field_indices(obj.dimensions(),num);
            u = U(indices);
        end
        
        function [e, d]= delta(obj, xs, ys, interp_order)
            % Source term
            dim  = obj.dimensions();
            d = delta2d(obj.grd.p.x,xs, obj.grd.p.y,ys, interp_order); % weights of discretization of the delta function.

            % Scale by material properties and grid spacing
            es = obj.M.K_t*(obj.op.p.Pxy2\d)/obj.M.w0;
            e = block_matrix(dim,1,0);
            e = block_matrix_insert(e,dim,1,1,1,es);
        end
        
        function [Ae, Ai] = discretization(obj)
            % unknown pxy, vx, vy.
            dim = obj.dimensions();

            % implicit part.
            w0     = obj.M.w0;
            rhoi = 1/obj.M.rho;
            mu  = obj.M.mu;
            alpha = 12*mu/w0^2; % drag coefficient.
            
            Ai = spdiags([zeros(dim(1), 1); ...
                                 -alpha*rhoi*ones(sum(dim(2:end)), 1)], 0, sum(dim), sum(dim));
                             % explicit part.
            Ae = block_matrix(dim,dim,0);
            
             % Wave propagation
            % A12 and A13 and coupling terms are addressed separatly.
            % to avoid assemble huge dense matrix.
%             Ae12 = -obj.M.K_t * obj.op.vx.Dx2;
%             Ae13 = -obj.M.K_t * obj.op.vy.Dy2;
            Ae21 = -rhoi*obj.op.p.Dx2;
            Ae31 = -rhoi*obj.op.p.Dy2;
            Ae = block_matrix_insert(Ae,dim,dim,2,1,  Ae21);
            Ae = block_matrix_insert(Ae,dim,dim,3,1,  Ae31); 
%             Ae = block_matrix_insert(Ae,dim,dim,1,2,  Ae12);
%             Ae = block_matrix_insert(Ae,dim,dim,1,3,  Ae13);
        end
        
        function A = getA(obj)
            % return the full A matrix.
            A_sparse = obj.Ai + obj.Ae;
            % insert the part has dense Kt matrix.
            % A12, A13
            dim = obj.dimensions();
            A_dense = block_matrix(dim,dim,0);
            Ae12 = -obj.M.K_t * obj.op.vx.Dx2;
            Ae13 = -obj.M.K_t * obj.op.vy.Dy2;
            A_dense = block_matrix_insert(A_dense,dim,dim,1, 2,  Ae12);
            A_dense = block_matrix_insert(A_dense,dim,dim,1, 3,  Ae13);
            A = A_sparse + A_dense;
        end
        
        function [cmax, hmin] = getCFL(obj)
            % get min(c) and max(dz) for CFL condition.
            cmax = max(obj.M.c);
            hmin  = min([obj.grd.p.hx, obj.grd.p.hy, ...
                                obj.grd.vx.hx, obj.grd.vx.hy,...
                                obj.grd.vy.hx, obj.grd.vy.hy]);
        end
        
        function [kernel, mask]= fft_kernel(obj, N)
            % build the fft kernel that takes dp = ifft(kernel.*fft(q));
            hx = obj.grd.p.hx;
            hy = obj.grd.p.hy;
            nx = obj.grd.p.nx;
            ny = obj.grd.p.ny;
            % fft is the fastest when nx, ny are power of 2.
            
            % N must be even
            if mod(N, 2) ~= 0
                error('max number of grid points after padding must be even!')
            end
            
            if N > nx && N>ny
                % zero padding nx and ny.
                % x -> ndx1 [nx] ndx2
                % y -> ndy1 [ny] ndy2
                %
                ndx1 = floor((N - nx)/2);
                ndy1 = floor((N - ny)/2);
                mask = zeros(N, N);
                mask(ndx1+1:ndx1+nx, ndy1+1:ndy1+ny) = 1;
                mask = logical(mask);
            end
            
            % create this fft kernel.
            % wave number.
            dkx = 2*pi/hx/N;
            kx   = zeros(1, N);
            kx(1:N/2 + 1) = [0: N/2]*dkx; 
            kx(N:-1:N/2+2) = -kx(2 :N/2);
            dky = 2*pi/hy/N;
            ky   = zeros(1, N);
            ky(1:N/2 + 1) = [0: N/2]*dky; 
            ky(N:-1:N/2+2) = -ky(2 :N/2);
            [KX, KY] = meshgrid(kx, ky);
            Kf = obj.M.K;
            G = obj.M.Gs;
            w0 = obj.M.w0;
            nu = obj.M.nu;
            kernel = 1./(1/Kf + (1-nu)/(G*w0)./sqrt(KX.^2 + KY.^2));
        end
        
        function U = eval_disp_p(obj, p, xq, yq)
            % evaluate the surface displacement a list of query points, xq, yq.
            % p should be the same size as pressure unknowns on this crack
            % xq, yq: [nq*1] vectors.
            % U: [3, nq]
            %
                U = obj.disloc3d.eval_disp_p(p, xq, yq);
        end
        
        function E = energy_norm(obj)
            % TO DO
            % Energy norm, E*A+A'*E
            dim = obj.dimensions();
            E    = block_matrix(dim,dim,1);
            E    = block_matrix_insert(E,dim,dim,1,1,obj.M.w0*inv(obj.M.K_t)*obj.op.p.Pxy2);
            E    = block_matrix_insert(E,dim,dim,2,2,obj.M.w0*obj.M.rho*obj.op.vx.Pxy2);
            E    = block_matrix_insert(E,dim,dim,3,3,obj.M.w0*obj.M.rho*obj.op.vy.Pxy2);
        end
        
    end
    
end
