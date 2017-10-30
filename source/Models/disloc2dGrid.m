classdef disloc2dGrid
    
    properties
        
        % Xc, Yc, Zc: the global coordinate of fracture centroid, in East, North, Depth.
        Xc, Yc, Zc;
        
        % L, D: the length and width of the entire fracture.
        
        % strike, dip: strike and dip of the fracture. 0<strike<360, 0<dip<90.
        strike, dip;
        
        x;
        % x: local x coordinates of each fracture cell on the grid in
        %         x in the direction of positive dip, y in the direction of
        %         strike, with the corner as the origin. 
        %        dimension: [Nx*1]
        
        L, W; % length (strike direction (y direction)) and width (dip direction) 

        %Xm Ym Zm: the global coordinate of reference point on each fracture cell of the grid in East, North, Depth
        % the reference point is defined same as disloc3d, the center point
        % of the bottom edge. For the fracture model.
        Xm, Ym, Zm;
        
        X, Y, Z; % global coordinates for fracture cell centers.
        
        % mu, nu: shear modulus and
        mu, nu;
        
        % K relative fracture opening to normal stress. p=  K * w
        K, K_inv;
        
    end
    
    methods
        function obj = disloc2dGrid(Xc, Yc, Zc, strike, dip, x, Ly, mu, nu)
            % construct dislocation model from grid.
            obj.Xc = Xc; obj.Yc = Yc; obj.Zc = Zc; 
            obj.strike = strike; obj.dip = dip; 
            obj.x = x;
            obj.L = Ly;
            obj.W = x(end);
            obj.mu = mu;
            obj.nu  = nu; 
            hx   = x(2)-x(1);

            [obj.Xm, obj.Ym, obj.Zm] = obj.local2global(x + hx/2, obj.W/2);
            [obj.X, obj.Y, obj.Z] = obj.local2global(x, obj.W/2);
            [obj.K, obj.K_inv] = obj.ddm_matrix();
        end
        
        function [X, Y, Z] = local2global(obj, x, xc)
            % convert local x,y coordinate (on fracture plane) to global coordinates X, Y, Z
            % xc, yc is the centroid of the fracture plane.
            
            phi = pi*obj.strike/180;
            theta = pi*obj.dip/180;
            
            X = (x-xc)*cos(phi)*cos(theta) + obj.Xc;
            Y = -(x-xc)*sin(phi)*cos(theta) + obj.Yc;
            Z = (x-xc)*sin(theta) + obj.Zc;
            
        end
        
        function [K, K_inv] = ddm_matrix(obj)
            %call disloc3d to build ddm matrix K, pressure p = K * w.
            n   = [sind(obj.strike + 90)*sind(obj.dip), ...
                         cosd(obj.strike+90)*sind(obj.dip), ...
                         cosd(obj.dip)]'; % unit normal of fault plane.
            
            N  = numel(obj.X);
            
            % the dimension for each fracture cell.
            nx = size(obj.X, 1) - 1;
            
            % the length and width for each fracture cell.
            Ly = obj.L;
            width  = obj.W/nx;
            
            % create observation at each cell center. 
            obs = [reshape(obj.X, N, 1)'; reshape(obj.Y, N, 1)'; - reshape(obj.Z, N, 1)'];
            
            K = zeros(N, N);
            for i = 1:N
                % create model for fracture cell i.
                depth = obj.Zm(i);
                east   = obj.Xm(i);
                north  = obj.Ym(i);
                % create a model with unit opening at fracture cell i
                mdl = [Ly, width, depth, obj.dip, obj.strike, east, north, 0, 0, 1]'; 
                % get the stress.
                [~, ~, S, ~] = disloc3d(mdl,obs, obj.mu,obj.nu);
                % Sxx, Sxy, Sxz, Syy, Syz, and Szz
                K(:, i) = -(S(1,:)*n(1)^2 + 2*S(2,:)*n(1)*n(2) + S(3,:)*n(1)*n(3) ...
                    + S(4,:)*n(2)^2 + S(5,:)*n(2)*n(3) + S(6,:)*n(3)^2); 
            end
            K_inv=inv(K);
        end
        
        function U = eval_disp_p(obj, p, xq, yq)
            % given solution of pressure distribution on the fault and query points location on the surface, evaluate the
            % surface displacement U = (Ux, Uy, Uz), [east, north, up].  disloc3d will sum up the contribution from each fracture
            % element for each query points.
            
            % U: 3*nq.
            
            % p, pressure at the center of each fracture cell.
            N = numel(p);
            Nq = numel(xq);
            
            if size(p,2) ~= 1
                p = p(:);
            end
            
            w = (obj.K_inv*p)'; % opening.
            
            % the dimension for each fracture cell.
            nx = size(obj.X, 1) - 1;
            ny = 1;
            
            % the length and width for each fracture cell.
            Ly = obj.L/ny;
            Lx  = obj.W/nx;
            
            % create fracture models:
            mdl = [Ly*ones(1,N); Lx*ones(1,N);
                        obj.Zm(:)'; obj.dip*ones(1,N); obj.strike*ones(1,N); 
                        obj.Xm(:)'; obj.Ym(:)'; zeros(1,N); zeros(1,N); w];
            
            obs = [xq(:)'; yq(:)'; zeros(1, Nq)];
            
            [U, ~, ~, ~] = disloc3d(mdl,obs, obj.mu,obj.nu);
           
        end
        
        function U = eval_disp_w(obj, w, xq, yq)
            % given solution of openning distribution on the fault and
            % query points location on the surface, evaluate the
            % surface displacement U = (Ux, Uy, Uz), [east, north, up]
                        % U: 3*nq.
            
            % p, pressure at the center of each fracture cell.
            N = numel(w);
            Nq = numel(xq);
            
            if size(w,1) ~= 1
                w = w(:)';
            end

            % the dimension for each fracture cell.
            nx = size(obj.X, 1) - 1;
            ny = 1;
            
            % the length and width for each fracture cell.
            Ly = obj.L/ny;
            Lx  = obj.W/nx;
            
            % create fracture models:
            mdl = [Ly*ones(1,N); Lx*ones(1,N);
                        obj.Zm(:)'; obj.dip*ones(1,N); obj.strike*ones(1,N); 
                        obj.Xm(:)'; obj.Ym(:)'; zeros(1,N); zeros(1,N); w];
            
            obs = [xq(:)'; yq(:)'; zeros(1, Nq)];
            
            [U, ~, ~, ~] = disloc3d(mdl,obs, obj.mu,obj.nu);
           
        end
        
    end
    
end
