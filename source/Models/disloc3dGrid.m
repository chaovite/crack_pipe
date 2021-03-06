classdef disloc3dGrid
    
    properties
        
        % Xc, Yc, Zc: the global coordinate of fracture centroid, in East, North, Depth.
        Xc, Yc, Zc;
        
        % L, D: the length and width of the entire fracture.
        
        % strike, dip: strike and dip of the fracture. 0<strike<360, 0<dip<90.
        strike, dip;
        
        % x, y: local x, y coordinates of each fracture cell on the grid in
        %         x in the direction of positive dip, y in the direction of
        %         strike, with the corner as the origin. 
        %        dimension: [Ny * Nx]
        x, y;
        
        L, W; % length (strike direction) and width (dip direction) 

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
        function obj = disloc3dGrid(Xc, Yc, Zc, strike, dip, x, y, mu, nu)
            % construct dislocation model from grid.
            obj.Xc = Xc; obj.Yc = Yc; obj.Zc = Zc; 
            obj.strike = strike; obj.dip = dip; 
            obj.x = x;
            obj.y = y;
            obj.L = y(end,1);
            obj.W = x(1,end);
            obj.mu = mu;
            obj.nu  = nu; 
            hx   = x(1,2)-x(1,1);

            [obj.Xm, obj.Ym, obj.Zm] = obj.local2global(x + hx/2, y, obj.W/2, obj.L/2);
            [obj.X, obj.Y, obj.Z] = obj.local2global( x, y, obj.W/2, obj.L/2);

            [obj.K, obj.K_inv] = obj.ddm_matrix(); 
            % this can be replaced with more efficient full space solution
            %TODO
%             [obj.K, obj.K_inv] = obj.ddm_matrix_fullspace();
        end
        
        function [X, Y, Z] = local2global(obj, x, y, xc, yc)
            % convert local x,y coordinate (on fracture plane) to global coordinates X, Y, Z
            % xc, yc is the centroid of the fracture plane.
            
            phi = pi*obj.strike/180;
            theta = pi*obj.dip/180;
            
            X = (x-xc)*cos(phi)*cos(theta) + (y-yc)*sin(phi) + obj.Xc;
            Y = -(x-xc)*sin(phi)*cos(theta) + (y-yc)*cos(phi)+ obj.Yc;
            Z = (x-xc)*sin(theta) + obj.Zc;
        end
        
        function [K, K_inv] = ddm_matrix(obj)
            % call disloc3d to build ddm matrix K, pressure p = K * w.
            % TODO: use full space elastic green's function for this part to speed up the calculation.
            
            n   = [sind(obj.strike + 90)*sind(obj.dip), ...
                         cosd(obj.strike+90)*sind(obj.dip), ...
                         cosd(obj.dip)]'; % unit normal of fault plane.
            
            N  = numel(obj.X);
            
            % the dimension for each fracture cell.
            nx = size(obj.X, 2) - 1;
            ny = size(obj.X, 1) - 1;
            
            % the length and width for each fracture cell.
            Ly = obj.L/ny;
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
                K(:, i) = -(S(1,:)*n(1)^2 + 2*S(2,:)*n(1)*n(2) + 2*S(3,:)*n(1)*n(3) ...
                    + S(4,:)*n(2)^2 + 2*S(5,:)*n(2)*n(3) + S(6,:)*n(3)^2); 
            end
            K_inv=inv(K);
        end
        
        function Hp = get_Hp(obj, xq, yq)
            % get the matrix Hp that relates the pressure to surface
            % displacement observed at xq and yq.
            % Hp size: [3*nq, np]
            Hw = obj.get_Hw(xq, yq);
            Hp = Hw*obj.K_inv;
        end
        
        
        function Hp = get_Hp_tilt(obj, xq, yq)
            % get the matrix Hp that relates the pressure to surface
            % tilt observed at xq and yq.
            % Hp size: [3*nq, np]
            Hw = obj.get_Hw_tilt(xq, yq);
            Hp = Hw*obj.K_inv;
        end
        
        function Hw_tilt = get_Hw_tilt(obj, xq, yq)
            % get the matrix Hw_tilt that relates the opening w to surface
            % tilt observed at xq and yq.
            % Hw_tilt size: [3*nq, nw]
            % Hw_tilt(i, j) is the amount of tilt component i given unit
            % dislocation at fault j.

            N = numel(obj.x); % number of fault elements.
            Nq = numel(xq);   % number of query points.
            Hw_tilt = zeros(Nq*3, N); 
            
            % we make use of the function obj.get_Hw
            Hwxp  = obj.get_Hw(xq+1, yq);  % x+1
            Hwxm = obj.get_Hw(xq-1, yq);   % x-1
            Hwyp  = obj.get_Hw(xq,yq+1);   % y+1
            Hwym = obj.get_Hw(xq,yq-1);    % y-1
            
            % compute tilt using 2nd order central difference.
            % tilt in x direction.
            Hw_tilt(1:3:end, :) =  (Hwxp(1:3:end, :) - Hwxm(1:3:end, :))/2;
            % tilt in y direction.
            Hw_tilt(2:3:end, :) = (Hwyp(2:3:end, :) - Hwym(2:3:end, :))/2;
            % tilt in z direction = 0
        end
        

        function Hw = get_Hw(obj, xq, yq)
            % get the matrix Hp that relates the opening w to surface
            % displacement observed at xq and yq.
            % Hp size: [3*nq, nw]
            % Hp(i, j) is the amount of displacement component i given unit
            % dislocation at fault j.
            
            N = numel(obj.x); % number of fault elements.
            Nq = numel(xq);
            
            w = ones(1, N);

            % the dimension for each fracture cell.
            nx = size(obj.X, 2) - 1;
            ny = size(obj.X, 1) - 1;
            
            % the length and width for each fracture cell.
            Ly = obj.L/ny;
            Lx  = obj.W/nx;
            Hw = zeros(Nq*3, N);
            
            % create fracture models:
            mdl = [Ly*ones(1,N); Lx*ones(1,N);
                        obj.Zm(:)'; obj.dip*ones(1,N); obj.strike*ones(1,N); 
                        obj.Xm(:)'; obj.Ym(:)'; zeros(1,N); zeros(1,N); w];
            
            obs = [xq(:)'; yq(:)'; zeros(1, Nq)];
            
            for i = 1: N
                [U, ~, ~, ~] = disloc3d(mdl(:, i), obs, obj.mu,obj.nu);
                Hw(:, i)      = reshape(U, 3*Nq, 1);
            end
            
        end  
        
        function Kb = get_Kb(obj, w0)
            % compute the bulk stiffness of the crack, ignore the
            % compressibility of the fluid.
            %
            % w0: the crack width
            % assume uniform pressure within the crack and caculate the
            
            N       =  numel(obj.x); % number of fault elements.
            p        =  ones(N, 1); % unit uniform pressure.
            w       =  obj.K_inv * p; % open dislocation.
            dw     = mean(w);
            s        = dw/w0; % strain.
            Kb      = 1/s; % bulk stiffness.
        end
        
        function [K, K_inv] = ddm_matrix_fullspace(obj)
            % implement the elastic kernel in 3D elastic fullspace.
            
        end
        
        function U = eval_disp_p(obj, p, xq, yq)
            % given solution of pressure distribution on the fault and query points location on the surface, evaluate the
            % surface displacement U = (Ux, Uy, Uz), [east, north, up].  disloc3d will sum up the contribution from each fracture
            % element for each query points.
            % This evaluation part is efficient.
            
            % U: [3, nq]
            
            % p, pressure at the center of each fracture cell.
            N = numel(p);
            Nq = numel(xq);
            
            if size(p,2) ~= 1
                p = p(:);
            end
            
            w = (obj.K_inv*p)'; % opening.
            
            % the dimension for each fracture cell.
            nx = size(obj.X, 2) - 1;
            ny = size(obj.X, 1) - 1;
            
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
            % This evaluation part is efficient.
            %
            
            N = numel(w);
            Nq = numel(xq);
            
            if size(w,1) ~= 1
                w = w(:)';
            end

            % the dimension for each fracture cell.
            nx = size(obj.X, 2) - 1;
            ny = size(obj.X, 1) - 1;
            
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
