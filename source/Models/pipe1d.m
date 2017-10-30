classdef pipe1d
    % Implement 1d pipe model, acoustic-gravity wave, velocity is cross-sectional averaged.
    % assume Poiseuille flow, du/dt + dp/dz = -alpha*u, alpha = 8*mu/R^2.
    % Don't consider non-equilibrium gas exsolution.
    % The pipe can have multiple interfaces, either 'm' (material
    % properties) or 'c' (coupled to a crack)
   
    properties
    M;              % material properties, including conduit geometry and grid configuration.
    op;             % discretization and width averaging operators
    grd;            % grid in z.
    u;               % unknowns in the conduit, uz, pz, h.
    indu;           % field index
    Ae;             % the part of A that treated explicitly PDE part.
    Ai;              % the part of A that treated implicitly.
    Fp;
    E;               % energy norm
    
    % material properties in M:
    % nr, R, nz, L, order. S (conduit cross-section area), 
    % rho, c, K, a, b.
    end
    
    methods
        function obj = pipe1d(M)
            [check, M] = pipe1d_precheck(M);
            if ~check
                return
            else
                obj.M = M;
            end
            % discretize:
            [obj.grd, obj.op] = grids_pipe1d(M);
            obj = obj.field_index();
            [obj.Ai, obj.Ae, obj.Fp] = obj.build();
            obj = obj.init();
        end
        
    function obj = init(obj)
        dim = obj.dimensions();
        obj.u = zeros(sum(dim), 1);
    end

   function obj = update(obj, u1)
       obj.u = u1;
   end
    
    function dim = dimensions(obj)
        % dimension for vz, pz, h.
      if strcmp(obj.M.bc.tp.type, 'h')
          % moving top surface equation for h.
          % Only this one is implemented.
          dim = [obj.grd.nz, obj.grd.nz, 1];
      else
          dim = [obj.grd.nz, obj.grd.nz];
      end
    end

    function  [vz, pz, h] = fields(obj, u)
%     dimension for vz, pz, h.
        indv  = obj.indu.indv;
        indp  = obj.indu.indp;
        indh  = obj.indu.indh;
        vz     = u(indv);
        pz     = u(indp);
        h       = u(indh);
    end
     
    function  obj = field_index(obj)
%     dimension for vz, pz, h.
        dim   = obj.dimensions();
        obj.indu.indv  = [1:dim(1)]';
        obj.indu.indp  = [(dim(1)+1): sum(dim(1:2))]';
        obj.indu.indh  = sum(dim);
    end
     
    function [Ai, Ae, Fp] = build(obj)
        % build part of the matrix that needs to be integrated implicitly.
        dim  = obj.dimensions();
        nz    = obj.grd.nz; % grid points for 
        rho    = obj.grd.rho;
        rhoi   = inv(rho);
        K       = obj.grd.K; 
        Ki      = inv(K);
        Dz     = obj.op.Dz;
        g        = obj.M.g;
        S       = obj.grd.S;
        Si      = inv(S);
        alpha = 8*obj.grd.mu*inv(obj.grd.R.^2);
        Ai11 = -alpha*rhoi;
        Ai = block_matrix(dim, dim, 0);
        Ai = block_matrix_insert(Ai, dim, dim, 1, 1, Ai11);
        
        % Ae is the PDE part of the operators.
        Ae11 = sparse(dim(1), dim(1));
        Ae12 = - rhoi*Dz - g*Ki;
        Ae13 = sparse(dim(1), dim(3));
        Ae21 = - K*Dz + rho*g;
        Ae22 = sparse(dim(2), dim(2));
        Ae23 = sparse(dim(2), dim(3));
        Ae31 = sparse(dim(3), dim(1));
        Ae32 = sparse(dim(3), dim(2));
        Ae33 = sparse(dim(3), dim(3));
        
        % build interface material interface conditions into Ae.
         for i = 1: length(obj.op.interfaces)
                face = obj.op.interfaces{i};
                type = face.type;
                % plus (1) and minus (2) index.
                e1    =  face.ep; %plus
                e2   =  face.em; % minus
                indp = face.indp;
                indm = face.indm;
                Z1    = obj.grd.rho(indp, indp) * obj.grd.c(indp, indp);% plus impedance
                Z2   = obj.grd.rho(indm, indm) * obj.grd.c(indm, indm);% minus impedance
                SAT1= face.SATp;
                SAT2= face.SATm;
                switch type 
                    case 'm'
                        % add terms related to material interfaces.
                        % these terms don't involve any coupling to the cracks.
                        Ae11 = Ae11- SAT1*Z2/(Z1+Z2)*e1*(e1' - e2')  ...
                                    - SAT2* Z1/(Z1+Z2)*e2*(e1' - e2');

                        Ae12 = Ae12+ 1/(Z1+Z2)*(SAT1*e1*(e2' - e1') + SAT2*e2*(e2' - e1'));
                        Ae22 = Ae22  - SAT1*Z1/(Z1+Z2)*(e1*(e1'-e2')) ...
                                    + SAT2*Z2/(Z1+Z2)*(e2*(e1'-e2'));
                        Ae21 = Ae21 - Z1*Z2/(Z1+Z2)*(SAT1*e1*(e1' - e2') + SAT2*e2*(e1' - e2'));
                    case 'c'
                        % add terms contributed by from pipe unknowns.
                        Ae12 = Ae12 + 1/(Z1*Z2)*(Z1*SAT2*(e2*e2') - Z2*SAT1*(e1*e1'));
                        Ae22 = Ae22 - SAT1*(e1*e1') - SAT2*(e2*e2');
                    otherwise
                end
                        
         end
         
         % Note that, if the interface is coupled to the crack, additional
         % terms need to added, but this requires knowing the information
         % of the crack, which will be added when building a coupled model.

         % check if the bottom is p=0 or v=0 boundary condition. if so, add
         % these terms into the Ae.
         bt = obj.op.bc.bt;
         type = bt.type;
         SAT = bt.SAT;
         e = bt.e;
         Z    = obj.grd.rho(1, 1) * obj.grd.c(1, 1);
         switch type
             case 'v'
                 % v0 = 0
                 Ae11 = Ae11 - SAT*(e*e');
                 Ae21 = Ae21 - SAT*Z*(e*e');
             case 'p'
                 % p0 = 0
                 Ae22 = Ae22 - SAT*(e*e');
                 Ae12 = Ae12 - SAT/Z*(e*e');
             case 'c'
                 % for crack bottom, build in the contribution from pipe
                 % pressures.
                 Ae22 = Ae22 - SAT*(e*e');
                 Ae12 = Ae12 - SAT/Z*(e*e');
         end
         
         % compute the source term vector from surface excitation.
         tp = obj.op.bc.tp;
         SAT = tp.SAT;
         e      = tp.e;
         Z      = obj.grd.rho(end, end) * obj.grd.c(end, end);
         rho_g = obj.grd.rho(end, end)*obj.M.g;
         
         % top moving surface.
         Ae12 = Ae12 + SAT/Z*(e*e');
         Ae13 = Ae13 - SAT/Z*rho_g*e;
         Fp1    = -SAT/Z*e;
         Ae22 = Ae22 - SAT*(e*e');
         Ae23 = Ae23 + SAT*rho_g*e;
         Fp2    = SAT*e;
         Ae31  = Ae31 + e';
         Ae32  = Ae32 + e'/Z;
         Ae33  = Ae33 - rho_g/Z;
         Fp3    = 1/Z;
         Ae = [Ae11, Ae12, Ae13; Ae21, Ae22, Ae23; Ae31, Ae32, Ae33];
         Fp = [Fp1; Fp2; Fp3];
    end
    
        function [cmax, hmin] = getCFL(obj)
            % get min(c) and max(dz) for CFL condition.
            cmax = max(obj.M.c);
            hmin  = min([obj.grd.dz]);
        end
    
    function [du]= fun_integrate(obj, u, t)
        % a function to incorporate update from interfaces and boundary condition.
        % crack_pressure a structure that contains crack pressure from each coupling cracks.
        % For example, 
        %           crack_pressure.interface_1 = pc means the crack
        %           pressure at interface_1 is pc.
        %
        % du is the explicit update.
        % Qc is the volumetric flow rate from the conduit into the crack, a structure.
        % For example:      
        %               Qc.interface_1 = qc means rate of qc leaks out of
        %               pipe into a crack through interface_1.
        
        % update from PDE
        du    =  obj.Ae*u + obj.Fp*obj.op.bc.tp.f(t);
    end
    end
    
end

