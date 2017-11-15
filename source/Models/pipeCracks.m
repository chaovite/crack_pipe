classdef pipeCracks
  % a coupled pipe and cracks model.
  % names of cracks match the names of interfaces
  % This class assemble the information of the pipe and cracks and couple
  % them together to form a complete problem, ready for time integration.
  %
  
  properties
    pipe;           %  a pipe object.
    
    cracks;    
    %     cracks:  a structure that contains the multiple crack objects linked to the
    %                  pipe interfaces.
    %                  name of crack matches the name of interfaces defined in pipe.
    %
    
    u;                 %  unknown fields: [vz; pz; h; pxy_i; vx_i; vy_i], i is the index of crack.
    indu;            %  the index of pipe, interface_i, bt...
    idx_blocks;   % block index of each component.
    
    Ai;                %  implicit part of matrix A.
    Ae;               %  explict part of matrix A less the terms involving dense matrix multiplification.
    Ap;                % used to update pressure before multiplie by dense matrix Kt.
    ft;                 %  a function handle ft(t) = Ae*u + Fp(t), 
                        %       with both explicit terms and source time function.
                        
    Fp;               %  operators for external source from both the conduit and fracture. 
    E;                 %  energy norm 
  end
  
  methods
      %%
    function obj = pipeCracks(pipe, cracks)
      % construct a coupled model with a conduit and frac
      % build the conduit discretization
      obj.pipe     = pipe;
      obj.cracks = cracks;
      obj.precheck();
      obj.indu     = obj.fields_index();
      obj.idx_blocks = obj.blocks_id();
      
      % build the discretization of the coupled model.
      [obj.Ai, obj.Ae, obj.Ap, obj.Fp] = obj.build();
      
      % build the time integration function handle
      obj.ft = @(u, t) fun_integrate(obj, u, t);
      
      obj = obj.init();
      
%       obj.E = obj.energy_norm();
    end
    %%
    function obj = init(obj)
        dim   = obj.dimensions();
        obj.u = zeros(sum(dim), 1);
    end
    
    function precheck(obj)
        % check if the names of cracks matches crack interfaces in the pipe
        %
        cracks_pipe = obj.pipe.M.crack_interfaces;
        for i = 1: length(cracks_pipe)
            name = cracks_pipe{i}.name;
            if ~isfield(obj.cracks, name)
                error(sprintf('pipe crack interface %s is not defined in cracks', name));
            end
        end
    end

    function [Ai, Ae, Ap, Fp] = build(obj)
        % this method build the Ai, Ae, Fp from both pipe and
        % fluid-filled cracks.
        
        % Note that all the terms involving multiplication of dense matrix
        % K_t are not included in Ae and are handled separately.
        
          dim = obj.dimensions();
          idx_block = obj.idx_blocks;
          
          Ae = block_matrix(dim, dim); 
          Ai = block_matrix(dim, dim);
          Fp = block_matrix(dim, ones(1, length(dim)));
          
          blocks = fields(idx_block);
          
          for i = 1: length(blocks)
              block = blocks{i};
              idx = idx_block.(block);
              switch block
                  case 'pipe'
                      Ae = block_matrix_insert(Ae,dim,dim,idx,idx, obj.pipe.Ae);
                      Ai  = block_matrix_insert(Ai,dim,dim,idx,idx, obj.pipe.Ai);
                      Fp = block_matrix_insert(Fp,dim, ones(1, length(dim)), idx, idx, obj.pipe.Fp);
                  otherwise
                      % cracks
                      Ae = block_matrix_insert(Ae,dim,dim,idx, idx, obj.cracks.(block).Ae);

                      % add contribution from the stiff coupling term
                      % Note that there's coupling term dpxy/dt = -Kt*S/w0/Z*Hinv*(d'*d)
                      
                      % some constants to use
                        Hinv    = inv(obj.cracks.(block).op.p.Pxy2);
                        dc       = obj.cracks.(block).dc;
                        w0      = obj.cracks.(block).M.w0;
                        Kt        = obj.cracks.(block).M.K_t;
                    switch block
                        case 'bt'
                            % a crack intersect pipe at the bottom.
                             % crack to pipe
                             c      = obj.pipe.grd.c(1,1);
                             rho   = obj.pipe.grd.rho(1,1);
                             S      = obj.pipe.grd.S(1, 1);
                             Z      = rho*c; 
                             % pxy -> pxy, 
                             % This term might be unstable and should be moved to the
                             % implicit side.
                             Ai_crack11 =  -S/w0/Z*(Kt*(Hinv*dc))*dc';
                        otherwise
                            strs   =strsplit(block,'_'); % interface_1 is 1.
                            nface = str2double(strs{2});
                            % a crack intersect pipe in the middle.
                            % get the index of the interface
                            % 
                            face = obj.pipe.op.interfaces{nface};
                            indm = face.indm;
                            indp = face.indp;
                            % plus:1, minus: 2
                            Z1      = obj.pipe.grd.c(indp,indp) * obj.pipe.grd.rho(indp,indp);
                            Z2      = obj.pipe.grd.c(indm,indm) * obj.pipe.grd.rho(indm,indm);
                            % Note that here assumes the area doesn't jump.
                            S        = obj.pipe.grd.S(indp, indp);  
                            % pxy -> pxy
                            Ai_crack11 =  -S/w0*(1/Z1+1/Z2)*(Kt*(Hinv*dc))*dc';
                    end
                    
                    % add Ai_crack_11 into the matrix.
                    dim_crack = obj.cracks.(block).dimensions();
                    Ai_crack    = block_matrix(dim_crack, dim_crack);
                    Ai_crack    = block_matrix_insert(Ai_crack, dim_crack, dim_crack, 1, 1, Ai_crack11);
                    Ai  = block_matrix_insert(Ai,dim,dim,idx, idx, obj.cracks.(block).Ai + Ai_crack);
                    Fp = block_matrix_insert(Fp,dim, ones(1, length(dim)), idx, idx, obj.cracks.(block).Fp);
              end
          end
          
         % adding the coupling terms, other coupling terms are explicit.
         
         % in this case.
         Ae_c =  block_matrix(dim, dim, 0);
         [Aec2p, Ap] = cracks2pipe_couple(obj.pipe, obj.cracks, dim); %explicit from crack to pipe.
         crack_names = fields(obj.cracks);
         for i = 1: length(crack_names)
             name = crack_names{i};
             idx_c  = idx_block.(name);
             Ae_c = block_matrix_insert(Ae_c,dim, dim, 1, idx_c, Aec2p.(name));
         end
         Ae = Ae + Ae_c;
    end
    
    function [A, Fp] = getA(obj)
        % This method returns the A matrix and source vector for inversion
        % in frequency domain.
        
        % return the source term only from the excitation from the pipe.
        Fp = obj.Fp(:,1);
        A = obj.Ai + obj.Ae; % the sparse part.
        % the dense part, which gives pressure in the crack and the other
        % fields. 
        crack_names = fields(obj.cracks);
        % update crack pressure.
        for i = 1: length(crack_names)
            name = crack_names{i};
            indp   = obj.indu.(name).p;
            A(indp, :) = A(indp, :) + obj.cracks.(name).M.K_t * obj.Ap.(name);
        end
    end
    
    function ft = fun_integrate(obj, u, t)
        % defines the integration function ft(u, t), additional parameters are
        % parsed in through obj.

        % contribution less the terms involving dense matrix.
        ft = obj.Ae*u;
        s = obj.Fp(:,1)*obj.pipe.M.bc.tp.f(t); % source term from the pipe
        
        % add sources from the crack.
        crack_names = fields(obj.cracks);
        for i = 1: length(crack_names)
            name = crack_names{i};
            s = s + obj.Fp(:, i + 1)*obj.cracks.(name).M.G(t);
        end
        ft = ft + s;
        
        % update crack pressure.
        for i = 1: length(crack_names)
            name = crack_names{i};
            indp   = obj.indu.(name).p;
            qp     = obj.Ap.(name) * u;
            
            if obj.cracks.(name).M.use_fft
                kernel = obj.cracks.(name).kernel_fft;
                mask  = obj.cracks.(name).mask_fft;
                qp_fft = zeros(size(kernel));
                qp_fft(mask) = qp;
                dp      = ifft2(kernel.*fft2(qp_fft));
                ft(indp) = ft(indp) + dp(mask);
            else
                ft(indp) = ft(indp) + obj.cracks.(name).M.K_t * qp;
            end
            
        end
       
    end

%     function E = energy_norm(obj)
%       dim = obj.dimensions();
%       E    = block_matrix(dim,dim,1);
%       % TODO: Add each term in the mechanical energy balance for the coupled conduit and crack system.
%       E    = block_matrix_insert(E, dim,dim, 1,1, obj.cond.E);
%       E    = block_matrix_insert(E, dim,dim, 2,2, obj.frac.E);
%     end

        function out = getfields(obj, u, field, subfields)
                % obtain a solution
                % field is 'pipe' or the names of crack interfaces, a string.
                % subfields is a list of strings.
                %
                indfield = obj.fields_index.(field);
                ufield    = u(indfield);
                switch field
                    case 'pipe'
                        for i = 1: length(subfields)
                            subfield = subfields{i};
                            switch subfield
                                case 'vz'
                                    indsub = obj.pipe.indu.indv;
                                case 'pz'
                                    indsub = obj.pipe.indu.indp;
                                case 'h'
                                    indsub = obj.pipe.indu.indh;
                                otherwise 
                                    error('Invalid subfield for pipe');
                            end
                                out.pipe.(subfield) = ufield(indsub);
                        end

                    otherwise
                        % a crack field.
                        for i = 1: length(subfields)
                            subfield = subfields{i};
                             switch subfield
                                case 'p'
                                    indsub = obj.cracks.(field).indu.indp;
                                case 'vx'
                                    indsub = obj.cracks.(field).indu.indvx;
                                case 'vy'
                                    indsub = obj.cracks.(field).indu.indvy;
                                otherwise
                                    error('Invalid subfields for crack')
                              end
                               out.(field).(subfield) = ufield(indsub);
                        end
                end
        end
        
        function fields_index = fields_index(obj)
                % store the index of all the field.
                 dim = obj.dimensions();
                 fields_index.pipe.all = [1: dim(1)]';
                 fields_index.pipe.vz = fields_index.pipe.all(obj.pipe.indu.indv);
                 fields_index.pipe.pz = fields_index.pipe.all(obj.pipe.indu.indp);
                 fields_index.pipe.h = fields_index.pipe.all(obj.pipe.indu.indh);
                 
                 crack_names = fields(obj.cracks);
                 for i = 1: length(crack_names)
                     name = crack_names{i};
                     s = sum(dim(1: i)) + 1;
                     e = sum(dim(1: i + 1));
                     indall = [s : e]';
                     indp  = obj.cracks.(name).indu.indp;
                     indvx  = obj.cracks.(name).indu.indvx;
                     indvy  = obj.cracks.(name).indu.indvy;
                     fields_index.(name).all=indall;
                     fields_index.(name).p = indall(indp);
                     fields_index.(name).vx = indall(indvx);
                     fields_index.(name).vy = indall(indvy);
                 end
        end
    
        function obj = update(obj, U)
               obj.u = U;
               % update each components.
               obj.pipe = obj.pipe.update(U(obj.indu.pipe.all));
               names   = fields(obj.cracks);
               for i = 1: length(names)
                   name = names{i};
                   obj.cracks.(name) = obj.cracks.(name).update(U(obj.indu.(name).all));
               end
        end

    function dim = dimensions(obj)
        %[conduit, frac].
        % index
        names = fields(obj.cracks);
        nc        = length(names);
        dim      = zeros(1, nc + 1);
        dim(1) = sum(obj.pipe.dimensions());
        for i = 2: (nc + 1)
            name = names{i-1};
            dim(i) = sum(obj.cracks.(name).dimensions());
        end
    end
    
    function idx_block = blocks_id(obj)
        
          names = fields(obj.cracks);
          ncrack   = length(names);
          
          % a structure that records the block index of each components,
          % pipe and cracks.
          %
          idx_block.pipe = 1;
          for i = 1: ncrack
              name = names{i};
              idx_block.(name) = i + 1;
          end
    end
    
        function [cmax, hmin] = getCFL(obj)
            % get min(c) and max(dz) for CFL condition.
            [cmax, hmin] = obj.pipe.getCFL();
            names = fields(obj.cracks);
            for i = 1: length(names)
                name = names{i};
                [c, h] = obj.cracks.(name).getCFL();
                if c>cmax
                    cmax=c;
                end
                if h<hmin
                    hmin=h;
                end
            end
        end
        
        function U = get_disp(obj, u, xq, yq)
            % get a surface displacement given the solution and a list of
            % query points, xq, and yq.
            % the total displacement is a sum of all the contribution from
            % all the cracks.
            names = fields(obj.cracks);
            nq        = length(xq);
            U          = zeros(3, nq);
            for i = 1: length(names)
                name = names{i};
                indp   = obj.indu.(name).p;
                p        = u(indp);
                U        = U + obj.cracks.(name).eval_disp_p(p, xq, yq);
            end
        end
        
        function Hp = get_Hp(obj, xq, yq)
            % get the observation matrix Hp, which relates u to surface
            % displacement.
            Nq  = length(xq);
            dim = obj.dimensions();
            N    = sum(dim);
            Hp = zeros(3*Nq, N);
            names = fields(obj.cracks);
            
            for i = 1: length(names)
                name = names{i};
                indp   = obj.indu.(name).p;
                Hp_i   =  obj.cracks.(name).disloc3d.get_Hp(xq, yq);
                Hp(:, indp) = Hp_i;
            end
            Hp = sparse(Hp);
        end
  end

end
