% test_conduit_internal_g
%
% test the new conduit model reduces to the Karlstrom and Dunham code when
% taking Mg = 0.
%
clear

sourcedir = '../source';
addpath(genpath(sourcedir));

Mc.R = 20; 
Mc.L = 1000; 
Mc.nz = 20; 
Mc.nr = 30; 
Mc.order = 6;
Mc.order_r = Mc.order;
Mc.S = pi*Mc.R^2*ones(Mc.nz+1, 1); 
Mc.g = 10;  
Mc.with_exsolution=false;
Mc.interface_split=false;
Mc.mu = 50*ones(Mc.nz+1, 1);
z = Mc.L/Mc.nz*[0:Mc.nz]';

% density profile
rho0 = 1000;
rho1 = 1500;
alpha = (log(rho1)-log(rho0))/Mc.L;
Mc.rho = rho0*exp(alpha*(Mc.L-z));
Mc.K    = Mc.rho*Mc.g/alpha;
Mc.Mg = alpha - Mc.rho*Mc.g./Mc.K;
Mc.c    = sqrt(Mc.K./Mc.rho);
Mc.epsilon = 1;
Mc.pT.A = 10e3; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
Mc.bottom_bc = 'p=0';
cond = conduit(Mc);
cond_g = conduit_internal_g(Mc);

%% check all the matrices agree when Mg = 0. 
dims = cond_g.dimensions;
indx = [1: sum(dims(1:2)), sum(dims)];
fields = {'A','Ai','Ae','Fp','E'};
for i = 1:length(fields)
    name = fields{i};
    prop = cond.(fields{i});
    prop_g = cond_g.(fields{i});
    switch name
        case {'A','Ai','Ae','E'}
            assert(norm(full(prop-prop_g(indx,indx)))<1e-16);
        case 'Fp'
            assert(norm(full(prop-prop_g(indx)))<1e-16);
    end
    disp([name, ' PASS!'])
end



