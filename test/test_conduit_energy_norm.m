% test conduit energy norm

source_dir = '../source';
addpath(genpath(source_dir));
%%
% conduit parameters.

Mc.R = 20; 
Mc.L = 1000; 
Mc.nz = 30; 
Mc.nr = 20; 
Mc.order = 6;
Mc.S = pi*Mc.R^2; 
Mc.g = 10;  
Mc.mu = 50; 
z = Mc.L/Mc.nz*[0:Mc.nz]'; 
Mc.tau = 1;
[Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex]=magma_st(z);

Mc.interface_split = false;
if Mc.interface_split 
    Mc.split_index = find(Mc.ex, 1) - 1;
    z = Mc.L/Mc.nz*[0:(Mc.split_index -1),Mc.split_index, Mc.split_index:Mc.nz]';
    [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
end

Mc.pT.A = 10e3; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
Mc.bottom_bc = 'p=0';

cond = conduit(Mc);
[A,G,F,~,~,~,~, Ai, Ae] = magma_2d(Mc.nr,Mc.R,Mc.nz,Mc.L,Mc, Mc.order);

% test the energy stability 
[is_stable, is_energy_stable,eig_s,eig_es] = test_energy_stability(cond.A, cond.E);
