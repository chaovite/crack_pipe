% test the conduit
% make sure the matrix Ai, Ae and A of the new code are the same with the
% old code.
source_dir = '../source';
addpath(genpath(source_dir));
%%
% conduit parameters.
% conduit parameters
Mc.R = 20; 
Mc.L = 1000; 
Mc.nz = 50; 
Mc.nr = 50; 
Mc.order = 6;
Mc.order_r = 6;
Mc.S = pi*Mc.R^2; 
Mc.g = 10;  
Mc.mu = 50; 
z = Mc.L/Mc.nz*[0:Mc.nz]'; 
Mc.tau = 1;
Mc.interface_split=false;
[Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex]=magma_st(z);
Mc.pT.A = 10e3; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
Mc.bottom_bc = 'p=0';

cond = conduit(Mc);

[A,G,F,~,~,~,~, Ai, Ae] = magma_2d(Mc.nr,Mc.R,Mc.nz,Mc.L,Mc, Mc.order);

 assert(norm(full(cond.Ai - Ai), 2) < 1e-10,'Fail the test of Ai');
 disp('Pass test for Ai');
 assert(norm(full(cond.Ae - Ae), 2) < 1e-10,'Fail the test of Ae');
 disp('Pass test for Ae');
 assert(norm(full(cond.Fp - F), 2) < 1e-10,'Fail the test of Fp');
 disp('Pass test for Fp');
