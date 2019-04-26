% test the frac2d.
% coupled conduit and bottom crack wave simulation
clear
source_dir = '../source';
addpath(genpath(source_dir));
%%
% fluid-filled fracture parameters.
Mf.w0 = 1; 
Mf.L   = 1000; 
Mf.Lz = 1000;
Mf.nx = 20; 
Mf.ny = 20; 
Mf.order = 2; 
Mf.interp_order = 4;
Mf.xs = 0; 
Mf.G  = @(t) 0; % no external forcing in the crack. 
Mf.xc = 0.5*Mf.L;% the coupling location to the conduit.

Mf.lc  = 200;%coupling length scale.
Mf.S = Mf.lc * Mf.w0;  
Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho = 1.8e3; 
Mf.c    = 2.5e3; 
Mf.K    = Mf.rho*Mf.c^2; 
Mf.mu = 50;
Mf.cp    = 5e3; 
Mf.cs    = 2.7e3; 
Mf.rhos = 3e3; 
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;

%
Mf.Xc = 0;
Mf.Yc  = 0;
Mf.Zc =1000;
Mf.strike =0;
Mf.dip =0;

frac = frac2d(Mf);

% construct the original fluid object.
metrics = struct();
f = fluid(frac.geom,metrics,frac.op,frac.M.w0, frac.M.rho, frac.M.K, frac.M.mu, frac.M.Gs, frac.M.nu, frac.M.isrigid);
e = f.source(frac.M.xs, frac.M.interp_order);
[f.Ae, f.Ai] = f.interior();

% compare results from frac2d with the originally fluid object. Operators
% Ai, Ae, Fp should match without coupling terms with the conduit.

assert(norm(full(frac.Ai - f.Ai), 2) < 1e-10,'Fail the test of Ai');
disp('Pass test for Ai');
assert(norm(full(frac.Ae - f.Ae), 2) < 1e-10,'Fail the test of Ae');
disp('Pass test for Ae');
assert(norm(full(frac.Fp - e), 2) < 1e-10,'Fail the test of Fp');
disp('Pass test for Fp');

% test energy stability:
disp('Test stability and energy stability');
[is_stable, is_energy_stable,eig_s,eig_es] = test_energy_stability(frac.A, frac.E);
is_stable
is_energy_stable













