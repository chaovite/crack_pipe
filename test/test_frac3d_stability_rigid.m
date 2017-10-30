%test frac3d stability, rigid.
% This test pass when the units are small.
source_dir = '../source';
addpath(genpath(source_dir));
%%
% fluid-filled fracture parameters.
Mf.w0 = 2; 
Mf.Lx   = 1; 
Mf.Ly   = 1; 
Mf.nx = 10; 
Mf.ny = 10; 
Mf.nz = 10; 
Mf.order = 2; 
Mf.interp_order = 2;
Mf.xs = 0.5; 
Mf.ys = 0.5; 
Mf.G  = @(t) 0; % no external forcing in the crack. 
Mf.xc = 0.5*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

% Mf.lc  = 2*pi*Mc.R;%coupling length scale.
% Mf.S = Mf.lc * Mf.w0;  
Mf.isrigid = true;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.3; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho = 1.8; 
Mf.c    = 2.5; 
Mf.K    = Mf.rho*Mf.c^2; 
Mf.mu = 50e-3;
Mf.cp   = 5e3; 
Mf.cs    = 2.7e3; 
Mf.rhos = 3e3; 
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu  = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
%
frac = frac3d(Mf);
%%
disp('Test stability and energy stability');
[is_stable, is_energy_stable,eig_s,eig_es] = test_energy_stability(frac.A, frac.E);
is_stable
is_energy_stable








