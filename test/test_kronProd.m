% test_kronecker product efficiency improvement

% build model.
source_dir = '../source';
addpath(genpath(source_dir));

% conduit parameters
Mc.R = 10;
Mc.L = 1000;
Mc.nz = 80;
Mc.nr = 30;
Mc.order = 6;
Mc.S = pi*Mc.R^2;
Mc.g = 10;
Mc.mu = 50;
Mc.interface_split = true;
z = Mc.L/Mc.nz*[0:Mc.nz]';
Mc.tau = 1;
[Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);

if Mc.interface_split
    Mc.split_index = find(Mc.ex, 1) - 1;
    z = Mc.L/Mc.nz*[0:(Mc.split_index -1),Mc.split_index, Mc.split_index:Mc.nz]';
    [Mc.rho, Mc.K, Mc.c, Mc.a, Mc.b, Mc.p0, Mc.ex] = magma_st(z);
end

Mc.pT.A = 5e3; % pressure perturbation amplitude
Mc.pT.T = 1; % pressure perturbation duration
Mc.pT.t = 2; % pressure perturbation center time
Mc.G = @(t) Mc.pT.A*exp(-0.5*((t-Mc.pT.t)/Mc.pT.T)^2); % external force from the conduit.
% Mc.G = @(t) 0; % external force from the conduit.

% fluid-filled fracture parameters.
Mf.w0 = 2;
Mf.Lx   = 2e3;
Mf.Ly   = 2e3;
Mf.nx = 50;
Mf.ny = 50;
Mf.nz = 30;
Mf.order = 6;
Mf.interp_order = 4;
Mf.xs = 0.75*Mf.Lx;
Mf.ys = 0.5*Mf.Ly;
Mf.G  = @(t) 0;
Mf.xc = 0.75*Mf.Lx;% the coupling location to the conduit.
Mf.yc = 0.5*Mf.Ly;% the coupling location to the conduit.

Mf.isrigid = false;
Mf.r_g  =  0.3;% ratio of grid points in boundary layer.
Mf.r_bl =  0.15; % estimated ration of boundary layer.

% fluid and solid properties.
Mf.rho = Mc.rho(1);
Mf.c    = Mc.c(1);
Mf.K    = Mf.rho*Mf.c^2;
Mf.mu = 20;

Mf.cp    = 5e3;
Mf.cs    = 2.7e3;
Mf.rhos = 3e3;
Mf.Gs   = Mf.cs^2*Mf.rhos;
Mf.nu   = ((Mf.cp/Mf.cs)^2-2)/((Mf.cp/Mf.cs)^2-1)/2;
Mf.Xc   = 0;
Mf.Yc   = 0;
Mf.Zc   = 1000 ;
Mf.strike = 0;
Mf.dip = 0;

x_obs = 0;
y_obs = 1000;
%% construct coupled models.
Model = coupledModel_NoMatrix(Mc, Mf);
%% test the performance improvement.
Dx2  =  Model.frac.op.ux.Dx2;
Wz3 =  Model.frac.op.vx.Wz3;
vx     = ones(size(Wz3,2), 1); 

% Dx2 = kron(Dxp, Iyp).
Iyp =  speye(Model.frac.geom.ux.ny);
Dxp = Model.frac.op.ux.Dx1;

% Wz3  = kron(kron(Ixm, Iyp), Wz1);
Wz1  = Model.frac.op.vx.Wz1;
Ixm   = speye(Model.frac.geom.vx.nx);

% METHOD I: based on kron()
Kt = sparse(ones(size(Dx2,1)));
%
disp('matlab kron()')
tic
for i =1:1000
    tmp1 = Kt*(Dx2* (Wz3*vx) + Dx2* (Wz3*vx));
end
toc;

%
% METHOD II: based on KronProd object.

% Dx2 = kron(Dxp, Iyp).
% Wz3  = kron(kron(Ixm, Iyp), Wz1) = kron(Ixm, kron(Iyp, Wz1));
Wz3_Kron = KronProd({ Wz1, Iyp, Ixm},[1,2,3]);
Dx2_Kron  = KronProd({Iyp, Dxp},[1 2]);

disp('KronProd class')
tic;
for i =1:1000
    tmp2 = Kt*(Dx2_Kron*(Wz3_Kron*vx));
end
toc;

% KronProd class saves a 2% of run time but dramatically reduce the memory.
% However, we cannot insert the operators directly into sparse matrix,
% which is kinda of a pain.



















































