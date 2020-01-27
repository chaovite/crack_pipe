% test_pennyshape_crack
clear all;
close all;
source_dir  = '../source';
mesh_dir    = '../meshes';

addpath(genpath(source_dir));
addpath(genpath(mesh_dir));

% load spherical mesh
mshfile = 'sphere.msh';
% mshfile = 'sphere_refine.mesh'; %refined mesh

tic;
[edge2el, el2edge, edgeBC, transT, edge_len, connDist, TR] = LoadGmsh2D(mshfile, true);

%% compute elasticity kernel
Ne = length(el2edge);
Np = size(TR.Points, 1);

% cell centers;
C = incenter(TR);
X = C(:, 1);
Y = C(:, 2);
Z = 0*ones(size(X));

%%
P1s = zeros(Ne, 3);
P2s = zeros(Ne, 3);
P3s = zeros(Ne, 3);

connList = TR.ConnectivityList;
P1s(:,1:2) = TR.Points(connList(:, 1), :);
P2s(:,1:2) = TR.Points(connList(:, 2), :);
P3s(:,1:2) = TR.Points(connList(:, 3), :);

Ss = 0;
Ds = 0;
Ts  = 1;
lambda = 1;
mu        = 1;
K0 = lambda + 2/3*mu;
nu = lambda/(2*(lambda+mu));

%%
K = zeros(Ne, Ne);

 for i = 1: Ne
    P1 = P1s(i, :);
    P2 = P2s(i, :); 
    P3 = P3s(i, :);
    Ts = 1;
    [Stress,~]=TDstressFS(X,Y,Z, P1,P2,P3,Ss,Ds,Ts,mu,lambda);
    K(:,i) = - Stress(:, 3);
 end
 
K_inv = inv(K);
toc

%% uniform pressure
figure;
P = ones(Ne, 1);
U = K_inv*P;
% plot the crack opening.
subplot(1,3,1);

trisurf(TR.ConnectivityList, TR.Points(:, 1),TR.Points(:, 2), zeros(Np,1),  U,'edgecolor','none');colorbar;
xlabel('X');
ylabel('Y');
title('Numerical');
colorbar;
view([0,0,1]);
daspect([1,1,1])
xlim([-0.5, 0.5]);
ylim([-0.5, 0.5])
set(gca,'fontsize',14);

subplot(1,3,2);

%
c    = 0.5;
r     =  sqrt(X.^2 + Y.^2);
Ua  = 2*4*(1-nu^2)*c*sqrt(1 - (r./c).^2)/3/pi/K0/(1-2*nu)*1;
trisurf(TR.ConnectivityList, TR.Points(:, 1),TR.Points(:, 2), zeros(Np,1), Ua,'edgecolor','none');
xlabel('X');
ylabel('Y');
title('Analytical');
suptitle('penny-shaped crack')
colorbar;
xlim([-0.5, 0.5]);
ylim([-0.5, 0.5])
view([0,0,1]);
daspect([1,1,1]);
set(gca,'fontsize',14);

subplot(1,3,3);
trisurf(TR.ConnectivityList, TR.Points(:, 1),TR.Points(:, 2), zeros(Np,1), (Ua-U)./Ua,'edgecolor','none');
colorbar;
view([0,0,1]);
xlim([-0.5, 0.5]);
ylim([-0.5, 0.5])
daspect([1,1,1]);
xlabel('X');
ylabel('Y');
title('Relative Diff');

% axis([-0.4])



