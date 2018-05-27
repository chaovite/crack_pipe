% OkadaGrid example
addpath(genpath('../source'))

% mid point of top edge.
x = 0; % east
y = 0; % north
z = 500; % depth
Lx = 1000; % crack dimension along dip.
Ly = 1000; % crack dimension along strike.
dip =  -30; % follow the notation of disloc3d to control the top edge
strike = 30;
mu  = 3e10; % shear modulus
nu   = 0.25; % poisson ratio.
xq   = [0, 500]';   % east, observation points. 
yq   = [1000, 0]'; % north, observation points.
N    = 16;  % number of segments in x and y direction. N*N fault elements.
plot_geometry = true;

tic
 [Gtrans, Gtilt, K, K_inv] = okadaGrid(x, y, z, Lx, Ly, strike, dip, mu, nu, xq, yq, N, plot_geometry);
toc
 