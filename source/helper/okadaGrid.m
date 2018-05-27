function [Gtrans, Gtilt, K, K_inv] = okadaGrid(x, y, z, Lx, Ly, strike, dip, mu, nu, xq, yq, N, plot_geometry)
% [Gtrans, Gtilt, K, K_inv] = okadaGrid(x, y, z, Lx, Ly, strike, dip, mu, nu, xq, yq, N, plot_geometry)
%
% Input:
%
%  x: east (m), mid-point upper edge
%  y: north (m), mid-point upper edge
%  z: depth (m), mid-point upper edge
%  Lx: crack dimension (m) along the dip.
%  Ly: crack dimension (m) along the strike.
%  strike: angle from the north (degree), range from 0 to 359
%  dip:  angle from the horizontal plane (degree), range from -90 to 0,
%             consistent with upper edge mid-point notation used in disloc3d.
%  mu: solid shear modulus (Pa)
%  nu: solid Poisson's ratio
%  xq: east (m), of seismic stations, dimension is [nq, 1].
%  yq: north (m), of seismic stations, dimension is [nq, 1]
%  N: number of small crack elements in each (x or y) direction.
% plot_geometry: True or False, If the geometry of the crack should be plotted.
%
%  Output:
%
%  Gtrans: the translation green's function of unit pressure on the crack.
%  Gtilt:     the tilt green's function of the unit pressure on the crack.
%  K      : stiffness matrix, p = K*w
%  Kinv : compliance matrix, w = Kinv*p, Kinv = inv(K).
%
% Note: 
% Given Gtrans, Gtilt and angular frequency, total Green's function can be constructed:
% Gtotal = Gtrans +  Gtilt*g/(i*omega)^2
% During later inversion when omega is specified, total Green's function can be constructed.
%
% Also note that these are Green's function for displacement, a factor of
% (-i*omega) or (i*omega) is needed to obtain velocity depending on the
% notation of Fourier transform.
%


if nargin < 13
    plot_geometry = false;
end


strike = strike/180*pi;
dip     = dip/180*pi;

% both Lx, Ly are divided into N segments.
dLx = Lx/N; % along the dip
dLy = Ly/N; % along the strike.

%% obtain the east, north, depth of the top edge middle point of each crack element.
yg_e   = (0.5:1:(N-0.5))*dLy - Ly/2;
xg_e   = (0:1:(N-1))*dLx;

yg_c   = (0.5:1:(N-0.5))*dLy - Ly/2;
xg_c   = (0.5:1:(N-0.5))*dLx;


% meshgrid of middle point of top edges
[xG_e, yG_e]= meshgrid(xg_e, yg_e);
xG_e = xG_e(:);
yG_e = yG_e(:);

% meshgrid of crack element centroids.
[xG_c, yG_c]= meshgrid(xg_c, yg_c);
xG_c = xG_c(:);
yG_c = yG_c(:);

% unit vector in dip and strike direction.
vec_x  = [-cos(dip)*cos(strike), sin(strike)*cos(dip), sin(dip)];
vec_y  = [sin(strike), cos(strike), 0];
XYZ_e    = [xG_e,yG_e]*[vec_x; vec_y];
XYZ_c    = [xG_c,yG_c]*[vec_x; vec_y];

% the coordinates of the top edge middle point for each crack element.
X_e = XYZ_e(:, 1) + x; % east
Y_e = XYZ_e(:, 2) + y; % north
Z_e = -XYZ_e(:, 3) + z; % depth

% the coordinates of the centroid point for each crack element.
X_c = XYZ_c(:, 1) + x; 
Y_c = XYZ_c(:, 2) + y;
Z_c = -XYZ_c(:, 3) + z;

%% create stiffness matrix
% create observation at each crack cell center.
obs = [X_c'; Y_c'; -Z_c'];

M = N*N;

% create the stiffness matrix to relate opening to pressure.
K = zeros(M, M);

% loop over each crack element.

len     = dLy; % along the strike
width = dLx; % along the dip

% unit vector normal to the fault plane.
nvec   = [sin(strike + pi/2)*sin(dip), cos(strike+pi/2)*sin(dip), cos(dip)]';

for i = 1: M
    depth = Z_e(i);
    east   = X_e(i);
    north  = Y_e(i);
    mdl = [len, width, depth, dip*180/pi, strike*180/pi, east, north, 0, 0, 1]';
    [~, ~, S, ~] = disloc3d(mdl,obs, mu, nu);
    
    % Stress: Sxx, Sxy, Sxz, Syy, Syz, and Szz.
    K(:, i) = -(S(1,:)*nvec(1)^2 + 2*S(2,:)*nvec(1)*nvec(2) + 2*S(3,:)*nvec(1)*nvec(3) ...
                + S(4,:)*nvec(2)^2 + 2*S(5,:)*nvec(2)*nvec(3) + S(6,:)*nvec(3)^2);
end
K_inv=inv(K);

% calculate opening for uniform pressure.
ws = K_inv*ones(M, 1);

%% calculate Green's function from unit pressure.

mdls = [len*ones(1,M); width*ones(1,M);
           Z_e'; dip*ones(1, M)*180/pi; strike*ones(1,M)*180/pi; 
           X_e'; Y_e'; zeros(1,M); zeros(1,M); ws'];

Nq = length(xq);
obs = [xq(:)'; yq(:)'; zeros(1, Nq)];
[Gtrans, ~, ~, ~] = disloc3d(mdls, obs, mu, nu);

obs_xp = [xq(:)'+1; yq(:)'; zeros(1, Nq)]; % x  - 1
obs_xm = [xq(:)'-1; yq(:)'; zeros(1, Nq)]; % x + 1
obs_yp = [xq(:)'; yq(:)'+1; zeros(1, Nq)]; % y  - 1
obs_ym = [xq(:)'; yq(:)'-1; zeros(1, Nq)]; % y + 1

[Gtrans_xp, ~, ~, ~] = disloc3d(mdls, obs_xp, mu, nu);
[Gtrans_xm, ~, ~, ~] = disloc3d(mdls, obs_xm, mu, nu);
[Gtrans_yp, ~, ~, ~] = disloc3d(mdls, obs_yp, mu, nu);
[Gtrans_ym, ~, ~, ~] = disloc3d(mdls, obs_ym, mu, nu);

Gtilt = zeros(Nq*3, 1);

% compute tilt using 2nd order central difference.
% tilt in x direction.
Gtilt(1:3:end, :) =  (Gtrans_xp(1:3:end, :) - Gtrans_xm(1:3:end, :))/2;
% tilt in y direction.
Gtilt(2:3:end, :) = (Gtrans_yp(2:3:end, :) - Gtrans_ym(2:3:end, :))/2;

%% plot the fault geometry and opening distribution if enabled.
if plot_geometry
    figure(1);
    plot3(X_e, Y_e, -Z_e,'k.','markersize',10);
    hold on;
    plot3(X_c, Y_c, -Z_c,'r.','markersize',10);xlabel('east');ylabel('north');zlabel('z');
    scatter3(X_c, Y_c, -Z_c, 50, ws);
    colorbar;
    grid minor;
    xlabel('East');
    ylabel('North');
    zlabel('Z');
    set(gca,'fontsize',16);
    legend({'top edge mid point','crack cell centroid','opening'});
    title('Fault geometry');
    daspect([1,1,1])
    hold off;
    
    figure(2);
    pcolor(xg_c, yg_c,reshape(ws,N,N));
    colorbar;
    shading interp
end

end

