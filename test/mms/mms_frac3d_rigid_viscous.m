% frac3d, manufactured solution with viscous fluid and rigid wall. 
% This solution leave the acoustic part out and only test the viscous
% damping part.
%
% The acoustic part is test using mms in mms_invicid_rigid_acoustics.m

clear
syms t x y z mu k;
vx = sin(k*x)*cos(k*y)*sin(k*z)*exp(-mu*k^2*t);
vy = -cos(k*x)*sin(k*y)*sin(k*z)*exp(-mu*k^2*t);
p  = sym(0);
s = sym(0); % the source term.

ux = int(vx, z, 0, 1);%(1-cos(k))/k*sin(k*x)*cos(k*y)*exp(-mu*k^2*t)
uy = int(vy, z, 0, 1);%-(1-cos(k))/k*cos(k*x)*sin(k*y)*exp(-mu*k^2*t)

disp('frac3d x y z directions:')
res_vx = simplify(diff(vx, t) + diff(p, x) - mu*diff(vx,z,2)) % x momentum
res_vy = simplify(diff(vy, t) + diff(p, y) - mu*diff(vy,z,2)) % y momentum
res_p  = simplify(diff(p,t) + diff(ux, x) + diff(uy, y) - s)     % mass balance

%% xz direction
vx = sin(k*x)*sin(k*z)*exp(-mu*k^2*t);
vy = sym(0);
p  = sym(0);
s = (1-cos(k))*cos(k*x)*exp(-mu*k^2*t); % the source term.

ux = int(vx, z, 0, 1);%(1-cos(k))/k*sin(k*x)*exp(-mu*k^2*t)
uy = int(vy, z, 0, 1);

disp('frac3d x z directions:')
res_vx = simplify(diff(vx, t) + diff(p, x) - mu*diff(vx,z,2)) % x momentum
res_vy = simplify(diff(vy, t) + diff(p, y) - mu*diff(vy,z,2)) % y momentum
res_p  = simplify(diff(p,t) + diff(ux, x) + diff(uy, y) - s)     % mass balance


%% yz direction
vx = sym(0);
vy = sin(k*y)*sin(k*z)*exp(-mu*k^2*t);
p  = sym(0);
s = (1-cos(k))*cos(k*y)*exp(-mu*k^2*t); % the source term.

ux = int(vx, z, 0, 1);
uy = int(vy, z, 0, 1);%(1-cos(k))/k*sin(k*y)*exp(-mu*k^2*t)

disp('frac3d y z directions:')
res_vx = simplify(diff(vx, t) + diff(p, x) - mu*diff(vx,z,2)) % x momentum
res_vy = simplify(diff(vy, t) + diff(p, y) - mu*diff(vy,z,2)) % y momentum
res_p  = simplify(diff(p,t) + diff(ux, x) + diff(uy, y) - s)     % mass balance
