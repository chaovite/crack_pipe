%frac3d, manufactured solution with inviscid fluid and rigid wall.
syms t x y z mu k;
% x y direction.
vx = - 1/sqrt(2)*sin(k*x)*cos(k*y)*cos(k*sqrt(2)*t);
vy =  - 1/sqrt(2)*cos(k*x)*sin(k*y)*cos(k*sqrt(2)*t);
p  = cos(k*x)*cos(k*y)*sin(k*sqrt(2)*t);
s = sym(0); % the source term.

disp('2D acoustics x y both directions:')
res_vx = simplify(diff(vx, t) + diff(p, x) - mu*diff(vx,z,2)) % x momentum
res_vy = simplify(diff(vy, t) + diff(p, y) - mu*diff(vy,z,2)) % y momentum
res_p  = simplify(diff(p,t) + diff(vx, x) + diff(vy, y) - s)     % mass balance

%% x direction
vx = - sin(k*x)*cos(k*t);
vy = sym(0);
p  = cos(k*x)*sin(k*t);
s = sym(0); % the source term.
disp('2D acoustics x direction:')
res_vx = simplify(diff(vx, t) + diff(p, x) - mu*diff(vx,z,2)) % x momentum
res_vy = simplify(diff(vy, t) + diff(p, y) - mu*diff(vy,z,2)) % y momentum
res_p  = simplify(diff(p,t) + diff(vx, x) + diff(vy, y) - s)     % mass balance

%% y direction
vx = sym(0);
vy = - sin(k*y)*cos(k*t);
p  = cos(k*y)*sin(k*t);
s = sym(0); % the source term.
disp('2D acoustics y direction:')
res_vx = simplify(diff(vx, t) + diff(p, x) - mu*diff(vx,z,2)) % x momentum
res_vy = simplify(diff(vy, t) + diff(p, y) - mu*diff(vy,z,2)) % y momentum
res_p  = simplify(diff(p,t) + diff(vx, x) + diff(vy, y) - s)     % mass balance