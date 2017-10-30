function y = imex_rk4i1_lu(y0,t,dt,g,A,L,U,p,q)
% y = imexrk4i1_lu(y0,t,dt,f,g)
% computes y' = f(t,y(t)) + A*y using low storage 4th-order runge-kutta and implicit
% Euler, 
% where g(y(t),t) is non-stiff and A*y is stiff.
% y0 : initial condition y0=y(t0)
% dt : time step
% g  : g(y,t) explicit part
% A  : A (matrix) implicit part
% L,U : LU-decompostion of A-a*h*I (see imex_rk4i1_get_lu)     
%
% Output arguments:
% y  : solution at time t + dt.

a = [0
     -567301805773/1357537059087
     -2404267990393/2016746695238 
     -3550918686646/2091501179385
     -1275806237668/842570457699];
b = [1432997174477/9575080441755
    5161836677717/13612068292357
    1720146321549/2090206949498
    3134564353537/4481467310338
    2277821191437/14882151754819 
    1];
c = [0
    1432997174477/9575080441755
    2526269341429/6820363962896
    2006345519317/3224310063776
    2802321613138/2924317926251 
    1];

% Set initial condition
y      = y0;



% Iterate in time
dy = 0;
I = eye(length(y));
B = I - dt*A;
for k=1:length(a)

    % Set rates
    dy = a(k)*dy + g(y,t+c(k)*dt); 

    % Update fields
    y = y + dt*b(k)*dy;
    B = I - dt*(c(k+1) - c(k))*A;
    %y = B\y;
    y = B\y;
    %y(q) = U\(L\(y(p))); 
    %y = fsolve(@(z) z - y - dt*f(z,t+c(k)*dt),y);

end   




end




