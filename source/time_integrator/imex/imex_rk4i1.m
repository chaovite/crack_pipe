function y = imex_rk4i1(y0,t,dt,f,g)
% y = imex_rk4i1(y0,t,dt,f,g)
% computes y' = f(t,y(t)) + g(t,y(t)) using 4th-order runge-kutta and implicit
% euler.
%
% Input arguments:
% y0    : Vector of initial conditions.
% t     : Current time.
% dt    : Time step.
% f     : Function f(y,t) stiff part (solved implicitly).
% g     : Function g(y,t) nonstiff part (solved explicitly).
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
    2277821191437/14882151754819];
c = [0
    1432997174477/9575080441755
    2526269341429/6820363962896
    2006345519317/3224310063776
    2802321613138/2924317926251];

% Set initial condition
y      = y0;



% Iterate in time
dy = 0;
I = eye(length(y));
for k=1:length(a)

    % Set rates
    dy = a(k)*dy + g(y,t+c(k)*dt); 

    % Update fields
    y = y + dt*b(k)*dy;

    y = fsolve(@(z) z - y - dt*f(z,t+c(k)*dt),y);

end   




end




