function [rho, K, c, a, b, p, n]=magma_st(x, p_ref)
%compute the depth dependent magma static properties given depth x 
%using ode45. All the parameters related to the equation 
%of state of magma is built into magma_eos.m.
%Input:
%       x:  the grid points (Note that x=0 means surface and increases as depth increases)
%
%Output:
%       n = gas exsolution mass fraction
%       rho = density
%       c = sound speed
%       K = bulk modulus
%       a = gas exsolution, -(drho/dn)/rho
%       b = gas exsolution, -dneq/dp
%
%Note: to be consistent with the sign notation of the rest of the code.
%first entry means the bottom while the last entry means the surface.
%Refer magma_eos.m for more parameters related to the eos of magma.

% Integrate dp/dz=rho*g with rho=rho(p,n_eq);

%reference pressure-atmosphere pressure in this case
if nargin < 2
     p_ref=1e5;
end

%Numeric integration using ode45.
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(@f,[x(end) x(1)],p_ref,options);
p=deval(sol,x)';

[rho,~,~, K, c, a, b,n]=magma_eos(p);

end

% define function dp/dx.
    function f=f(x, p)
        g=10;
        rho=magma_eos(p);
        f=-rho*g;
    end

