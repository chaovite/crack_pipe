function [rho, rho_g, rho_l, K, c, a, b,n_eq]=magma_eos(p,n)
% magma_eos(p,n) returns the general values whil magma_eos(p) return
% equibrium values.
%
% p: pressure
% n: gas exsolution
% rho = density of mixture
% rho_g: gas density
% rho_l: liquid density
% c = sound speed
% K = bulk modulus
% a = gas exsolution, -(drho/dn)/rho
% b = gas exsolution, -dneq/dp

% equation of state:
% a mixture of ideal gas and melt with constant compressibility.
% gas: rho_g=p/(R*T_m).
% melt:rho_l=rho_l_ref*(1+(p-p)*beta_l).

n_t =0.5e-2;% total volative content 0.1% w.t of water

R = 461.5; %[J/kg/K] specific gas constant for water vapor
T_m =1000+272.15;% Temperature of magma about 1000 deg C.

beta_l = 1e-10; %[Pa^-1] compressibility constant for magma melt
K_l = 1/beta_l;% [Pa] bulk modulus of magma melt
rho_l_ref = 2800;%[kg/m^3], reference density, typical density for basalt
p_ref = 1e5;%[Pa], reference pressure, atmosphere pressure

% Constants for Henry's law. n_eq=n_o-s*p^m;
s_h = 4e-6;%[Pa^(-1/2)]. solubility constant.
m = 0.5;% common for Rhyolite magma.
p_ex = (n_t/s_h)^(1/m);%[Pa] exsolution pressure.

if nargin<2 % if magma_eos(p) return equilibrium value
    n=n_t-s_h*p.^m;
    n(p>p_ex)=0;
    n_eq=n;
end

%compute density
rho_g = p./(R*T_m);
rho_l = rho_l_ref*(1+(p-p_ref)*beta_l);
rho = 1./(n./rho_g+(1-n)./rho_l);

K = 1./(rho.*(n./(rho_g.*p)+(1-n)./(rho_l*K_l)));
a = rho.*(1./rho_g-1./rho_l);
b = s_h*m*p.^(m-1);
% if p>p_ex, b=0; 
b(p>p_ex)=0;

%instantaneous wave speed
c=sqrt(K./rho);

end 
