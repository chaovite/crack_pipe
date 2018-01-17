function [p, v, t] = acoustics2D_axial_sym(a, rho, c0, g, dt, r, d, mu, btype)
% outgoing wave for axial-symmetrical 2d acoustics. velocity boundary
% condition is prescribed at r=a. 
%

if nargin<9
    btype='p';
end

beta = 12*mu/d^2;
[ghat, f] = fft_dim(g,dt);
omega = 2*pi*f;
c  = c0./sqrt(1+1i*beta./(rho*omega));
k  = omega./c;

H01a = besselh(0,1,k*a);
H01r = besselh(0,1,k*r);
H11a = besselh(1,1,k*a);
H11r = besselh(1,1,k*r);

switch btype
    case 'p'
        G  = ghat./H01a;
    case 'v'
        G  = ghat.*rho*c0^2./(1i*c.*H11a);
end

phat  = G.*H01r;
vhat  = 1i*c/rho/c0^2.*G.*H11r;
phat(1) = 0;
vhat(1) = 0;
df = f(2) - f(1);

[p, t] = ifft_dim(phat, df);
[v, ~] = ifft_dim(vhat, df);
end