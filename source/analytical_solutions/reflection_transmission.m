function [r, t] = reflection_transmission(f, a, d, c)
% reflection and transmission coefficient of tube wave at a infinite
% inviscid crack.
%
%   f: frequencies.
%   a: pipe radius
%   d: crack total opening.
%   c: fluid wave speed.
%        c is constant for non-dispersive crack
%        c is function of frequency for dispersive crack.
%
% Solutions based on formula in:
%  
%  Hornby et al., 1989, fracture evaluation using reflected Stoneley-wave
%  arrivals.
%
% Notation of Fourier transform: 
%   ghat = int_{-inf}^{+inf} g(t)*exp(i*omega*t) dt, consistent with Hornby et al. 1989.
%

k = f*2*pi/c; % wave number.

% hankel's function
h01 = besselh(0,1,k*a);
h11 = besselh(1,1,k*a);

% reflection and transmission:
r     =  - 1i*d*h11./(a*h01)./(1 + 1i*d*h11./(a*h01));
t     =  1./(1 + 1i*d*h11./(a*h01));
end

