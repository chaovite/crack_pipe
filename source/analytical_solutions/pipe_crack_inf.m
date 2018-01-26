function [p, v, t, R, omega] = pipe_crack_inf(L, a, d, rho, cp0, cc0, ...
                                                                          g, dt, z, r, mu)
% [p, v, t, R, omega] = pipe_crack_inf(L, a, d, rho, c, ...
%                                                         g, dt, z, r, mu)
% calculate the pressure and velocity response of a pipe-crack system.
% Infinite pipe starts from [-L, +inf) and penetrates a crack at z = 0.
% The crack is radially symmetric and intersects the pipe at r = 0. The
% fluid is inviscid and crack wall is rigid (non-dispersive)
%
% Solutions based on formula in:
%  
%  Hornby et al., 1989, fracture evaluation using reflected Stoneley-wave
%  arrivals.
%
%  boundary condition, p(x = -L, t) = g(t).
%
% Notation of Fourier transform: 
%   ghat = int_{-inf}^{+inf} g(t)*exp(i*omega*t) dt, consistent with Hornby et al. 1989.
%
% argins:
%   L: distance from pipe top to crack interface.
%   a: pipe radius
%   d: crack total opening.
%   rho: fluid density
%   c: fluid wave speed.
%   g: pressure source time function at the top
%   dt: time step of the source.
%   z: z coordinate of the query point, range: [-L, +inf)
%   r: r coordinate from the crack center. range [a, +inf)
%  mu: viscosity
% 
% argouts:
%   p: pressure
%   v: velocity in +z or +r direction depending on the location of the query
%       point.
%   t : time vector
%   R: reflection coefficient.
%   omega: angular frequency.
%

if nargin<11
    mu = 0;
end

if z>0
    error('z must be [-L, 0]');
end

% location of query point.
if z ==0 && r >= a
    loc_q = 'crack';
else
    loc_q = 'pipe';
end

[ghat, f] = fft_dim(g, dt);
omega = f*2*pi;

alpha  = 8*mu/a^2; % drag coefficient in the pipe.
beta    = 12*mu/d^2; % drag coefficient in the crack.

if d==0
    beta = 0;
end

cp = cp0*sqrt(1./(1+1i*alpha./(rho*omega))); % phase velocity in pipe.
cc = cc0*sqrt(1./(1+1i*beta./(rho*omega))); % phase velocity in crack.

% k   = omega./c;
kp = omega./cp;
kc = omega./cc;

%%  reflection and transmission coefficient:
% solution in the pipe:
% p = [exp(i*k*z) +  R*exp(-i*k*z)]*A, when z<0
% v  = c/(rho*c0^2)*[exp(i*k*z) -  R*exp(-i*k*z)]*A, when z<0
% where A is to be determined by b.c.
%
% solution in the crack:
% p(r) = GH01(k*r), r>a
% v(r) = i*G*H11(k*a)c/(rho*c0^2), r>a 
%

% hankel's function.
h01 = besselh(0,1,kc*a);
h11 = besselh(1,1,kc*a);
% reflection and transmission:
F      = (2*1i*d*cc.*h11) ./(a*h01.*cp)*cp0^2/cc0^2;
R     =  (1 - F)./(1 + F);
A     = ghat./(exp(1i*kp*(-L)) + R.*exp(-1i*kp*(-L)));
G     = A.*(1+R)./h01;

switch loc_q
   case 'pipe'
        phat =  (exp(1i*kp*z) +  R.*exp(-1i*kp*z)).*A;
        vhat =  (exp(1i*kp*z)  -  R.*exp(-1i*kp*z)).*A.*cp/(rho*cp0^2);
    case 'crack'
        phat = G.*besselh(0,1,r*kc);
        vhat = 1i*G.*besselh(1,1,kc*r).*cc./(rho*cc0^2);
end

df = f(2) - f(1);
phat(1) = 0; % remove dc component.
vhat(1)  = 0; % remove dc component.
[p, t] = ifft_dim(phat, df);
[v, ~] = ifft_dim(vhat, df);
end

function [G,f] = fft_dim(g,dt)
% Notation of Fourier transform: 
%   ghat = int_{-inf}^{+inf} g(t)*exp(i*omega*t) dt  

% Fourier transform real g(t) to G(f), where f=frequency
% (f is NOT natural frequency, omega)
  
  N = length(g);
  % drop the last element if the signal is odd.
   if mod(N,2)~=0, g=g(1:N-1); N=N-1; end
   
  fN = 1/(2*dt); % Nyquist
  f  = fN*linspace(0,1,N/2+1);

  G = dt*conj(fft(g));
  G = G(1:N/2+1);
end

function [g,t] = ifft_dim(G,df)
% Inverse Fourier transform: 
% g= int_{-inf}^{+inf} ghat(omega)*exp(-i*omega*t) domega
% 
% check if the G is a column vector or a row vector

  if size(G,1) == 1
    g = [G conj(G(end-1:-1:2))];%row vector
  else
    g = [G; conj(G(end-1:-1:2))];%column vector
  end

  N = length(g);
  
  t = [0:N-1]/(N*df);
  g = ifft(conj(g))/t(2);

end
