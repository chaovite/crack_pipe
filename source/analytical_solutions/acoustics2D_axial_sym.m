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
% g= int_{-inf}^{+inf} ghat(omega)*exp(i*omega*t) domega
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



