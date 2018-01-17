function [p, t] = acoustics2D_pointsource(r, c, g, dt)
% 2d linear acoustics with subject to a point source at xs, ys.
% Greens' function G = -i/4*H_0^{2}.
%
% dp+ux+uy = g(t)*delta(r).
%
% 1/c^2*d^2p/dt^2 + Lap(p) = g'(t)*delta(r)
%
% Notation of fourier transform ghat(omega) = int exp(-i*omega*t)*g(t)dt.
% this notation is different from Eric's notation.
%

[ghat, f] = fft_dim(g,dt);
omega = 2*pi*f;
k = omega/c;
ghat = (1i*omega).*ghat;
H = besselh(0, 2, k*r);
phat = -1i/4*H.*ghat;
phat(1) = 0;
[p,t] = ifft_dim(phat,f(2));

end

function [G,f] = fft_dim(g,dt)
  
% Fourier transform real g(t) to G(f), where f=frequency
% (f is NOT natural frequency, omega)
  
  N = length(g);
  % drop the last element if the signal is odd.
   if mod(N,2)~=0, g=g(1:N-1); N=N-1; end
  % zeropad the signal if the signal is odd;
   
  fN = 1/(2*dt); % Nyquist
  f  = fN*linspace(0,1,N/2+1);

  G = dt*(fft(g));
  G = (G(1:N/2+1));
end

function [g,t] = ifft_dim(G,df)

% check if the G is a column vector or a row vector
% handle the nyquist.

  if size(G,1) == 1
    g = [G conj(G(end-1:-1:2))];%row vector
  else
    g = [G; conj(G(end-1:-1:2))];%column vector
  end

  N = length(g);
  
  t = [0:N-1]/(N*df);
  
  g = ifft((g))/t(2);

end

