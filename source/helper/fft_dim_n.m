function [G,f] = fft_dim_n(g,dt)
% Notation of Fourier transform: 
%   ghat = int_{-inf}^{+inf} g(t)*exp(-i*omega*t) dt  

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
