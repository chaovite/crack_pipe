function [g,t] = ifft_dim_n(G,df)
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