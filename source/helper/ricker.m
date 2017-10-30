function [rw] = ricker(t, f, t0)
%RICKER creates an causal ricker wavelet signal
%
%   RICKER creates and plots a default causal ricker wavelet with:
%
%       peak frequency   = 20 Hz
%       sampling time    = 0.001 seconds
%       number of points = 100;
%       peak location    = 1/F = 1/20Hz
%
tau = t-t0;
rw = (1-2*tau.*tau*f^2*pi^2).*exp(-tau.^2*pi^2*f^2);
end
