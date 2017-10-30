function [T_out] = ...
  compute_transfer_function(data_in,T,dt,wlen,run_tests)
% [T_out status] = compute_transfer_function(data_in,T,dt,wlen)
% 
% Computes the dimensionless transfer function 'F' relating pressure 'p' (input) to velocity 'u'
% (output). Thus,  'u = F*p', where * is the convolution operation. 
%
% Input:
%   data_in: A struct that contains time series data t (time), p (pressure), and u (velocity). 
%            This struct also contains 'Z', which is the impedance (see
%            load_transfer_function_data).           
%         T: maximum signal duration (seconds). 
%            Signals longer than T are truncated and signals shoer than T are zero-padded.
%        dt: Size of constant sample spacing in time, which is used to control
%            the frequency content in the transfer function. Provide 0 or leave
%            parameter values as empty to disable.
%      wlen: window length to use when tapering the signal (unit in time).
% run_tests: Run verification tests to ensure that the transfer function satisfies certain 
%            mathematical properties.
%
% Output: 
%  T_out: A struct that contains the frequency vector 'f' and transfer
%         function 'F' as well as the settings used to construct the transfer
%         function.
% status: A struct that contains the status of the verifications tests if they
% ran.
% 
% See below for more details.
%
% Duration and resampling
% The parameter 'T' (seconds) determines the length of the time series. If 'T' is
% shorter than the total duration of the time series (computed from data.t(end))
% then the time series is truncated. If 'T' is longer than the total duration of
% the time series, then the time series is augmented by zero padding at the end.
% The longer the duration is, the more resolved the transfer function becomes.
%
% The time series can be resampled by changing the sampling size 'dt' in time
% (seconds). A smaller size yields a transfer function with more frequency
% content.
% Windowing
% The parameter 'wlen' controls the window length which will apply a cosine window
% tapering of length 'wlen' (units in time) to the end of the time series.

% Verification
% The use of padding and windowing can destroy important properties of the
% transfer function that are needed to expect reasonable physical behavior.
% Therefore, the following tests are performed:
% 1. Stability
% 2. Consistency  
% The status of each test is stored in the output argument 'status'.

u = data_in.u;
p = data_in.p;
t = data_in.t;
Z = data_in.Z;
 
% Window signal
u=coswindow(u,t,[0 t(2) t(end)-wlen t(end)]);
p=coswindow(p,t,[0 t(2) t(end)-wlen t(end)]);

% Pad signal
if ~isempty(T) & T>t(end)
    u = zpad(u,round(T/t(2)));
    p = zpad(p,round(T/t(2)));
end

% resample the signal in time domain
if (~isempty(dt) && dt ~= 0)
  tol = 0.1*min(t(2),dt);
  [a,b] = rat(t(2)/dt,tol);
  u = resample(u,a,b);
  p = resample(p,a,b);
  t = dt*[0:length(u)-1];
else
  dt = t(2);
end

% Make the signal have an even number of samples
% Remove the last sample if the number of samples of the signal is odd
if mod(numel(u),2)==1 
  u(end) = [];
  p(end) = [];
  t(end) = [];
end

% Construct transfer function
[uhat f] = fft_dim(u,dt); 
[phat f] = fft_dim(p,dt); 
F = Z*uhat./phat; 

% Verification
if run_tests
  status.is_stable = transfer_function_stability_test(f,F,min(f),max(f));
  status.is_consistent = transfer_function_dc_test(f,F);
  status.is_energy_stable = transfer_function_energy_test(F);
  check('Stability',status.is_stable);
  check('Consistency',status.is_consistent);
  check('Energy',status.is_energy_stable);     
end

T_out.f = f;
T_out.F = F;
T_out.t = t;
T_out.u = u;
T_out.p = p;
T_out.uhat = uhat;
T_out.phat = phat;
T_out.data = data_in;
T_out.status = status;
end
