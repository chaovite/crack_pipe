f0                  = 0;           % frequency for input source
fn                  = 0.1*500^2;
tmax                = 0.5*2;
A                   = 1;
t0                  = 0.4;
tc                  = 0.2;
phi                 = 0;
k                   = fn/tmax;   
shutoff_time        = -1;          % time when input source is turned off
bcLs                = 300/tmax;
dt=1e-4;

% Plot chirp/delta

    % Time specifications:
   Fs = 1e3/dt;                  % samples per second
   deltat = dt;                     % seconds per sample
   Fs = 1/dt;
   StopTime = 10;
   t = (0:deltat:StopTime-deltat)';
   N = size(t,1);
   % Sine wave:
   Fc = 12;                       % hertz
   phi = 0;
   p = A*sin(phi + 2*pi*(f0*t + k*t.^2/2)).*exp(-bcLs*t);
   % Fourier Transform:
   X = fftshift(fft(p));
   % Frequency specifications:
   dF = Fs/N;                      % hertz
   f = -Fs/2:dF:Fs/2-dF;           % hertz
   % Plot the spectrum:
   figure(1);
   plot(f,abs(X)/N);
   %xlim([0 0.1])
   xlabel('Frequency (in hertz)');
   title('Magnitude Response');
   figure(2)
   plot(t,p)
   %xlim([0 1]);
   xlabel('s');  