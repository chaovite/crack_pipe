function y = imex_ark4(y0,t,h,f,A)
% [y t] = imex_ark4(y0,t0,n,f,g)
% 4th-order additive Runge-kutta time integrator.
% This function computes the solution of the ordinary differential equation:
% y' = f(y(t),t) + A*y
% where f(y(t),t) is non-stiff and A*y is stiff.
% y0 : initial condition y0=y(t0)
% t0 : initial time
% f  : f(t,y) explicit part
% A  : A (matrix) implicit part

% IMEX 4
% Implementation of
% ARK4(3)6L[2]SA–ERK
% ARK4(3)6L[2]SA–ESDIRK
% Christopher A. Kennedy and Mark H. Carpenter (2003)
% Additive Runge–Kutta schemes for convection–diffusion–reaction equations  

c1 = 0;
c2 = 1/2;
c3 = 83/250;
c4 = 31/50;
c5 = 17/20;
c6 = 1;

b1 = 82889/524892;
b2 = 0;
b3 = 15625/83664;
b4 = 69875/102672;
b5 = -2260/8211;
b6 = 1/4;

a21 = 1/2;

a31 = 13861/62500;
a32 = 6889/62500;

a41 = -116923316275/ 2393684061468;
a42 = -2731218467317/15368042101831;
a43 = 9408046702089/11113171139209;

a51 = -451086348788/ 2902428689909;
a52 = -2682348792572/7519795681897;
a53 = 12662868775082/11960479115383;
a54 = 3355817975965/11060851509271;

a61 = 647845179188/3216320057751;
a62 = 73281519250/8382639484533;
a63 = 552539513391/3454668386233;
a64 = 3354512671639/8306763924573;
a65 = 4040/17871;

ia21=1/4;

ia31=8611/62500;
ia32=-1743/31250;

ia41=5012029/ 34652500;
ia42=-654441/2922500;
ia43=174375/388108;

ia51=15267082809/155376265600;
ia52=-71443401/120774400;
ia53=730878875/902184768;
ia54=2285395/8070912;

ia61=82889/524892;
ia62=0;
ia63=15625/83664;
ia64=69875/102672;
ia65=-2260/8211;     



I = speye(size(A));
a = 1/4;
B = I - a*h*A;

  k1 = y0;
  k2 = B\(y0 + a21*h*f(k1,t+c1*h) + ia21*h*A*k1);
  k3 = B\(y0 + a31*h*f(k1,t+c1*h) + a32*h*f(k2,t+c2*h) + ...
       ...
               ia31*h*A*k1 + ia32*h*A*k2);
  k4 = B\(y0 + a41*h*f(k1,t+c1*h) +  a42*h*f(k2,t+c2*h) + ...
               a43*h*f(k3,t+c3*h) + ...
               ...
               ia41*h*A*k1 + ia42*h*A*k2 + ia43*h*A*k3);
  k5 = B\(y0 + a51*h*f(k1,t+c1*h)   + a52*h*f(k2,t+c2*h) + ...
               a53*h*f(k3,t+c3*h)   + a54*h*f(k4,t+c4*h) + ...
                ... 
               ia51*h*A*k1 + ia52*h*A*k2 + ia53*h*A*k3 + ia54*h*A*k4);
  k6 = B\(y0 + a61*h*f(k1,t+c1*h)   + a62*h*f(k2,t+c2*h) + ...
               a63*h*f(k3,t+c3*h)   + a64*h*f(k4,t+c4*h) + ...
               a65*h*f(k5,t+c5*h)   + ...
                ...
               ia61*h*A*k1 + ia63*h*A*k3 + ia64*h*A*k4 + ...
               ia65*h*A*k5);

  y =  y0 + b1*h*f(k1,t+c1*h) + b3*h*f(k3,t+c3*h) + ...
            b4*h*f(k4,t+c4*h) + b5*h*f(k5,t+c5*h) + ...
            b6*h*f(k6,t+c6*h) + ...
            ...
            b1*h*A*k1 + b3*h*A*k3 + ...
            b4*h*A*k4 + b5*h*A*k5 + b6*h*A*k6;




end
