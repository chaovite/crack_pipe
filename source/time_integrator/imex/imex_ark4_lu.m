function y = imex_ark4_lu(y0,t,dt,f,A,L,U,p,q)
% [y t] = imex_ark4(y0,t,dt,f,A,L,U)
% 4th-order additive Runge-kutta time integrator.
% This function computes the solution of the ordinary differential equation:
% y' = f(y(t),t) + A*y
% where f(y(t),t) is non-stiff and A*y is stiff.
% y0 : initial condition y0=y(t0)
% dt : time step
% f  : f(t,y) explicit part
% A  : A (matrix) implicit part
% L,U : LU-decompostion of A-a*h*I (see imex_ark4_get_lu)

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

h = dt;
k1 = y0;


f1 = f(k1,t+c1*h);

hAk1 = A*(h*k1);

b = (y0 + a21*h*f1 + ia21*hAk1);
k2(q) = U\(L\b(p));
k2 = k2';

f2 = f(k2,t+c2*h);
hAk2=A*(h*k2);
b  =(y0 + a31*h*f1 + a32*h*f2 + ...
     ...
             ia31*hAk1 + ia32*hAk2);
k3(q) = U\(L\b(p));
k3 = k3';

f3 = f(k3,t+c3*h);
hAk3=A*(h*k3);
b = (y0 + a41*h*f1 +  a42*h*f2 + ...
             a43*h*f3 + ...
             ...
             ia41*hAk1 + ia42*hAk2 + ia43*hAk3);
k4(q) = U\(L\b(p));
k4 = k4';

f4 = f(k4,t+c4*h);
hAk4=A*(h*k4);
b = (y0 + a51*h*f1   + a52*h*f2 + ...
             a53*h*f3   + a54*h*f4 + ...
              ... 
             ia51*hAk1 + ia52*hAk2 + ia53*hAk3 + ia54*hAk4);
k5(q) = U\(L\b(p));
k5 = k5';
f5 = f(k5,t+c5*h);

hAk5=A*(h*k5); 
b =(y0 + a61*h*f1   + a62*h*f2 + ...
             a63*h*f3   + a64*h*f4 + ...
             a65*h*f5   + ...
              ...
             ia61*hAk1 + ia63*hAk3 + ia64*hAk4 + ...
             ia65*hAk5);
k6(q) = U\(L\b(p));
k6 = k6';
f6 = f(k6,t+c6*h);
hAk6=A*(h*k6);
y =  y0 + b1*h*f1 + b3*h*f3 + ...
          b4*h*f4 + b5*h*f5 + ...
          b6*h*f6 + ...
          ...
          b1*hAk1 + b3*hAk3 + ...
          b4*hAk4 + b5*hAk5 + b6*hAk6;

end
