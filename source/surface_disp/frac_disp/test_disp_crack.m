% Compare different model for horizontal openning crack

x=zeros(1,1000);% East
y=linspace(0,1e3,1000); %North
station_loc = [x;y]; % station location
H=1000; % depth of crack
R=100; % crack radius
W=1; % width of the crack
strike=0; % strike
dip=0;  % dip
p=1e4; % crack hydrostatic pressure
mu=1e10; % shear modulus
nu=0.25; % Poisson's ratio. 

[U1_okada, U2_okada, U3_okada] = disp_crack(H, R, W, strike, dip, station_loc, p, mu, nu, 'Okada');
[U1_fialko, U2_fialko, U3_fialko] = disp_crack(H, R, W, strike, dip, station_loc, p, mu, nu, 'Fialko2001');
[U1_sun, U2_sun, U3_sun] = disp_crack(H, R, W, strike, dip, station_loc, p, mu, nu, 'Fialko2001');
%%
close all;
plot(y'/R, U2_okada,'k',y'/R, U2_fialko,'r*',y'/R, U2_sun,'b');
h=legend({'okada','fialko2001','sun69'});
set(h,'fontsize',16);
xlabel('distance/R','fontsize',18);
ylabel('Ur (m)','fontsize',18);
set(gca,'fontsize',18);
title(['Compare Ur H/R=',num2str(H/R)]);

figure;
plot(y'/R, U3_okada,'k',y'/R, U3_fialko,'r*',y'/R, U3_sun,'b');
h=legend({'okada','fialko2001','sun69'});
set(h,'fontsize',16);
xlabel('distance/R','fontsize',18);
ylabel('Ur (m)','fontsize',18);
set(gca,'fontsize',18);
title(['Compare Ur H/R=',num2str(H/R)]);



