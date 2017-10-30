function [xp,xm,Pp,Pm,Qp,Qm] = sbp_half_4th(n,h,test)
% [xp,xm,Pp,Pm,Qp,Qm] = sbp_half_4th(n,h,x)
% n : number of grid points  (n+1) nodal grid (n+2) cell-centered grid
% h : grid spacing

if nargin < 3
  test = false;
end


% Unknown coefficients
x(1) = -0.742829358845154;
x(2) = 1.2121;
pm4 = x(2);
qm43 = x(1);


% Coefficients determined such that the SBP property is satisfied
qp22 = -317*pm4/21 + 9*qm43 + 37447/1512;
qm31 = 110*pm4/7 - 9*qm43 - 26333/1008;
pp3 = -16*pm4/21 + 2783/1512;
qp02 = -26*pm4/7 + 3*qm43 + 1853/252;
qm30 = -110*pm4/21 + 3*qm43 + 26375/3024;
qp03 = 110*pm4/21 - 3*qm43 - 26375/3024;
qp31 = -61*pm4/21 + qm43 + 12739/3024;
qm21 = -317*pm4/21 + 9*qm43 + 37447/1512;
qp14 = 5*pm4 - 3*qm43 - 67/8;
pm2 = 83*pm4/21 - 4093/1512;
qm10 = -10*pm4/21 + qm43 + 3799/3024;
qp30 = 8*pm4/7 - 151/126;
pm3 = -23*pm4/7 + 2125/504;
qm20 = 26*pm4/7 - 3*qm43 - 1853/252;
qm11 = 44*pm4/7 - 3*qm43 - 3253/336;
pp0 = 16*pm4/21 - 88/189;
qm42 = 3*pm4 - 3*qm43 - 19/3;
qm03 = -8*pm4/7 + 151/126;
qp01 = 10*pm4/21 - qm43 - 3799/3024;
pp1 = -16*pm4/7 + 1859/504;
qp12 = 317*pm4/21 - 9*qm43 - 37447/1512;
qm02 = 64*pm4/21 - 1271/378;
qm22 = 317*pm4/21 - 9*qm43 - 37447/1512;
qp11 = -44*pm4/7 + 3*qm43 + 3253/336;
qm23 = -26*pm4/7 + 3*qm43 + 1853/252;
qm13 = 61*pm4/21 - qm43 - 12739/3024;
qp20 = -64*pm4/21 + 1271/378;
qm12 = -61*pm4/7 + 3*qm43 + 12739/1008;
qp10 = 40*pm4/21 - 1007/378;
qm40 = 2*pm4 - qm43 - 25/8;
qp04 = -2*pm4 + qm43 + 25/8;
qm32 = -87*pm4/7 + 9*qm43 + 7333/336;
qm41 = -5*pm4 + 3*qm43 + 67/8;
qp21 = 61*pm4/7 - 3*qm43 - 12739/1008;
qm33 = 41*pm4/21 - 3*qm43 - 13247/3024;
qp33 = -41*pm4/21 + 3*qm43 + 13247/3024;
qm01 = -40*pm4/21 + 1007/378;
qp24 = -3*pm4 + 3*qm43 + 19/3;
qp32 = 26*pm4/7 - 3*qm43 - 1853/252;
pm1 = -17*pm4/7 + 745/252;
qm00 = -1/2;
qp34 = -qm43;
qp00 = -1/2;
qp13 = -110*pm4/7 + 9*qm43 + 26333/1008;
pm0 = 16*pm4/21 - 88/189;
pp2 = 16*pm4/7 - 197/126;
qp23 = 87*pm4/7 - 9*qm43 - 7333/336;



% Number of coefficients
b = 4;

% Q+ and Q-, top-left corner
QpL = [...
qp00, qp01, qp02, qp03, qp04;
 qp10, qp11, qp12, qp13, qp14;
 qp20, qp21, qp22, qp23, qp24;
 qp30, qp31, qp32, qp33, qp34
];
QmL = [...
qm00, qm01, qm02, qm03;
 qm10, qm11, qm12, qm13;
 qm20, qm21, qm22, qm23;
 qm30, qm31, qm32, qm33;
 qm40, qm41, qm42, qm43
];

% Q+ and Q-
w = b; 
s = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Qp = spdiags(repmat(-s(end:-1:1)',[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qp(end,:) = [];
Qm(:,end) = [];

% Add SBP boundary closures
Qp(1:b,1:b+1) = QpL;
Qp(end-b+1:end,end-b:end) = -fliplr(flipud(QpL));
Qm(1:b+1,1:b) = QmL;
Qm(end-b:end,end-b+1:end) = -fliplr(flipud(QmL));

% P+ and P-
Pp = ones(n+1,1);
Pm = ones(n+2,1);

Pp(1:b) = [pp0,  pp1,  pp2,  pp3]; 
Pp(end-b+1:end) = Pp(b:-1:1);
Pm(1:b+1) = [pm0,  pm1,  pm2,  pm3,  pm4];
Pm(end-b:end) = Pm(b+1:-1:1);
Pp = spdiags(Pp,0,n+1,n+1);
Pm = spdiags(Pm,0,n+2,n+2);

Pp = h*Pp;
Pm = h*Pm;

% nodal and cell-centered grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  


% Test operators
if test
for j=0:b/2
  disp([ 'Dp, j = ' num2str(j) ' Error max = ' ...
  num2str(max(abs(Qp*xm.^j-j*Pp*xp.^max([j-1,0]))))]);
  disp([ 'Dm, j = ' num2str(j) ' Error max = '...
  num2str(max(abs(Qm*xp.^j-j*Pm*xm.^max([j-1,0]))))]);
end  
end

