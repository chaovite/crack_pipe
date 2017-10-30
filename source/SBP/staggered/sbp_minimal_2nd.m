function [spectral_radius,Pp,Pm,Ppi,Pmi,Qp,Qm,Dp,Dm,x,y] = sbp_staggered_2nd(n,h,x)

% Coefficients determined such that the SBP property is satisfied

q12 = x(1);

q11 = -q12;
q22 = -q12;
q21 =  q12 - 1;
 
% P+
pp1 =  q12/2;
pp2 = -q12/2 + 3/2;
 
% P-
pm1 = -q12 + 1;
pm2 =  q12;

% Interior stencil
s1 = -1;
s2 =  1;

% Q+ top-left corner
QpL = [
      q11 q12 0   0;
      q21 q22 s2  0;
      0   0   s1 s2;
      ];

% Q- top-left corner
QmL = -[
      q11 q21 0 ;
      q12 q22 0 ;
      0   s2  s1;
      0   0   s2;
      ];
QmL(1,1) = QmL(1,1) -1;

% Construct (n+2 x n+1) matrix Qp with interior stencil
w = 2; 
s = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Qp = spdiags(repmat(-s(end:-1:1)',[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qp(end,:) = [];
% Add SBP boundary closures
Qp(1:3,1:4) = QpL;
Qp(end-2:end,end-3:end) = -fliplr(flipud(QpL));
% Repeat for Q-
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qm(:,end) = [];
Qm(1:4,1:3) = QmL;
Qm(end-3:end,end-2:end) = -fliplr(flipud(QmL));

% P+ and P+^-1
Pp  =   h*spdiags(   [pp1 pp2 ones(1,n-3) pp2 pp1]',0,n+1,n+1);
Ppi = 1/h*spdiags(1./[pp1 pp2 ones(1,n-3) pp2 pp1]',0,n+1,n+1); 

% P- and P-^-1
Pm  =   h*spdiags(   [pm1 pm2 ones(1,n-2) pm2 pm1]',0,n+2,n+2);
Pmi = 1/h*spdiags(1./[pm1 pm2 ones(1,n-2) pm2 pm1]',0,n+2,n+2);

% D+ and D-
Dp = Ppi*Qp;
Dm = Pmi*Qm;

% nodal and cell-centered grids
x = h*[0:n]';
y = h*[0 1/2+0:n n]';  

% Test operators
test = true;
if test
for j=0:2
  disp([ 'Dp, j = ' num2str(j) ' Error max = ' ...
  num2str(max(abs(Qp*y.^j-j*Pp*x.^max([j-1,0]))))]);
  disp([ 'Dm, j = ' num2str(j) ' Error max = '...
  num2str(max(abs(Qm*x.^j-j*Pm*y.^max([j-1,0]))))]);
end  
end
  
disp('Leading Error Term');
disp(max([abs(Pp\Qp*y.^j-j*x.^max([j-1,0]));...
          abs(Pm\Qm*x.^j-j*y.^max([j-1,0]))])) 

% Minimize spectral radius
zp = spalloc(n+1,1,0);
zm = spalloc(n+2,1,0);
A = [zp*zp', Dp;
         Dm, zm*zm'];

lambda = eig(full(h*A));
disp('Spectral radius');
spectral_radius = max(abs(lambda));
disp([spectral_radius]);
