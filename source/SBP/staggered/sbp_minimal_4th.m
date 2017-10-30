function [xp,xm,Pp,Pm,Qp,Qm,alpha,bet] = ...
         sbp_staggered_4th(n,h)




% Coefficients determined such that the SBP property is satisfied

q11 = - 14/15;
q12 =   1;
q13 = - 1/24;
q14 = - 1/40;

q21 = - 4/45;
q22 = - 1;
q23 =   19/18;
q24 =   1/30;

q31 =   1/45;
q32 =   0;
q33 = - 19/18;
q34 =   43/40;
 
% P+
pp1 = 3/8;
pp2 = 7/6;
pp3 = 23/24;
 
% P-
pm1 = 2/45;
pm2 = 1;
pm3 = 67/72;
pm4 = 41/40;

% Interior stencil
s1 = -  1/24;
s2 =    9/8;
s3 = -  9/8;
s4 =    1/24;

% Q+ top-left corner
QpL = [
      q11 q12 q13 q14 0  0 ;
      q21 q22 q23 q24 0  0 ;
      q31 q32 q33 q34 s1 0 ;
      ];

% Q- top-left corner
QmL = -[
      q11 q21 q31 0  0 ;
      q12 q22 q32 0  0 ;
      q13 q23 q33 s4 0 ;
      q14 q24 q34 s3 s4;
      ];
QmL(1,1) = QmL(1,1) -1;

% Construct (n+2 x n+1) matrix Qp with interior stencil
w = 4; 
s = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Qp = spdiags(repmat(-s(end:-1:1)',[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qp(end,:) = [];
% Add SBP boundary closures
Qp(1:3,1:6) = QpL;
Qp(end-2:end,end-5:end) = -fliplr(flipud(QpL));
% Repeat for Q-
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qm(:,end) = [];
Qm(1:4,1:5) = QmL;
Qm(end-3:end,end-4:end) = -fliplr(flipud(QmL));

% P+ and P+^-1
Pp  =   h*spdiags(   [pp1 pp2 pp3 ones(1,n-5) pp3 pp2 pp1]',0,n+1,n+1);
Ppi = 1/h*spdiags(1./[pp1 pp2 pp3 ones(1,n-5) pp3 pp2 pp1]',0,n+1,n+1); 

% P- and P-^-1
Pm  =   h*spdiags(   [pm1 pm2 pm3 pm4 ones(1,n-6) pm4 pm3 pm2 pm1]',0,n+2,n+2);
Pmi = 1/h*spdiags(1./[pm1 pm2 pm3 pm4 ones(1,n-6) pm4 pm3 pm2 pm1]',0,n+2,n+2);

% D+ and D-
Dp = Ppi*Qp;
Dm = Pmi*Qm;

% nodal and cell-centered grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  

% Test operators
test = false;
if test
  for j=0:2
    disp([ 'Dp, j = ' num2str(j) ' Error max = ' ...
    num2str(max(abs(Qp*xm.^j-j*Pp*xp.^max([j-1,0]))))]);
    disp([ 'Dm, j = ' num2str(j) ' Error max = '...
    num2str(max(abs(Qm*xp.^j-j*Pm*xm.^max([j-1,0]))))]);
  end  
    
  disp('Leading Error Term');
  disp(max([abs(Pp\Qp*xm.^j-j*xp.^max([j-1,0]));...
            abs(Pm\Qm*xp.^j-j*xm.^max([j-1,0]))])) 
end
