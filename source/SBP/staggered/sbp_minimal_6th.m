function [xp,xm,Pp,Pm,Qp,Qm] = ...
         sbp_staggered_6th(n,h,r,x)

test = false;

% Free parameters 
if nargin < 4
 % Optimized for spectral radius
 if r == 0 
   q55 = -1.080520271717843;
   q56 = 1.133105138515982; 
 elseif r == 1
   q55 = -1.067407245522082;
   q56 = 1.134379933657703; 
  else
    q55 = -1.008115141072967;
    q56 = 1.120566024128776; 
  end
else
  q55 = x(1);
  q56 = x(2);
end

% Coefficients determined such that the SBP property is satisfied
% free parameters: q55, q56
q11 = 16*q55/5 + 64*q56/5 - 12176791/1020600; 
q12 = -7*q55 - 27*q56 + 15921763/663552; 
q13 = 7*q55 + 24*q56 - 48888389/2488320; 
q14 = -21*q55/5 - 54*q56/5 + 21040763/2764800; 
q15 = q55 + 6396463/5806080; 
q16 = q56 - 33696889/29859840; 

q21 = -64*q55/5 - 256*q56/5 + 4996231/113400; 
q22 = 28*q55 + 108*q56 - 102978463/1105920; 
q23 = -28*q55 - 96*q56 + 7319477/92160; 
q24 = 84*q55/5 + 216*q56/5 - 28102619/921600; 
q25 = -4*q55 - 8513639/1935360; 
q26 = -4*q56 + 15023011/3317760; 

q31 = 96*q55/5 + 384*q56/5 - 15041203/226800; 
q32 = -42*q55 - 162*q56 + 152928983/1105920; 
q33 = 42*q55 + 144*q56 - 2186833/18432; 
q34 = -126*q55/5 - 324*q56/5 + 43285411/921600; 
q35 = 6*q55 + 12649207/1935360; 
q36 = 6*q56 - 22637531/3317760; 

q41 = -64*q55/5 - 256*q56/5 + 45176119/1020600; 
q42 = 28*q55 + 108*q56 - 305991767/3317760; 
q43 = -28*q55 - 96*q56 + 195283253/2488320; 
q44 = 84*q55/5 + 216*q56/5 - 87800627/2764800; 
q45 = -4*q55 - 19153951/5806080;
q46 = -4*q56 + 135415537/29859840; 

q51 = 16*q55/5 + 64*q56/5 - 4429/400; 
q52 = -7*q55 - 27*q56 + 14763/640; 
q53 = 7*q55 + 24*q56 - 12551/640; 
q54 = -21*q55/5 - 54*q56/5 + 2303/300; 
 
% P+
pp1 =   95/288; 
pp2 =  317/240; 
pp3 =   23/30; 
pp4 =  793/720;
pp5 =  157/160; 
 
% P-
pm1 =    5251/68040; 
pm2 =  133801/138240; 
pm3 =   90101/103680; 
pm4 =   26369/23040;
pm5 =  224201/241920; 
pm6 = 1262401/1244160; 

% Interior stencil
s1 = -  3/640;
s2 =   25/384;
s3 = - 75/64;
s4 =   75/64;
s5 = - 25/384;
s6 =    3/640;  

% Q+ top-left corner
QpL = [
      q11 q12 q13 q14 q15 q16 0  0 ;
      q21 q22 q23 q24 q25 q26 0  0 ;
      q31 q32 q33 q34 q35 q36 0  0 ;
      q41 q42 q43 q44 q45 q46 s6 0 ;
      q51 q52 q53 q54 q55 q56 s5 s6;
      ];

% Q- top-left corner
QmL = -[
      q11 q21 q31 q41 q51 0  0  0 ;
      q12 q22 q32 q42 q52 0  0  0 ;
      q13 q23 q33 q43 q53 0  0  0 ;
      q14 q24 q34 q44 q54 s1 0  0 ;
      q15 q25 q35 q45 q55 s2 s1 0 ;
      q16 q26 q36 q46 q56 s3 s2 s1;
      ];
QmL(1,1) = QmL(1,1) -1;


% Construct (n+2 x n+1) matrix Qp with interior stencil
w = 6; 
s = rot90(vander(1:w))\((0:(w-1)).*(w/2-1/2+1).^([0 0:w-2]))';  
Qp = spdiags(repmat(-s(end:-1:1)',[n+2 1]), -(w/2-1):w/2, n+2, n+2); 
Qp(end,:) = [];
% Add SBP boundary closures
Qp(1:5,1:8) = QpL;
Qp(end-4:end,end-7:end) = -fliplr(flipud(QpL));
% Repeat for Q-
Qm = spdiags(repmat(s(:)',[n+2 1]), -(w/2-1)-1:w/2-1, n+2, n+2);
Qm(:,end) = [];
Qm(1:6,1:8) = QmL;
Qm(end-5:end,end-7:end) = -fliplr(flipud(QmL));

% P+ and P+^-1
Pp  =   h*spdiags(   [pp1 pp2 pp3 pp4 pp5 ones(1,n-9) pp5 pp4 pp3 pp2 pp1]',0,n+1,n+1);
Ppi = 1/h*spdiags(1./[pp1 pp2 pp3 pp4 pp5 ones(1,n-9) pp5 pp4 pp3 pp2 pp1]',0,n+1,n+1); 

% P- and P-^-1
Pm  =   h*spdiags(   [pm1 pm2 pm3 pm4 pm5 pm6 ones(1,n-10) pm6 pm5 pm4 pm3 pm2 pm1]',0,n+2,n+2);
Pmi = 1/h*spdiags(1./[pm1 pm2 pm3 pm4 pm5 pm6 ones(1,n-10) pm6 pm5 pm4 pm3 pm2 pm1]',0,n+2,n+2);

% D+ and D-
Dp = Ppi*Qp;
Dm = Pmi*Qm;

% nodal and cell-centered grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  

% Test operators
test = false;
if test
  for j=0:4
    disp([ 'Dp, j = ' num2str(j) ' Error max = ' ...
    num2str(max(abs(Qp*xm.^j-j*Pp*xp.^max([j-1,0]))))]);
    disp([ 'Dm, j = ' num2str(j) ' Error max = '...
    num2str(max(abs(Qm*xp.^j-j*Pm*xm.^max([j-1,0]))))]);
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
end


