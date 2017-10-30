function [U D S flag] = disloc3dpm(mdl,stat,mu,nu)
%[U D S flag] = disloc3dpm(m,x,mu,nu) [pure Matlab version of disloc3d]
%
% Returns the deformation at point 'x', given dislocation model 'm'.  'mu'
% specifies the shear modulus and 'nu' specifies Poisson's ratio.
%
% Both 'm' and 'x' can be matrices, holding different models and observation
% coordinates in the columns.  In this case, the function returns the
% deformation at the points specified in the columns of 'x' from the sum of
% all the models in the columns of 'm'.  'x' must be 3xi (i = number of
% observation coordinates) and 'm' must be 10xj (j = number of models).
%
% The coordinate system is as follows:
%   east = positive X, north = positive Y, up = positive Z.
% Observation coordinates with positive Z values return a warning.
%
% The outputs are 'U', the three displacement components: east, north, and up
% (on the cols); 'D', the nine spatial derivatives of the displacement: Uxx,
% Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, and Uzz (on the cols); and 'S', the 6
% independent stress components: Sxx, Sxy, Sxz, Syy, Syz, and Szz (on the
% cols). All these outputs have the same number of columns as 'x'.
%
% Output 'flag' is set for a singularity.
%
% The dislocation model is specified as:
%   length width depth dip strike east north strike-slip dip-slip opening.
% The coordinates (depth, east, north) specify a point at the middle of the
% bottom edge of the fault for positive dips and the middle of the top edge
% for negative dips. The coordinates (east, north) are specified just as for
% obs; however, note carefully that depth is a positive number in 'm', while
% in 'x' the equivalent coordinate is z and so is negative.
%
% Mind your units: for example, if lengths are given in km, then so should be
% slips.
%
% ---------------
% This is a pure-Matlab version of disloc3d based on the Matlab translation
% of Y. Okada's dc3d.f found in Coulomb 3.1. The citations are as follows:
%   Toda, S., R. S. Stein, K. Richards-Dinger and S. Bozkurt (2005),
%     Forecasting the evolution of seismicity in southern California:
%     Animations built on earthquake stress transfer, J. Geophys. Res.,
%     B05S16, doi:10.1029/2004JB003415.
%   Lin, J., and R.S. Stein, Stress triggering in thrust and subduction
%     earthquakes, and stress interaction between the southern San Andreas
%     and nearby thrust and strike-slip faults, J. Geophys. Res., 109,
%     B02303, doi:10.1029/2003JB002607, 2004.
%   Okada, Y., Internal deformation due to shear and tensile faults in a
%     half-space, Bull. Seismol. Soc. Amer., 82 (2), 1018-1040, 1992.

% 20 Oct 2010. AMB. Initial version.
  
  global N_CELL
  N_CELL = 1;
  
  lambda = 2*mu*nu/(1-2*nu);
  
  n = size(stat,2);
  nm = size(mdl,2);
  
  U = zeros(3,n);
  D = zeros(9,n);
  S = zeros(6,n);
  flag = zeros(1,n);

  for(j = 1:nm) % dislocation elements
    m = mdl(:,j);

    alpha = (lambda + mu)/(lambda + 2*mu);
    
    strike = m(5) - 90;
    cs = cosd(strike);
    ss = sind(strike);
    R = [ cs ss 0
	 -ss cs 0
	  0  0  1];

    dip = m(4);
    cd = cosd(dip);
    sd = sind(dip);
    
    disl1 = m(8);
    disl2 = m(9);
    disl3 = m(10);
    al1 = m(1)/2;
    al2 = al1;
    aw1 = m(2)/2;
    aw2 = aw1;

    for(i = 1:n) % observation points
      s = stat(:,i);
      
      % disloc3d coords -> dc3d coords
      x = cs*(-m(6) + s(1)) - ss*(-m(7) + s(2));
      y = -cd*m(2)/2 + ss*(-m(6) + s(1)) + cs*(-m(7) + s(2));
      z = s(3);
      depth = m(3) - 0.5*m(2)*sd;
    
      [Ux Uy Uz Uxx Uyx Uzx Uxy Uyy Uzy Uxz Uyz Uzz iret] = Okada_DC3D(...
	  alpha,x,y,z,depth,dip,...
	  al1,al2,aw1,aw2,disl1,disl2,disl3);

      % dc3d coords -> disloc3d coords and accumulate
      flag(i) = flag(i) | iret;
      U(:,i) = U(:,i) + R*[Ux Uy Uz]';
      ud = R*[Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz]'*R';
      D(:,i) = D(:,i) + ud(:);
    end
  end
        
  % stress
  for(i = 1:n) % observation points
    Ud = reshape(D(:,i),3,3);
    Ud = (Ud + Ud')/2; % symmetrize
    s = 2*mu*Ud + lambda*sum(diag(Ud))*eye(3); % isotropic form of Hooke's law
    S(:,i) = s([1:3 5:6 9]); % six unique elements due to symmetry
  end

  clear global N_CELL ...
      ALE    ALP4   CD     ET2    FZ     HZ     R3     SDCD   X32    Y32 ...
      ALP1   ALP5   CDCD   EY     GY     Q2     R5     SDSD   XI2    ...
      ALP2   ALX    D      EZ     GZ     R      S2D    TT     Y      ...
      ALP3   C2D    DUMMY  FY     HY     R2     SD     X11    Y11;

% ------------------------------------------------------------------------------
% The rest of the code is from Coulomb 3.1.
% ------------------------------------------------------------------------------
function[UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET] = Okada_DC3D(ALPHA,...
                X,Y,Z,DEPTH,DIP,...
                AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3)
            
%       IMPLICIT REAL*8 (A-H,O-Z)                                         04640000
%       REAL*4   ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3, 04650000
%      *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             04660000
% C                                                                       04670000
% C********************************************************************   04680000
% C*****                                                          *****   04690000
% C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   04700000
% C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   04710000
% C*****                         CODED BY  Y.OKADA ... SEP 1991   *****   04720002
% C*****                         REVISED   Y.OKADA ... NOV 1991   *****   04730002
% C*****                                                          *****   04740000
% C********************************************************************   04750000
% C                                                                       04760000
% C***** INPUT                                                            04770000
% C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           04780000
% C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          04790000
% C*****   DEPTH : SOURCE DEPTH                                           04800000
% C*****   DIP   : DIP-ANGLE (DEGREE)                                     04810000
% C*****   AL1,AL2   : FAULT LENGTH (-STRIKE,+STRIKE)                     04820000
% C*****   AW1,AW2   : FAULT WIDTH  ( DOWNDIP, UPDIP)                     04830000
% C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              04840000
% C                                                                       04850000
% C***** OUTPUT                                                           04860000
% C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               04870000
% C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             04880000
% C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )04890000
% C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     04900000
% C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  04910002
% C                                                                       04920000
    global DUMMY SD CD
    global XI2 ET2 Q2 R
    global N_CELL
    
%       COMMON /C0/DUMMY(5),SD,CD                                         04930000
%       COMMON /C2/XI2,ET2,Q2,R                                           04940000
%       DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12)                   04950000
%       DATA  F0/0.D0/                                                    04960000

%    F0 = double(0.0);
F0 = zeros(N_CELL,1,'double');
%C----- 
if Z>0.0
    disp(' ** POSITIVE Z WAS GIVEN IN SUB-DC3D');
end
%for I=1:1:12
        U    = zeros(N_CELL,12,'double');
        DUA  = zeros(N_CELL,12,'double');
        DUB  = zeros(N_CELL,12,'double');
        DUC  = zeros(N_CELL,12,'double');
        IRET = zeros(N_CELL,1,'double');
%         U  (1:N_CELL,1:12)=0.0;
%         DUA(1:N_CELL,1:12)=0.0;
%         DUB(1:N_CELL,1:12)=0.0;
%         DUC(1:N_CELL,1:12)=0.0;
%         IRET(1:N_CELL,1:12)=0.0;
%end
      AALPHA=ALPHA;
      DDIP=DIP;
% % NEED TO CHECK THE CODE CAREFULLY but here is a temporal solution!
% % high dip gives really unstable solutions
%        if DDIP>=89.999
%            DDIP = double(89.999);
%        end
      
DCCON0(AALPHA,DDIP);
% C======================================                                 05080000
% C=====  REAL-SOURCE CONTRIBUTION  =====                                 05090000
% C======================================                                 05100000
      D=DEPTH+Z;
      P=Y.*CD+D.*SD;
      Q=Y.*SD-D.*CD;
%       JXI=0;
%       JET=0;
        JXI = int8(zeros(N_CELL,1));
        JET = int8(zeros(N_CELL,1));
      aa = (X+AL1).*(X-AL2);
        cneg = aa <= 0.0;
        JXI = JXI + int8(cneg);
        jxi_sum = sum(rot90(sum(JXI)));
      bb = (P+AW1).*(P-AW2);
        cneg = bb <= 0.0; 
        JET = JET + int8(cneg);
        jet_sum = sum(rot90(sum(JET)));
%       if aa<=0.0
%           JXI=1;
%       end
%       if bb<=0.0
%           JET=1;
%       end
      DD1=DISL1;
      DD2=DISL2;
      DD3=DISL3;
%C-----                                                                  05210000
for K=1:2
      if(K==1)
          ET=P+AW1;
      end
      if(K==2)
          ET=P-AW2;
      end
    for J=1:2
        if(J==1)
            XI=X+AL1;
        end
        if(J==2)
            XI=X-AL2;
        end
            
DCCON2(XI,ET,Q,SD,CD);

% To detect singular point
    cjxi1 = JXI == 1;
    cjxi2 = JXI ~= 1;
    cjet1 = JET == 1;
    cjet2 = JET ~= 1;
    cq1   = abs(Q)   <= 1.0e-12;
    cq2   = abs(Q)   > 1.0e-12;
    cet1  = abs(ET)  <= 1.0e-12;
    cet2  = abs(ET)  > 1.0e-12;
    cxi1  = abs(XI)  <= 1.0e-12;
    cxi2  = abs(XI)  > 1.0e-12;
%     cq1   = Q  == 0.0;
%     cq2   = Q  ~= 0.0;
%     cet1  = ET  == 0.0;
%     cet2  = ET  ~= 0.0;
%     cxi1  = XI  == 0.0;
%     cxi2  = XI  ~= 0.0;
    cc1 = cjxi1.*cq1.*cet1; cc2 = (cc1 - 1.0).*(-1.0);
    cc3 = cjet1.*cq1.*cxi1; cc4 = (cc3 - 1.0).*(-1.0);
    cc0 = (cc1 + cc3) >= 1;
    cc5 = (cc1 + cc3) < 1;
	IRET = IRET + cc0;
% sum(rot90(sum(cc3)))

%         q_sum = sum(rot90(sum(Q)));
%         et_sum = sum(rot90(sum(ET)));
%         xi_sum = sum(rot90(sum(ET)));
% %        if ((jxi_sum>=1 & Q==F0 & ET==F0) | (jet_sum>=1 & Q==F0 & XI==F0))
%         if ((jxi_sum>=1 & q_sum==0.0 & et_sum==0.0) | (jet_sum>=1 & q_sum==0.0 & xi_sum==0.0))
% %        if ((jxi_sum>=1) | (jet_sum>=1))
% % C=======================================                                06030000
% % C=====  IN CASE OF SINGULAR (R=0)  =====                                06040000
% % C=======================================                                06050000
%          UX=zeros(N_CELL,1,'double');
%          UY=zeros(N_CELL,1,'double');
%          UZ=zeros(N_CELL,1,'double');
%          UXX=zeros(N_CELL,1,'double');
%          UYX=zeros(N_CELL,1,'double');
%          UZX=zeros(N_CELL,1,'double');
%          UXY=zeros(N_CELL,1,'double');
%          UYY=zeros(N_CELL,1,'double');
%          UZY=zeros(N_CELL,1,'double');
%          UXZ=zeros(N_CELL,1,'double');
%          UYZ=zeros(N_CELL,1,'double');
%          UZZ=zeros(N_CELL,1,'double');
%          IRET=ones(N_CELL,1,'double');
%             disp('error');
% %         return
%             break;
%          end
DUA = UA(XI,ET,Q,DD1,DD2,DD3);
%C-----                                                                  05320000
        for I=1:3:10
          DU(:,I)  =-DUA(:,I);
          DU(:,I+1)=-DUA(:,I+1).*CD+DUA(:,I+2).*SD;
          DU(:,I+2)=-DUA(:,I+1).*SD-DUA(:,I+2).*CD;
          if I<10.0
            continue;
          end
          DU(:,I)  =-DU(:,I);
          DU(:,I+1)=-DU(:,I+1);
          DU(:,I+2)=-DU(:,I+2);
        end
%        for I=1:1:12
          if(J+K)~=3
              U(:,1:12)=U(:,1:12)+DU(:,1:12);
          end
          if(J+K)==3
              U(:,1:12)=U(:,1:12)-DU(:,1:12);
          end
%        end
    end
end
% C=======================================                                05490000
% C=====  IMAGE-SOURCE CONTRIBUTION  =====                                05500000
% C=======================================                                05510000
      ZZ=Z;
      D=DEPTH-Z;
      P=Y.*CD+D.*SD;
      Q=Y.*SD-D.*CD;
%      JET=0;
    JET = int8(ones(N_CELL,1));
      
      cc=(P+AW1).*(P-AW2);
      
      c1 = cc <= 0.0;
      c2 = cc >  0.0;
      JET = int8(c1).*JET;
%       if cc<=0.0
%           JET=1;
%       end
%C-----                                                                  05580000
for K=1:2
      if K==1 
          ET=P+AW1;
      end
      if K==2
          ET=P-AW2;
      end
      for J=1:2
        if J==1
          XI=X+AL1;
        end
        if J==2
          XI=X-AL2;
        end
        DCCON2(XI,ET,Q,SD,CD);
        DUA = UA(XI,ET,Q,DD1,DD2,DD3);
        DUB = UB(XI,ET,Q,DD1,DD2,DD3);
        DUC = UC(XI,ET,Q,ZZ,DD1,DD2,DD3);        
%C-----                                                                  05690000
        for I=1:3:10
          DU(:,I)=DUA(:,I)+DUB(:,I)+Z.*DUC(:,I);
          DU(:,I+1)=(DUA(:,I+1)+DUB(:,I+1)+Z.*DUC(:,I+1)).*CD...
                    -(DUA(:,I+2)+DUB(:,I+2)+Z.*DUC(:,I+2)).*SD;
          DU(:,I+2)=(DUA(:,I+1)+DUB(:,I+1)-Z.*DUC(:,I+1)).*SD...
                    +(DUA(:,I+2)+DUB(:,I+2)-Z.*DUC(:,I+2)).*CD;
                if I<10.0
                    continue;
                end
          DU(:,10)=DU(:,10)+DUC(:,1);
          DU(:,11)=DU(:,11)+DUC(:,2).*CD-DUC(:,3).*SD;
          DU(:,12)=DU(:,12)-DUC(:,2).*SD-DUC(:,3).*CD;
        end
%        for I=1:1:12
          if(J+K~=3)
              U(:,1:12)=U(:,1:12)+DU(:,1:12);
          end
          if(J+K==3)
              U(:,1:12)=U(:,1:12)-DU(:,1:12);
% end
        end
%C-----                                                                  05850000
      end
end

%C=====                                                                  05880000
      UX=U(:,1);
      UY=U(:,2);
      UZ=U(:,3);
      UXX=U(:,4);
      UYX=U(:,5);
      UZX=U(:,6);
      UXY=U(:,7);
      UYY=U(:,8);
      UZY=U(:,9);
      UXZ=U(:,10);
      UYZ=U(:,11);
      UZZ=U(:,12);
      cc5 = IRET >= 1;
      IRET = cc5;
%      IRET=0;
%      IRET = zeros(N_CELL,1,'double');
%      IRET = IRET + cc0;
      sum(rot90(sum(IRET)));
%       isa(UX,'double')
%       isa(UXX,'double')
%       isa(UYZ,'double')
%       isa(IRET,'double')
%       RETURN                                                            06020000
% C=======================================                                06030000
% C=====  IN CASE OF SINGULAR (R=0)  =====                                06040000
% C=======================================                                06050000
%    99 UX=F0                                                             06060000
%       UY=F0                                                             06070000
%       UZ=F0                                                             06080000
%       UXX=F0                                                            06090000
%       UYX=F0                                                            06100000
%       UZX=F0                                                            06110000
%       UXY=F0                                                            06120000
%       UYY=F0                                                            06130000
%       UZY=F0                                                            06140000
%       UXZ=F0                                                            06150000
%       UYZ=F0                                                            06160000
%       UZZ=F0                                                            06170000
%       IRET=1                                                            06180002
%       RETURN                                                            06190000
%       END
  
% ------------------------------------------------------------------------------
function DCCON0(ALPHA,DIP)
% Okada 92 code subroutine DCCON0
%
global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D
global N_CELL
%       DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/             %09430000
%       DATA EPS/1.D-6/ 
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F2 = ones(N_CELL,1,'double').*2.0;
PI2 = ones(N_CELL,1,'double').*6.283185307179586;
EPS = ones(N_CELL,1,'double').*1.0e-6;

      ALP1=(F1-ALPHA)./F2;
      ALP2= ALPHA./F2;
      ALP3=(F1-ALPHA)./ALPHA;
      ALP4= F1-ALPHA;
      ALP5= ALPHA;

      P18=PI2./double(360.0);                                                    %09520000
      SD=double(sin(DIP.*P18));                                                  %09530000
      CD=double(cos(DIP.*P18)); 
        c1 = abs(CD) < EPS;
        c2 = abs(CD) >= EPS;
        s1 = SD > F0;
        s2 = SD == F0;
        s3 = SD < F0;

        CD = F0.*c1 + CD.*c2;

% CAUTION ************ modified by Shinji Toda (CD = 0.0 produces 'nan')
%                      in MATLAB
%         c3 = abs(CD) < EPS;
%         c4 = abs(CD) <= EPS;
%         CD = c3.*EPS + c4.*CD;
% CAUTION ***************************************************************

%09560000
%     if SD>F0
%        SD= F1;
%     end
%     if SD<F0
%         SD=-F1;                                             %09580000
%     end
    SD = c1.*(F1.*s1 + SD.*s2 + (-1.0).*F1.*s3) + c2.*SD;
%end
                                                            %09590000
      SDSD=SD.*SD;                                                        %09600000
      CDCD=CD.*CD;                                                        %09610000
      SDCD=SD.*CD;                                                        %09620000
      S2D=F2.*SDCD;                                                       %09630000
      C2D=CDCD-SDSD;                                                     %09640000
%       RETURN                                                            %09650000
%       END                                                               %09660000

% ------------------------------------------------------------------------------
function DCCON11(X,Y,D)
%       SUBROUTINE  DCCON1(X,Y,D)                                         09670000
%       IMPLICIT REAL*8 (A-H,O-Z)                                         09680000
% C                                                                       09690000
% C********************************************************************** 09700000
% C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** 09710000
% C********************************************************************** 09720000
% C                                                                       09730000
% C***** INPUT                                                            09740000
% C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    09750000
% C### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO  09760000
% C                                                                       09770000
%       COMMON /C0/DUMMY(5),SD,CD                                         09780000
%       COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     09790000
%      *           UY,VY,WY,UZ,VZ,WZ                                      09800000
global DUMMY SD CD
global P Q S T XY X2 Y2 D2 R R2 R3 R5 QR QRX A3 A5 B3 C3 UY VY WY UZ VZ WZ
global N_CELL
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F3 = ones(N_CELL,1,'double').*3.0;
F5 = ones(N_CELL,1,'double').*5.0;
EPS = ones(N_CELL,1,'double').*0.000001;

%       DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,1.D-6/                  09810000
% C-----                                                                  09820000
c1 = abs(X) < EPS;
c2 = abs(X) >= EPS;
X = F0.*c1 + X.*c2;
c1 = abs(Y) < EPS;
c2 = abs(Y) >= EPS;
Y = F0.*c1 + Y.*c2;
c1 = abs(D) < EPS;
c2 = abs(D) >= EPS;
D = F0.*c1 + D.*c2;
%       IF(DABS(X).LT.EPS) X=F0                                           09830000
%       IF(DABS(Y).LT.EPS) Y=F0                                           09840000
%       IF(DABS(D).LT.EPS) D=F0                                           09850000
	P=Y.*CD+D.*SD;
	Q=Y.*SD-D.*CD;
	S=P.*SD+Q.*CD;
	T=P.*CD-Q.*SD;
	XY=X.*Y;
	X2=X.*X;
	Y2=Y.*Y;
	D2=D.*D;
	R2=X2+Y2+D2;
	R =sqrt(R2);
%       IF(R.EQ.F0) RETURN                                                09960000
c1 = R == F0;
if sum(rot90(sum(c1))) > 0
    return
end
    R3=R .*R2;
	R5=R3.*R2;
	R7=R5.*R2;
% C-----                                                                  10000000
	A3=F1-F3.*X2./R2;
	A5=F1-F5.*X2./R2;
	B3=F1-F3.*Y2./R2;
	C3=F1-F3.*D2./R2;
% C-----                                                                  10050000
	QR=F3.*Q./R5;
	QRX=F5.*QR.*X./R2;
% C-----                                                                  10080000
	UY=SD-F5.*Y.*Q./R2;
	UZ=CD+F5.*D.*Q./R2;
	VY=S -F5.*Y.*P.*Q./R2;
	VZ=T +F5.*D.*P.*Q./R2;
	WY=UY+SD;
	WZ=UZ+CD;
%       RETURN                                                            10150000
%       END                                                               10160000

% ------------------------------------------------------------------------------
function DCCON2(XI,ET,Q,SD,CD)
% Okada 92 code subroutine DCCON2
%
global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32
global EY EZ FY FZ GY GZ HY HZ
global N_CELL

%disp('DCCON2');
% DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F2 = ones(N_CELL,1,'double').*2.0;
EPS = ones(N_CELL,1,'double').*0.000001;

c1 = abs(XI) < EPS;
c2 = abs(XI) >= EPS;
% if abs(XI)<EPS
%     XI=F0;
% end
XI = F0.*c1 + XI.*c2;
% if abs(ET)<EPS
%     ET=F0;
% end
c1 = abs(ET) < EPS;
c2 = abs(ET) >= EPS;
ET = F0.*c1 + ET.*c2;
% if abs( Q)<EPS
%     Q=F0;
% end
c1 = abs(Q) < EPS;
c2 = abs(Q) >= EPS;
Q = F0.*c1 + Q.*c2;

      XI2=XI.*XI;
      ET2=ET.*ET;
      Q2=Q.*Q;
      R2=XI2+ET2+Q2;
      R =double(sqrt(R2));
c1 = R==F0;
c1_sum = sum(rot90(sum(c1)));
if c1_sum > 0
    return;
end
      R3=R .*R2;
      R5=R3.*R2;
      Y =ET.*CD+Q.*SD;
      D =ET.*SD-Q.*CD;
%C----- 
c1 = Q == F0;
c2 = Q ~= F0;
s1 = Q.*R == F0;
s2 = Q.*R ~= F0;
% if Q==F0                                                  %10480000
%         TT=F0;                                                           %10490000
% else
% %    	if (Q.*R)==0.0                  % modified by Shinji Toda
% %        TT=double(atan(XI.*ET./EPS));  % modified by Shinji Toda 
% %        else                            % modified by Shinji Toda
%         TT=double(atan(XI.*ET./(Q.*R)));  
% %        end                             % modified by Shinji Toda
% end

% TT = c1.*F0 + c2.*(double(atan(XI.*ET./EPS)).*s1+double(atan(XI.*ET./(Q.*R))).*s2);
TT = c1.*F0 + c2.*double(atan(XI.*ET./(Q.*R)));

%C-----  
c1 = XI < F0; c2 = Q == F0; c3 = ET == F0;
c4 = c1.*c2.*c3;
c5 = zeros(N_CELL,1,'double'); c5 = (c5 - c4)+1.0;
        RXI=R+XI;
        ALX = (-double(log(R-XI))).*c4 + double(log(RXI)).*c5;
        X11 = F0.*c4 + (F1./(R.*RXI)).*c5;
        X32 = F0.*c4 + ((R+RXI).*X11.*X11./R) .*c5;
%       if(XI<F0 & Q==F0 & ET==F0)                    %10540002
%         ALX=-double(log(R-XI));                                                 %10550000
%         X11=F0;                                                          %10560000
%         X32=F0;                                                          %10570000
%       else                                                              %10580000
%         RXI=R+XI;                                                        %10590002
%         ALX=double(log(RXI));                                                   %10600000
%         X11=F1./(R.*RXI);                                                  %106%10000
%         X32=(R+RXI).*X11.*X11./R;                                           %10620002
%       end                                                             %10630000
%C----- 
c1 = ET < F0; c2 = Q == F0; c3 = XI == F0;
c4 = c1.*c2.*c3;
c5 = zeros(N_CELL,1,'double'); c5 = (c5 - c4)+1.0;
        RET=R+ET;
        ALE = (-double(log(R-ET))).*c4 + double(log(RET)).*c5;
        Y11 = F0.*c4 + (F1./(R.*RET)).*c5;
        Y32 = F0.*c4 + ((R+RET).*Y11.*Y11./R).*c5;
%         
%       if(ET<F0 & Q==F0 & XI==F0)                   %10650002
%         ALE=-double(log(R-ET));                                                 %10660000
%         Y11=F0;                                                          %10670000
%         Y32=F0;                                                          %10680000
%       else                                                              %10690000
%         RET=R+ET;                                                        %10700002
%         ALE=double(log(RET));                                                   %107%10000
%         Y11=F1./(R.*RET);                                                  %10720000
%         Y32=(R+RET).*Y11.*Y11./R;                                          %10730002
%       end                                                             %10740000
%C-----                                                                  %10750000
      EY=SD./R-Y.*Q./R3;                                                    %10760000
      EZ=CD./R+D.*Q./R3;                                                    %10770000
      FY=D./R3+XI2.*Y32.*SD;                                                %10780000
      FZ=Y./R3+XI2.*Y32.*CD;                                                %10790000
      GY=F2.*X11.*SD-Y.*Q.*X32;                                              %10800000
      GZ=F2.*X11.*CD+D.*Q.*X32;                                              %108%10000
      HY=D.*Q.*X32+XI.*Q.*Y32.*SD;                                            %10820000
      HZ=Y.*Q.*X32+XI.*Q.*Y32.*CD;                                            %10830000
%      RETURN                                                            %10840000
%      END                                                               %10850000

% ------------------------------------------------------------------------------
function [U] = UA(XI,ET,Q,DISL1,DISL2,DISL3)
%   DIMENSION U(12),DU(12)                                            06230000
% C                                                                       06240000
% C********************************************************************   06250000
% C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   06260000
% C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   06270000
% C********************************************************************   06280000
% C                                                                       06290000
% C***** INPUT                                                            06300000
% C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  06310000
% C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              06320000
% C***** OUTPUT                                                           06330000
% C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     06340000
% C                                                                       06350000
%       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  06360000
%       COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  06370000
%      *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                06380000
     
global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D
global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32
global EY EZ FY FZ GY GZ HY HZ
global N_CELL

%       DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/                     06390000
F0 = zeros(N_CELL,1,'double');
F2 = ones(N_CELL,1,'double').*2.0;
PI2 = ones(N_CELL,1,'double').*6.283185307179586;

DU = zeros(N_CELL,12,'double');
du1 = zeros(N_CELL,12,'double');
du2 = zeros(N_CELL,12,'double');
du3 = zeros(N_CELL,12,'double');

%C----- 
%for I=1:1:12
    U(1:N_CELL,1:12)=0.0;
%end
      XY=XI.*Y11;
      QX=Q .*X11;
      QY=Q .*Y11;
% C======================================                                 06460000
% C=====  STRIKE-SLIP CONTRIBUTION  =====                                 06470000
% C======================================                                 06480000
%      if DISL1~=F0
    c1 = DISL1 ~= F0;
        du1(:,1)=    TT./F2 +ALP2.*XI.*QY;
        du1(:,2)=           ALP2.*Q./R;
        du1(:,3)= ALP1.*ALE -ALP2.*Q.*QY;
        du1(:,4)=-ALP1.*QY  -ALP2.*XI2.*Q.*Y32;
        du1(:,5)=          -ALP2.*XI.*Q./R3;
        du1(:,6)= ALP1.*XY  +ALP2.*XI.*Q2.*Y32;
        du1(:,7)= ALP1.*XY.*SD        +ALP2.*XI.*FY+D./F2.*X11;
        du1(:,8)=                    ALP2.*EY;
        du1(:,9)= ALP1.*(CD./R+QY.*SD) -ALP2.*Q.*FY;
        du1(:,10)= ALP1.*XY.*CD        +ALP2.*XI.*FZ+Y./F2.*X11;
        du1(:,11)=                    ALP2.*EZ;
        du1(:,12)=-ALP1.*(SD./R-QY.*CD) -ALP2.*Q.*FZ;
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL1./PI2,1,12).*du1(1:N_CELL,1:12)...
                .*repmat(c1,1,12);
%        end
%      end
% C======================================                                 06650000
% C=====    DIP-SLIP CONTRIBUTION   =====                                 06660000
% C======================================                                 06670000
%      if DISL2~=F0
    c2 = DISL2 ~= F0;
        du2(:,1)=           ALP2.*Q./R;
        du2(:,2)=    TT./F2 +ALP2.*ET.*QX;
        du2(:,3)= ALP1.*ALX -ALP2.*Q.*QX;
        du2(:,4)=        -ALP2.*XI.*Q./R3;
        du2(:,5)= -QY./F2 -ALP2.*ET.*Q./R3;
        du2(:,6)= ALP1./R +ALP2.*Q2./R3;
        du2(:,7)=                      ALP2.*EY;
        du2(:,8)= ALP1.*D.*X11+XY./F2.*SD +ALP2.*ET.*GY;
        du2(:,9)= ALP1.*Y.*X11          -ALP2.*Q.*GY;
        du2(:,10)=                      ALP2.*EZ;
        du2(:,11)= ALP1.*Y.*X11+XY./F2.*CD +ALP2.*ET.*GZ;
        du2(:,12)=-ALP1.*D.*X11          -ALP2.*Q.*GZ;
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL2./PI2,1,12).*du2(1:N_CELL,1:12)...
                .*repmat(c2,1,12);
%      end
% C========================================                               06840000
% C=====  TENSILE-FAULT CONTRIBUTION  =====                               06850000
% C========================================                               06860000
%      if DISL3~=F0
    c3 = DISL3 ~= F0;
        du3(:,1)=-ALP1.*ALE -ALP2.*Q.*QY;
        du3(:,2)=-ALP1.*ALX -ALP2.*Q.*QX;
        du3(:,3)=    TT./F2 -ALP2.*(ET.*QX+XI.*QY);
        du3(:,4)=-ALP1.*XY  +ALP2.*XI.*Q2.*Y32;
        du3(:,5)=-ALP1./R   +ALP2.*Q2./R3;
        du3(:,6)=-ALP1.*QY  -ALP2.*Q.*Q2.*Y32;
        du3(:,7)=-ALP1.*(CD./R+QY.*SD)  -ALP2.*Q.*FY;
        du3(:,8)=-ALP1.*Y.*X11         -ALP2.*Q.*GY;
        du3(:,9)= ALP1.*(D.*X11+XY.*SD) +ALP2.*Q.*HY;
        du3(:,10)= ALP1.*(SD./R-QY.*CD)  -ALP2.*Q.*FZ;
        du3(:,11)= ALP1.*D.*X11         -ALP2.*Q.*GZ;
        du3(:,12)= ALP1.*(Y.*X11+XY.*CD) +ALP2.*Q.*HZ;
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL3./PI2,1,12).*du3(1:N_CELL,1:12)...
                .*repmat(c3,1,12);
%        end
%      end
% 
%       RETURN                                                            07030000
%       END  

% ------------------------------------------------------------------------------
function [U] = UB(XI,ET,Q,DISL1,DISL2,DISL3)
% DIMENSION U(12),DU(12)

% C                                                                       07080000
% C********************************************************************   07090000
% C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   07100000
% C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   07110000
% C********************************************************************   07120000
% C                                                                       07130000
% C***** INPUT                                                            07140000
% C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  07150000
% C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              07160000
% C***** OUTPUT                                                           07170000
% C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     07180000
% C                                                                       07190000
%       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  07200000
%       COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  07210000
%      *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                07220000
     
global ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D
global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32
global EY EZ FY FZ GY GZ HY HZ
global N_CELL

%       DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/            07230000
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F2 = ones(N_CELL,1,'double').*2.0;
PI2 = ones(N_CELL,1,'double').*6.283185307179586;

DU = zeros(N_CELL,12,'double');

%C-----                                                                  07240000
      RD=R+D;
      D11=F1./(R.*RD);
      AJ2=XI.*Y./RD.*D11;
      AJ5=-(D+Y.*Y./RD).*D11;
      
% if CD~=F0
   c1 = CD ~= F0;
   c2 = CD == F0;
   s1 = XI == F0;
   s2 = XI ~= F0;

% ----- To avoid 'Inf' and 'nan' troubles (divided by zero) ------
  tempCD = CD;
  tempCDCD = CDCD;
   CD = c1.*CD + c2.*1.0e-12;
   CDCD = c1.*CDCD + c2.*1.0e-12;
   
   X=double(sqrt(XI2+Q2));
   RD2=RD.*RD;
	AI4 = c1.*(s1.*F0 + s2.*(F1./CDCD.*( XI./RD.*SDCD...
                +F2.*atan((ET.*(X+Q.*CD)+X.*(R+X).*SD)./(XI.*(R+X).*CD)))))...
                +c2.*(XI.*Y./RD2./F2);
	AI3 = c1.*((Y.*CD./RD-ALE+SD.*double(log(RD)))./CDCD)...
            +c2.*((ET./RD+Y.*Q./RD2-ALE)./F2);
	AK1 = c1.*(XI.*(D11-Y11.*SD)./CD)+c2.*(XI.*Q./RD.*D11);
	AK3 = c1.*((Q.*Y11-Y.*D11)./CD)+c2.*(SD./RD.*(XI2.*D11-F1));
	AJ3 = c1.*((AK1-AJ2.*SD)./CD)+c2.*(-XI./RD2.*(Q2.*D11-F1./F2));
	AJ6 = c1.*((AK3-AJ5.*SD)./CD)+c2.*(-Y./RD2.*(XI2.*D11-F1./F2));

   CD = tempCD;
   CDCD = tempCDCD;
% -----
      XY=XI.*Y11;
      AI1=-XI./RD.*CD-AI4.*SD;
      AI2= double(log(RD))+AI3.*SD;
      AK2= F1./R+AK3.*SD;
      AK4= XY.*CD-AK1.*SD;
      AJ1= AJ5.*CD-AJ6.*SD;
      AJ4=-XY-AJ2.*CD+AJ3.*SD;
%C=====                                                                  07590000
%    for I=1:1:12
       U(1:N_CELL,1:12) = 0.0; 
%    end
      QX=Q.*X11;
      QY=Q.*Y11;
% C======================================                                 07640000
% C=====  STRIKE-SLIP CONTRIBUTION  =====                                 07650000
% C======================================                                 07660000
%      if DISL1~=F0
          c1 = DISL1 ~= F0;
        DU(:,1)=-XI.*QY-TT -ALP3.*AI1.*SD;
        DU(:,2)=-Q./R      +ALP3.*Y./RD.*SD;
        DU(:,3)= Q.*QY     -ALP3.*AI2.*SD;
        DU(:,4)= XI2.*Q.*Y32 -ALP3.*AJ1.*SD;
        DU(:,5)= XI.*Q./R3   -ALP3.*AJ2.*SD;
        DU(:,6)=-XI.*Q2.*Y32 -ALP3.*AJ3.*SD;
        DU(:,7)=-XI.*FY-D.*X11 +ALP3.*(XY+AJ4).*SD;
        DU(:,8)=-EY          +ALP3.*(F1./R+AJ5).*SD;
        DU(:,9)= Q.*FY        -ALP3.*(QY-AJ6).*SD;
        DU(:,10)=-XI.*FZ-Y.*X11 +ALP3.*AK1.*SD;
        DU(:,11)=-EZ          +ALP3.*Y.*D11.*SD;
        DU(:,12)= Q.*FZ        +ALP3.*AK2.*SD;
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL1./PI2,1,12).*DU(1:N_CELL,1:12)...
                .*repmat(c1,1,12);
%        end
%      end
% C======================================                                 07830000
% C=====    DIP-SLIP CONTRIBUTION   =====                                 07840000
% C======================================                                 07850000
%      if DISL2~=F0
           c2 = DISL2 ~= F0;
        DU(:,1)=-Q./R      +ALP3.*AI3.*SDCD;
        DU(:,2)=-ET.*QX-TT -ALP3.*XI./RD.*SDCD;
        DU(:,3)= Q.*QX     +ALP3.*AI4.*SDCD;
        DU(:,4)= XI.*Q./R3     +ALP3.*AJ4.*SDCD;
        DU(:,5)= ET.*Q./R3+QY  +ALP3.*AJ5.*SDCD;
        DU(:,6)=-Q2./R3       +ALP3.*AJ6.*SDCD;
        DU(:,7)=-EY          +ALP3.*AJ1.*SDCD;
        DU(:,8)=-ET.*GY-XY.*SD +ALP3.*AJ2.*SDCD;
        DU(:,9)= Q.*GY        +ALP3.*AJ3.*SDCD;
        DU(:,10)=-EZ          -ALP3.*AK3.*SDCD;
        DU(:,11)=-ET.*GZ-XY.*CD -ALP3.*XI.*D11.*SDCD;
        DU(:,12)= Q.*GZ        -ALP3.*AK4.*SDCD;
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL2./PI2,1,12).*DU(1:N_CELL,1:12)...
                .*repmat(c2,1,12);
%        end
%      end
% C========================================                               08020000
% C=====  TENSILE-FAULT CONTRIBUTION  =====                               08030000
% C========================================                               08040000
%      if DISL3~=F0
           c3 = DISL3 ~= F0;
        DU(:,1)= Q.*QY           -ALP3.*AI3.*SDSD;
        DU(:,2)= Q.*QX           +ALP3.*XI./RD.*SDSD;
        DU(:,3)= ET.*QX+XI.*QY-TT -ALP3.*AI4.*SDSD;
        DU(:,4)=-XI.*Q2.*Y32 -ALP3.*AJ4.*SDSD;
        DU(:,5)=-Q2./R3     -ALP3.*AJ5.*SDSD;
        DU(:,6)= Q.*Q2.*Y32  -ALP3.*AJ6.*SDSD;
        DU(:,7)= Q.*FY -ALP3.*AJ1.*SDSD;
        DU(:,8)= Q.*GY -ALP3.*AJ2.*SDSD;
        DU(:,9)=-Q.*HY -ALP3.*AJ3.*SDSD;
        DU(:,10)= Q.*FZ +ALP3.*AK3.*SDSD;
        DU(:,11)= Q.*GZ +ALP3.*XI.*D11.*SDSD;
        DU(:,12)=-Q.*HZ +ALP3.*AK4.*SDSD;
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL3./PI2,1,12).*DU(1:N_CELL,1:12)...
                .*repmat(c3,1,12);
%        end
%      end
%       RETURN                                                            08210000
%       END                                                               08220000

% ------------------------------------------------------------------------------
function [U] = UC(XI,ET,Q,Z,DISL1,DISL2,DISL3)
%
%      DIMENSION U(12),DU(12)                                            08250000
% C                                                                       08260000
% C********************************************************************   08270000
% C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   08280000
% C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   08290000
% C********************************************************************   08300000
% C                                                                       08310000
% C***** INPUT                                                            08320000
% C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              08330000
% C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              08340000
% C***** OUTPUT                                                           08350000
% C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     08360000
% C                                                                       08370000
%       COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  08380000
%       COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  08390000
%      *           EY,EZ,FY,FZ,GY,GZ,HY,HZ                                08400000
     
global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D
global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32
global EY EZ FY FZ GY GZ HY HZ
global N_CELL

%      DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/     08410000
%F0 = double(0.0);
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F2 = ones(N_CELL,1,'double').*2.0;
F3 = ones(N_CELL,1,'double').*3.0;
PI2 = ones(N_CELL,1,'double').*6.283185307179586;

DU = zeros(N_CELL,12,'double');

%C-----                                                                  08420000
      C=D+Z;                                                             %08430000
      X53=(double(8.0).*R2+double(9.0).*R.*XI+F3.*XI2).*X11.*X11.*X11./R2;
      Y53=(double(8.0).*R2+double(9.0).*R.*ET+F3.*ET2).*Y11.*Y11.*Y11./R2;
      H=Q.*CD-Z;
      Z32=SD./R3-H.*Y32;
      Z53=F3.*SD./R5-H.*Y53;
      Y0=Y11-XI2.*Y32;
      Z0=Z32-XI2.*Z53;
      PPY=CD./R3+Q.*Y32.*SD;
      PPZ=SD./R3-Q.*Y32.*CD;
      QQ=Z.*Y32+Z32+Z0;
      QQY=F3.*C.*D./R5-QQ.*SD;
      QQZ=F3.*C.*Y./R5-QQ.*CD+Q.*Y32;
      XY=XI.*Y11;
      QX=Q.*X11;
      QY=Q.*Y11;
      QR=F3.*Q./R5;
      CQX=C.*Q.*X53;
      CDR=(C+D)./R3;
      YY0=Y./R3-Y0.*CD;
%C=====
%    for I=1:1:12
        U(1:N_CELL,1:12)=0.0;
%    end
% C======================================                                 08660000
% C=====  STRIKE-SLIP CONTRIBUTION  =====                                 08670000
% C======================================                                 08680000
%      if DISL1~=F0
           c1 = DISL1 ~= F0;
        DU(:,1)= ALP4.*XY.*CD           -ALP5.*XI.*Q.*Z32;
        DU(:,2)= ALP4.*(CD./R+F2.*QY.*SD) -ALP5.*C.*Q./R3;
        DU(:,3)= ALP4.*QY.*CD           -ALP5.*(C.*ET./R3-Z.*Y11+XI2.*Z32);
        DU(:,4)= ALP4.*Y0.*CD                  -ALP5.*Q.*Z0;
        DU(:,5)=-ALP4.*XI.*(CD./R3+F2.*Q.*Y32.*SD) +ALP5.*C.*XI.*QR;
        DU(:,6)=-ALP4.*XI.*Q.*Y32.*CD            +ALP5.*XI.*(F3.*C.*ET./R5-QQ);
        DU(:,7)=-ALP4.*XI.*PPY.*CD    -ALP5.*XI.*QQY;
        DU(:,8)= ALP4.*F2.*(D./R3-Y0.*SD).*SD-Y./R3.*CD...
                -ALP5.*(CDR.*SD-ET./R3-C.*Y.*QR);
        DU(:,9)=-ALP4.*Q./R3+YY0.*SD  +ALP5.*(CDR.*CD+C.*D.*QR-(Y0.*CD+Q.*Z0).*SD);
        DU(:,10)= ALP4.*XI.*PPZ.*CD    -ALP5.*XI.*QQZ;
        DU(:,11)= ALP4.*F2.*(Y./R3-Y0.*CD).*SD+D./R3.*CD -ALP5.*(CDR.*CD+C.*D.*QR);
        DU(:,12)=         YY0.*CD    -ALP5.*(CDR.*SD-C.*Y.*QR-Y0.*SDSD+Q.*Z0.*CD);
 %       for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL1./PI2,1,12).*DU(1:N_CELL,1:12)...
                .*repmat(c1,1,12);
 %       end
 %     end
% C======================================                                 08860000
% C=====    DIP-SLIP CONTRIBUTION   =====                                 08870000
% C======================================                                 08880000
%      if DISL2~=F0
           c2 = DISL2 ~= F0;
        DU(:,1)= ALP4.*CD./R -QY.*SD -ALP5.*C.*Q./R3;
        DU(:,2)= ALP4.*Y.*X11       -ALP5.*C.*ET.*Q.*X32;
        DU(:,3)=     -D.*X11-XY.*SD -ALP5.*C.*(X11-Q2.*X32);
        DU(:,4)=-ALP4.*XI./R3.*CD +ALP5.*C.*XI.*QR +XI.*Q.*Y32.*SD;
        DU(:,5)=-ALP4.*Y./R3     +ALP5.*C.*ET.*QR;
        DU(:,6)=    D./R3-Y0.*SD +ALP5.*C./R3.*(F1-F3.*Q2./R2);
        DU(:,7)=-ALP4.*ET./R3+Y0.*SDSD -ALP5.*(CDR.*SD-C.*Y.*QR);
        DU(:,8)= ALP4.*(X11-Y.*Y.*X32) -ALP5.*C.*((D+F2.*Q.*CD).*X32-Y.*ET.*Q.*X53);
        DU(:,9)=  XI.*PPY.*SD+Y.*D.*X32 +ALP5.*C.*((Y+F2.*Q.*SD).*X32-Y.*Q2.*X53);
        DU(:,10)=      -Q./R3+Y0.*SDCD -ALP5.*(CDR.*CD+C.*D.*QR);
        DU(:,11)= ALP4.*Y.*D.*X32       -ALP5.*C.*((Y-F2.*Q.*SD).*X32+D.*ET.*Q.*X53);
        DU(:,12)=-XI.*PPZ.*SD+X11-D.*D.*X32-ALP5.*C.*((D-F2.*Q.*CD).*X32-D.*Q2.*X53);
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL2./PI2,1,12).*DU(1:N_CELL,1:12)...
                .*repmat(c2,1,12);
%        end
 %     end
% C========================================                               09050000
% C=====  TENSILE-FAULT CONTRIBUTION  =====                               09060000
% C========================================                               09070000
%      if DISL3~=F0
          c3 = DISL3 ~= F0;
        DU(:,1)=-ALP4.*(SD./R+QY.*CD)   -ALP5.*(Z.*Y11-Q2.*Z32);
        DU(:,2)= ALP4.*F2.*XY.*SD+D.*X11 -ALP5.*C.*(X11-Q2.*X32);
        DU(:,3)= ALP4.*(Y.*X11+XY.*CD)  +ALP5.*Q.*(C.*ET.*X32+XI.*Z32);
        DU(:,4)= ALP4.*XI./R3.*SD+XI.*Q.*Y32.*CD+ALP5.*XI.*(F3.*C.*ET./R5-F2.*Z32-Z0);
        DU(:,5)= ALP4.*F2.*Y0.*SD-D./R3 +ALP5.*C./R3.*(F1-F3.*Q2./R2);
        DU(:,6)=-ALP4.*YY0           -ALP5.*(C.*ET.*QR-Q.*Z0);
        DU(:,7)= ALP4.*(Q./R3+Y0.*SDCD)   +ALP5.*(Z./R3.*CD+C.*D.*QR-Q.*Z0.*SD);
        DU(:,8)=-ALP4.*F2.*XI.*PPY.*SD-Y.*D.*X32...
                +ALP5.*C.*((Y+F2.*Q.*SD).*X32-Y.*Q2.*X53);
        DU(:,9)=-ALP4.*(XI.*PPY.*CD-X11+Y.*Y.*X32)...
                +ALP5.*(C.*((D+F2.*Q.*CD).*X32-Y.*ET.*Q.*X53)+XI.*QQY);
        DU(:,10)=  -ET./R3+Y0.*CDCD -ALP5.*(Z./R3.*SD-C.*Y.*QR-Y0.*SDSD+Q.*Z0.*CD);
        DU(:,11)= ALP4.*F2.*XI.*PPZ.*SD-X11+D.*D.*X32...
                -ALP5.*C.*((D-F2.*Q.*CD).*X32-D.*Q2.*X53);
        DU(:,12)= ALP4.*(XI.*PPZ.*CD+Y.*D.*X32)...
                +ALP5.*(C.*((Y-F2.*Q.*SD).*X32+D.*ET.*Q.*X53)+XI.*QQZ);
%        for I=1:1:12
            U(1:N_CELL,1:12)=U(1:N_CELL,1:12)...
                +repmat(DISL3./PI2,1,12).*DU(1:N_CELL,1:12)...
                .*repmat(c3,1,12);
%        end
%      end
%       RETURN                                                            09280000
%       END                                                               09290000
