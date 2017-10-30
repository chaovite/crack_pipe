function [U1, U2, U3] = disp_crack(H, R, W, strike, dip, station_loc, p, mu, nu, model)
%[U1, U2, U3] = disp_crack(H, R, W,strike, dip, station_loc, p, mu, nu, model)
%   disp_crack calculate the surface displacement from an openning crack in
%   the elastic half space. The centroid of the crack is located at [0, 0, -H]
%
% This function take p (pressure) as a time series (a row vector [1, nt])
%
%Input:
%       H               : depth of the crack center [m]
%       R               : radius of the crack [m]
%       W              : width of the crack center [m]
%       strike         : strike angle [Deg C]
%       dip             : dip angle [Deg C]
%       station_loc : location of the receiver stations [2, n_st],
%       p                : crack pressure [1, nt]
%       mu             : shear modulus [Pa]
%       nu              : Poisson's ratio
%       model        : 'Okada', 'Fialko2001', 'Sun1969'
%
%Output:
%      U1: displacement in x direction (East)
%      U2: displacement in y direction (North)  
%      U3: displacement in up direction (Up)
%
% 1. Okada        : Convert penny shape crack to rectangular crack with same moment
%                          and use Okada?s solution (disloc.m from Paul Segall?s group). 
%                          Easy to use, can have arbitrary orientation but may not be accurate.
% 
% 2. Sun1969    : Horizontal penny shape crack in elastic half space. 
%                         Easy to implement and code ready to use. But crack must be horizontal. 
%                         Accurate (2~3% off) when crack radius to depth is about 2. 
%                         This approximate solution because of ignoring the image source?s contribution to cavity.
%
% 3. Fialko2001: Horizontal penny shape crack in elastic half space. More exact solution, more accurate than Sun1969 paper. 
%                        The code is ready to use from Fialko?s website. Must be horizontal crack.
% 
% 4. Davis1986: Surface deformation due to inflation of an arbitrarily oriented triaxial ellipsoidal cavity in an elastic half-space. 
%                        Crack can be arbitrary orientation. But complicated to implement. I haven?t fully understand how to implement this. 
%                        This solution also neglect the image source?s contribution to the cavity. 
%                        So, this solution also is only accurate when the depth to radius is large enough.
%
% see also sun69.m, disloc.m, calc_crack.m

% add subfolders
addpath(genpath(pwd));

%calculate bulk modulus and Young's modulus from mu and nu
K      = 2*mu/3*(1+nu)/(1-2*nu);
E      = 2*mu*(1+nu);

% storage for U1, U2, U3
n_st = size(station_loc,2);
nt   = length(p);
x1 =  station_loc(1,:); % x1 coordinates
x2 =  station_loc(2,:); % x2 coordinates
r = sqrt(x1.^2 + x2.^2); % horizontal distance to crack


U1  = zeros(n_st,nt);
U2  = zeros(n_st,nt);
U3  = zeros(n_st,nt);

% calculate displacement according to different models

if dip~=0
    warning(['Fialko2001 and Sun1969 only account for horizontal crack, dip=0!!',...
        ' Okada''s solution handle arbitrary orientation but suffer from inaccuracy']);
end


switch model
    
    case 'Okada'
        % calculate equivalent dimension of the crack
        D     = sqrt(pi * R^2);
        
        % crack stiffness that can be obtained from Sneddon 1946
        Kc   = 3*pi*K*(1-2*nu)/(4*R/(W/2)*(1-nu^2));
        
        % equivalent openning. 
        w     = 3/4*p/Kc*(W/2);  % dimension [1 nt].
        
        % use disloc.m to calculate the surface displacement.
        
        % calculate the displacement from unit openning.
        source =  [D, D, H, strike, dip, D/2,0, 0, 0, 1]';
        U_0 = disloc(source, station_loc, nu);
            
        for i = 1: n_st;
            U1(i,:)=U_0(1,i)*w;
            U2(i,:)=U_0(2,i)*w;
            U3(i,:)=U_0(3,i)*w;
        end
        
    case 'Fialko2001'
        
        % The following parameters need to be user-supplied:   
        nis = 2;       % number of sub-intervals on [0,1] on which integration is done
                     % using a 16-point Gauss quadrature (i.e., total of nis*16 points)
        eps = 1e-6;    % solution accuracy for Fredholm integral equations (stop 
                     % iterations when relative change is less than eps)
        h = H/R;         % dimensionless crack depth (Depth/Radius ratio)

        % Solve a coupled system of Fredholm eqs. of 2nd kind for basis functions fi,psi
        % t and Wt are nodes and weights of the num. integration quadrature
         
        [fi,psi,t,Wt]=fredholm(h,nis,eps);

        [Uv,Ur]=intgr(r/R,fi,psi,h,Wt,t);  %calculate vertical and radial displ.
        Uv = - Uv;% change notation such that positive upward
        
        for i=1:n_st
            Uv_i =  Uv(i);
            Ur_i  =  Ur(i);
            const= 2*(1-2*nu)*R*p/mu; % the dimension of displacement
            
            U1(i,:) =x1(i)./r(i) *Ur_i*const;
            U2(i,:) =x2(i)./r(i) *Ur_i*const;
            U3(i,:) =Uv_i*const;
            
        end
           
    case 'Sun1969'
        
        r = sqrt(x1.^2 + x2.^2);
        
        for i = 1: n_st
            
            [Ur, Uv] = sun69(r(i),H,R,p,E,nu);

            U1(i,:) =x1(i)./r(i) *Ur;
            U2(i,:) =x2(i)./r(i) *Ur;
            U3(i,:) =Uv;
        end
                     
end


end

