function [U, T] = disp_pipe( pipe_loc, station_loc, p, z, R, mu, nu, model)
%[U, T]= disp_pipe( pipe_loc, station_loc, p, z, R, mu, nu)
%   [U, T] = disp_pipe( pipe_loc, station_loc, p, z, R, mu, nu) calculates the
%   surface displacement U(1), U(2), U(3) (East, North, Up) from
%   a pressurized vertical conduit or pipe with depth varying pressure
%   distribution. T(1), T(2) are the tilt components, east, north.
%   The receiver locations are defined in station_loc.
%  
%
%Input:
%     pipe_loc      : vector of dimension [2,1], horizontal coordinates: East,
%                           North specifying the location of the conduit
%     station_loc  : matrix of dimension [2,N], specify the (x, y)
%                          coordinates of N receiver stations.
%     p                 : vector of dimension [nz+1,1], pressure distribution along the depth 
%     z                 : depth [nz+1,1], positive downward 
%     R                : radius of the conduit [nz+1, 1] or a scalar.
%     mu              : shear modulus of the wall rock
%     nu               : Poisson's ratio of the wall rock
%     model         : 'open_pipe' or 'close_pipe'. ('open_pipe' by default). Based on the solution from Bonaccorso and Davis
%
%Output:
%     U                : surface displacement U(1), U(2), U(3) (East, North, Up), Dimension [3,N]
%     T                : surface tilt T(1), t(2) (East, North), Dimension [2,N]
%
%
% This function is based on the paper Bonaccorso and Davis (1999):
%
%       Models of ground deformation from vertical volcanic conduits with
%       application to eruptions of Mount St. Helens and Mount Etna. 
%   
%       JGR, VOL 104, NO.B5, PAGES 10,531-10,542, MAY 10, 1990
%
%Implemented by: CHAO LIANG, NOV. 24, 2015,@ Stanford University
%Modified by: CHAO LIANG, SEP. 19, 2018,@Stanford University, implementing
%tilt.
%
%Example: Please run test_disp_pipe.m to see the comparison of numeric
%                integration and analytic solution.
% 
% See also: test_disp_pipe.m

%check input:
if size(pipe_loc,1) ~=2;
    error('Dimension of pipe_loc must be [2,1]');
end

if size(station_loc,1) ~=2;
    error('Dimension of station_loc must be [2,N]');
end

if size(p,2) ~=1 || size(z,2) ~=1;
    error('Please set p and z to be a column vector');
end

%default model is set to be 'open_pipe'
if nargin<8
    model = 'open_pipe';
end

%
N  = size(station_loc,2); % number of stations;

nz = length(p) - 1; % number of elements in z direction

dz = z(2:end) - z(1:end-1);% grid spacing. It's possible that the grid are not equally spaced.

z_c = 1/2 * (z(2:end) + z(1:end-1));% depth at the central points of two grid points

% calculate the pressure and radius at the central points between grid points.
p_c = 1/2 * (p( 2 : end) + p(1 : end - 1));

if length(R) == 1;%scalar
    R_c = R*ones(nz,1);
else
    R_c = 1/2 * (R( 2 : end) + R(1 : end - 1));
end 

% put conduit at the origin and calculate the coordinates of receiver locations
x = station_loc;
x(1,:) = x(1,:) - pipe_loc(1);% Dimension [1 N]
x(2,:) = x(2,:) - pipe_loc(2);% Dimension [1 N]

U = zeros(3,N); % displacements
T = zeros(2,N); % tilts.

% loop over stations
for i=1:N;
     x_i = x(:,i); %coordinate of the i-th station [2, 1].
     r_i  = sqrt(x_i(1).^2 + x_i(2).^2 + z_c.^2);% distance from the station to source
     
     % loop over displacement components
     for j = 1:3
         
         if j<3 % for U1 and U2
             % equation (6a), page 10534
             U11_U22_U33_j = (1 - 2*nu)/(2*pi*mu)*x_i(j)./(r_i.^3);%U1,1^j+U2,2^j+U3,3^j
             % equation (6c), page 10534
             U33_j = x_i(j)./(4*pi*mu).*(3*z_c.^2./r_i.^5- 2*nu./r_i.^3);%U33^j.
         else  % for U3
             % equation (6b), page 10534
             U11_U22_U33_j =   -(1 - 2*nu)/(2*pi*mu).*z_c./r_i.^3;
             % equation (6d), page 10534
             U33_j = z_c./(4*pi*mu).*(-3*z_c.^2./r_i.^5+2*nu./r_i.^3);
         end
             
             switch model
                 
                 case 'open_pipe'
                    % equation (10c), page 10535
                    U(j,i) = 2*pi*sum(R_c.^2.*p_c.*(1/(1-2*nu)*U11_U22_U33_j-U33_j).*dz);
                    
                    U(3,i)= -U(3,i);% different sign notation from the paper for vertical component
                    
                 case 'close_pipe'
                    % equation (5b), page 10534
                     U(j,i) = 3*pi*sum(R_c.^2.*p_c.*(U11_U22_U33_j-1/3*U33_j).*dz);
                     
                     U(3,i)= -U(3,i);% different sign notation from the paper for vertical component
             end
     end
     
     % calculate tilt.
     for j = 1: 2
         d_ujj3_dxi    = 3*(1-2*nu)/(2*pi*mu)*z_c./r_i.^4*x_i(j)./r_i;
         d_u333_dxi  = z_c./(4*pi*mu).*(5*3*z_c.^2./(r_i.^7) - 3*2*nu./(r_i.^5)).*x_i(j);
         switch model
             case 'open_pipe'
                 T(j,i)  = 2*pi*sum(R_c.^2.*p_c.*(1/(1-2*nu)*d_ujj3_dxi - d_u333_dxi ).*dz);
             case 'close_pipe'
                 T(j,i)  = 3*pi*sum(R_c.^2.*p_c.*(d_ujj3_dxi - 1/3*d_u333_dxi).*dz);
         end
         T(j,i) = -T(j,i); % different sign notation for vertical component.
     end
end

end

