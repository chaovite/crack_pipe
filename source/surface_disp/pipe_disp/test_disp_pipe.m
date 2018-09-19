% Check the validity of disp_pipe.m
%
% Set constant pressure from c_1 (depth 1) to c_2 (depth 2) and test if my code disp_pipe.m gives the same
% solution with the integration calculated analytically in the paper. 
%
% Refer to Bonaccorso and Davis 1999 for detailed information
%
% Note that in Bonaccorso and Davis 1999, positive vertical displacement
% point downward while here we adopt the notation such that positive
% displacement point upward.

model='close_pipe'; % model 'open_pipe' or 'close pipe'
N_st=1000; % number of station
N_z=100; % number of grid points in z direction
R=5; % conduit radius
c1 = R; % depth 1
c2 = 20*R; % depth 2
P=1e4; % magnitude pressure
pipe_loc=[0 0]'; % conduit location
x=zeros(1,N_st); % x coordinates (East)
y = linspace(0.1*R,10*R,N_st); %y coordinates (North)
station_loc = [x; y]; % station location
p=P*ones(1,N_z)'; % pressure along the depth
z = linspace(c1,c2,N_z)';% depth of each grid point
mu = 1e10; % shear modulus
nu = 0.25; % Poisson's ratio

%%  calculate solution from numeric integration
[U, T]= disp_pipe( pipe_loc, station_loc, p, z, R, mu, nu, model);
%% calculate the solution from analytic integration from equation given by 

% close_pipe
a = sqrt(x.^2+y.^2); r1 = sqrt(x.^2+y.^2+c1^2); r2 = sqrt(x.^2+y.^2+c2^2);

switch model
    case 'close_pipe'
        %equation (7a) and (7b), page 10,534
        m  = 3/4*P/mu*R^2; n = 1/3;
        U1 = m*x.*(2*(1-2*nu)*(c2./(a.^2.*r2)-c1./(a.^2.*r1))...
                 -n*((c2^3./(a.^2.*r2.^3)-2*nu*c2./(a.^2.*r2))-(c1^3./(a.^2.*r1.^3)-2*nu*c1./(a.^2.*r1))));

        U2 = m*y.*(2*(1-2*nu)*(c2./(a.^2.*r2)-c1./(a.^2.*r1))...
                 -n*((c2^3./(a.^2.*r2.^3)-2*nu*c2./(a.^2.*r2))-(c1^3./(a.^2.*r1.^3)-2*nu*c1./(a.^2.*r1))));

        U3 =- m* (   2*(1-2*nu).*(1./r2-1./r1) ... 
                           -n*(3-2*nu).*(1./r2-1./r1) + n * (a.^2./r2.^3 - a.^2./r1.^3)   );
                       
        % tilt
        T1 = -m*x.*(1./r2.^3.*(n*(5-2*nu) -2 + 4*nu - 3*n*a.^2./r2.^2)  -  1./r1.^3.*(n*(5-2*nu) -2 + 4*nu - 3*n*a.^2./r1.^2));
        T2 = -m*y.*(1./r2.^3.*(n*(5-2*nu) -2 + 4*nu - 3*n*a.^2./r2.^2)  -  1./r1.^3.*(n*(5-2*nu) -2 + 4*nu - 3*n*a.^2./r1.^2));
    case 'open_pipe'
        %equation (11a) and (11b), page 10,535
        m =  P/mu*R^2;
        U1 = m*x.*((1+nu)*(c2./(a.^2.*r2)-c1./(a.^2.*r1))...
                 -(1/2*c2^3./(a.^2.*r2.^3)-1/2*c1^3./(a.^2.*r1.^3)));
             
        U2 = m*y.*((1+nu)*(c2./(a.^2.*r2)-c1./(a.^2.*r1))...
                 - (1/2*c2^3./(a.^2.*r2.^3)-1/2*c1^3./(a.^2.*r1.^3)));
             
        U3 = -m* (   (2*nu-1)/2.*(1./r2-1./r1) ... 
                           + 1/2 * (a.^2./r2.^3 - a.^2./r1.^3)   );   
                       
       % tilt
       T1  = -m * x.* (1./r2.^3.*(3/2-nu - 3*a.^2./2./r2.^2) - 1./r1.^3.*(3/2-nu - 3*a.^2./2./r1.^2));
       T2  = -m * y.* (1./r2.^3.*(3/2-nu - 3*a.^2./2./r2.^2) - 1./r1.^3.*(3/2-nu - 3*a.^2./2./r1.^2));
end

%% compare results from numeric integration to analytic results to make sure the implementation is corrrect
figure(1);
plot(y/R,U(2,:)*1e6,'r',y/R,U(3,:)*1e6,'b',y/R,U2*1e6,'r*',y/R,U3*1e6,'b*');
%  ylim([-4,4])
h = legend({'radial numeric','up numeric','radia analytic','up analytic'});
set(h,'fontsize',18);
xlabel('distance / R','fontsize',18);
ylabel('Displacement (\mum)','fontsize',18);
set(gca,'fontsize',18);
shg

figure(2);
plot(y/R,T(2,:)*1e6,'r', y/R, T2*1e6,'k.');
xlabel('distance / R','fontsize',18);
ylabel('Tilt (micron)','fontsize',18);
set(gca,'fontsize',18);
