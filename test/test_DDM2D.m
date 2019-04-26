% DDM2D test
% Bench mark DDM2D to the constant pressure analytical solution.
%
clear;
close all;
source_dir = '../source';
addpath(genpath(source_dir));

mu = 10e9;                        % shear modulus
nu  = 0.25;                         % Poison's ratio.
L    = 10;                            % crack length;
N   =  100;                         % number of crack elements
p    = 1e6;       % uniform pressure 1 MPa.
flag_convergence_test = 'False';

% Numerical solution;
h = L/N;
x = [0: h : (N*h)];
% x = [h/2: h : (N*h-h/2)];
p_vec = ones(length(x),1)*p; 
[K, K_inv] = KernelBEM2D(x, mu, nu);
w_n = K_inv*p_vec; 

% Analytical solution;
c    = L/2;  
w_a = 2*c*p*(1-nu)*sqrt(1-((x-c)/c).^2)/mu;

% compare results:
figure(1);
plot(x, w_n, 'b--',x, w_a,'k-','linew',2);
legend({'DDM2D','Analytical solution'},'Location','South');
xlabel('x along the crack (m)');
ylabel('Opening (m)');
title(['Test DDM2D, N = ',num2str(N)]);
set(gca, 'fontsize',20);


    







