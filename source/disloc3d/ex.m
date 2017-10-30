function ex()
% Example showing how to use disloc3d.
  mu = 1;
  nu = 0.25;
  n = 200;
  d = -100;
  x = linspace(-3,1,n);
  z = linspace(d-2,d+2,n);
  [xm zm] = meshgrid(x,z);
  obs = [xm(:)'
	 zeros(1,n^2)
	 zm(:)'];
  % I sent length to a very large number to simulate a 2D (in-plane) problem.
  length = 1e5; % N-S
  width = 2;    % E-W
  depth = -d;
  dip = 12;
  strike = 0;
  east = 0;
  north = 0;
  ss = 0;
  ds = 1;
  op = 0;
  mdl = [length width depth dip strike east north ss ds op]';

  % You can use any of these three versions. The 'pm' version is very slow.
  tic
  [U D S flag] = disloc3d(mdl,obs,mu,nu);
%   [U D S flag] = disloc3dpm(mdl,obs,mu,nu);
%   [U D S flag] = disloc3domp(mdl,obs,mu,nu,4);
  toc
  
  % Displacements
  figure; clf;
  c = 'xyz';
  for(i=1:3)
    subplot(2,2,i); im(x,z,reshape(U(i,:),n,n));
    title(sprintf('U_{%s}',c(i)));
  end
  figure; clf;
  Ux = reshape(U(1,:),n,n);
  Uz = reshape(U(3,:),n,n);
  k = ceil(n/20);
  xs = x(1:k:end); zs = z(1:k:end);
  Ux = Ux(1:k:end,1:k:end); Uz = Uz(1:k:end,1:k:end);
  quiver(xs,zs,Ux,Uz); axis tight;
  title('Displacements');
  
  % Stresses in global coordinates
  figure; clf;
  subplot(221); im(x,z,reshape(S(1,:),n,n));
  title('S_{xx}');
  subplot(222); im(x,z,reshape(S(3,:),n,n));
  title('S_{xz}');
  subplot(223); im(x,z,reshape(S(6,:),n,n));
  title('S_{zz}');
  
  % Stress resolved along and normal to the fault
  cd = cosd(dip);
  sd = sind(dip);
  normal = [sd 0 cd]';
  along = [cd 0 -sd]';
  [shr_stress nml_stress] = Stresses(S,along,normal);
  figure;
  subplot(211); im(x,z,reshape(shr_stress,n,n)); title('Shear stress');
  subplot(212); im(x,z,reshape(nml_stress,n,n)); title('Normal stress');
  
function im(x,y,I)
  n = length(x);
  imagesc(x,y,reshape(I,n,n));
  colorbar;
  caxis(1e1*median(abs(I(:)))*[-1 1]);
  axis xy;

function [ss ns] = Stresses(S,along,normal)
% From S, the output of disloc3d, compute the shear and normal stresses. along
% and normal are vectors pointing along and normal to the fault.
  n = size(S,2);
  ss = zeros(n,1);
  ns = zeros(n,1);
  % Vectorized code that is equivalent to
  %   sigma = [Sxx(i) Sxy(i) Sxz(i)
  %            Sxy(i) Syy(i) Syz(i)
  %            Sxz(i) Syz(i) Szz(i)];
  %   ss(i) = along' *sigma*normal;
  %   ns(i) = normal'*sigma*normal;
  normal = normal(:);
  S = S([1 2 3 2 4 5 3 5 6],:);
  nml = repmat(normal(repmat(1:3,1,3)),1,n);
  a = S.*nml;
  a = [sum(a(1:3,:))
       sum(a(4:6,:))
       sum(a(7:9,:))];
  ss = along(:)'*a;
  ns = normal(:)'*a;
