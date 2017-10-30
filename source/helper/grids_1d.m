function [xp xm] = grids_1d(n,operator_type)
  if nargin < 2
    operator_type = 'weak';
  end
  %n = n -1;
  h = 1/n;
  % Grids
  xp = [0:n]'*h;
  xm = ([0,1/2+(0:n-1),n])'*h;

  if strcmp(operator_type,'strong')==1
    xm = xm(2:end-1);
  end
end
