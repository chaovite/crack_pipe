function [xp xm] = grids_1d_periodic(n)
  % Grids
  h = 1/n;
  xp = [0:n-1]'*h;
  xm = ([0:n-1]-1/2)'*h;
end
