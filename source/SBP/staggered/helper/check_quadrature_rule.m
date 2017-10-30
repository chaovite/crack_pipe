function pass = check_quadrature_rule(Pp,Pm,order,h,verbose)
% pass = check_quadrature_rule(Pp,Pm,order,h,verbose)
% Input arguments:
% Pp and Pm : SBP staggered grid quadrature rule matrices
% order     : Order of accuracy at the boundary of the SBP operator
% h         : Grid spacing
% verbose   : Show output (optional)
%
% Output arguments:
% pass      : Return true if quadrature rule is 2*s-1 order accurate (s = order)

  pass = false;

  if nargin < 5
    verbose = true;
  end

  n = size(Pp,1)-1;
  xp = h*[0:n]';
  xm = h*[0 1/2+0:n n]';  

  if verbose
    fprintf(' Test polynomial to integrate : x^(2*s-1), exact integral I = x(end)^(2s)/2s - x(1)^(2s)/2s \n');
    fprintf('   \t s  \t error  \n');
  end
  for j=0:order
      a = sum(2*j*Pp*xp.^max([2*j-1,0]));
      b = sum(2*j*Pm*xm.^max([2*j-1,0]));
      err_p = abs(a-(xp(end).^max([2*j,0])-xp(1).^max([2*j,0])));
      err_m = abs(b-(xm(end).^max([2*j,0])-xm(1).^max([2*j,0])));
      if verbose
        fprintf('Pp \t %d  \t %e \n',j,err_p);
        fprintf('Pm \t %d  \t %e \n',j,err_m);
      end
  end

  tol  = eps;
  err  = @(e) max(abs(e));
  pass = all(and(err(err_p)<eps,err(err_m)<eps));


end
