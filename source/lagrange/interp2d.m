function l = interp2d(x,y,xs,ys)
  % Returns the Lagrange basis functions l_j(xs) that constructs an interpolant and evaluates at it at (xs,ys):
  % f(xs,ys) = sum ( f_jl_j(xs,ys) ) 
  lambda_x = lambda_coeffs(x);
  lambda_y = lambda_coeffs(y);
  n = length(x);
  m = length(y);
  l = zeros(n,m);
  for j=1:n
    for k=1:m
      l(k,j) = basis(j,x,xs,lambda_x)*basis(k,y,ys,lambda_y);
    end
  end

end
