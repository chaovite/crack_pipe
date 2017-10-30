function l = interp(x,xs)
  % Returns the Lagrange basis functions l_j(xs) that constructs an interpolant and evaluates at it at xs:
  % f(xs) = sum ( f_jl_j(xs) ) 
  lambda = lambda_coeffs(x);
  n = length(x);
  l = zeros(n,1);
  for j=1:n
    l(j) = basis(j,x,xs,lambda);
  end

end
