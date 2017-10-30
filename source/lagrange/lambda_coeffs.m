function lambda = lambda_coeffs(x)
  n = length(x);
  lambda = ones(n,1);
  for j=1:n
    for k=1:n
      if j == k
        continue
      end
      lambda(j) = lambda(j)*(x(j)-x(k));
    end
  end
  lambda = 1./lambda;

end
