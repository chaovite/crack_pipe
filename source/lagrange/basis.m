function lj = basis(j,x,xs,lambda)
  if abs(xs-x(j)) < eps
    lj = 1;
    return
  end
  lj = l_fun(x,xs)*lambda(j)/(xs-x(j));

end
