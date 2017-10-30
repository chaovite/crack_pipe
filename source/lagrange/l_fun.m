function l = l_fun(x,xs)
  l = 1;
  for k=1:length(x)
    l = l*(xs-x(k));
  end
end
