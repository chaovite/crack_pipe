function v = zpad(u,n)

  if length(u) >= n
      v = u;
      warning('nothing to pad');
      return
  end
  
  z = zeros(n - length(u),1);
  if size(u,1) > 1
    v = [u;z];
  else
    v = [u,z'];
  end
end
