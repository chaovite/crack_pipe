    function indices = field_indices(dim,num)
      scan = 1;
      for i=1:num
        scan(i+1) = scan(i) + dim(i);
      end
      indices = scan(num):(scan(num+1)-1);
    end
