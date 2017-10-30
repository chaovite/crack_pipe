function d = delta1d(x,xs, order)
%discretization of delta function without scaled by the grid spacing.
%The sum of d is equal to 1.
d = spalloc(length(x),1,1);
[~, index] = min(abs(x-xs));
range = index  + (-order/2:order/2);
max_range = min(max(range), length(x));
min_range  = max(min(range), 1);
range = min_range: max_range;
l = interp(x(range), xs);
d(index) = 1;
d(range) = double(l);
end

