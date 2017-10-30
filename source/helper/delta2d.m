function d = delta2d(x,xs,y,ys, order)
%discretization of delta function without scaled by the grid spacing.
%The sum of d is equal to 1.
d_x = delta1d(x,xs, order);
d_y = delta1d(y,ys, order);
d = kron(d_x, d_y);
end
