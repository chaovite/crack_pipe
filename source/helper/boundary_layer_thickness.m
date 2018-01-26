function [y_stretch, J] = boundary_layer_thickness(y_unit, r_g, r_bl)
% [y_stretch, J] = boundary_layer_thickness(y_unit, r_g, r_bl)
% Returns a stretched grid given the uniform unit grid y_unit (from 0 to 1),  the fraction of grid points r_g in the
% boundary layer and the ratio of the boundary layer thickness to
% crack width r_bl.
% y_stretch: the stretched grid
% J             : the Jacobian, dz/deta (z is stretched and eta is standard grid)


y = @(eta,r) 0.5 + 0.5*tanh(r*(2*eta-1))/tanh(r);
f = @(r) r_bl - y(r_g,r); 
r = fsolve(f,2,optimoptions('fsolve','Display','off'));
n = length(y_unit);
if r_g~=r_bl
    y_stretch = y(y_unit, r);
    J             =  1.0*r*(-tanh(r*(2*y_unit - 1)).^2 + 1)/tanh(r);
    J             =  spdiags(J', 0, n, n); 
else 
    % no stretch
    y_stretch = y_unit;
    J              = speye(n);
end

end

    

