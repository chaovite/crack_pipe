function [y_stretch] = boundary_layer_thickness(y_unit, r_g, r_bl)
% [y_stretch] = boundary_layer_thickness(y_unit, r_g, r_bl)
% Returns a stretched grid given the uniform unit grid y_unit (from 0 to 1),  the fraction of grid points r_g in the
% boundary layer and the ratio of the boundary layer thickness to
% crack width r_bl.

y = @(eta,r) 0.5 + 0.5*tanh(r*(2*eta-1))/tanh(r);
f = @(r) r_bl - y(r_g,r); 
r = fsolve(f,2,optimoptions('fsolve','Display','off'));
if r_g~=r_bl
    y_stretch = y(y_unit, r);
else 
    % no stretch
    y_stretch = y_unit;
end
end

    

