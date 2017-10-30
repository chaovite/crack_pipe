function [hx hy] = compute_grid_spacing(metric)
% h_min = compute_grid_spacing(metric)
% 
% Computes the grid spacing in a curvilinear grid

  grad_xi  = min(sqrt(diag(metric.xi_x).^2  + diag(metric.xi_y).^2));
  grad_eta = min(sqrt(diag(metric.eta_x).^2 + diag(metric.eta_y).^2));
  hx       = metric.h1*diag(metric.J).*grad_xi;
  hy       = metric.h2*diag(metric.J).*grad_eta;


end
