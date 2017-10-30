function nrm = normal(metrics,restrictions,side)

if (strcmp(side,'west')==1)
  b = -1;
  E = restrictions.W;
end
if (strcmp(side,'east')==1)
  b = 1;
  E = restrictions.E;
end
if (strcmp(side,'north')==1)
  b = 1;
  E = restrictions.N;
end
if (strcmp(side,'south')==1)
  b = -1;
  E = restrictions.S;
end

% Construct the normal based on the side
if (strcmp(side,'west')==1 || strcmp(side,'east')==1 )
  n = metrics.n2;
  gx = b*E*diag(metrics.xi_x);
  gy = b*E*diag(metrics.xi_y);
end

if (strcmp(side,'north')==1 || strcmp(side,'south')==1 )
  n = metrics.n1;
  gx = b*E*diag(metrics.eta_x);
  gy = b*E*diag(metrics.eta_y);
end

ds  = E*metrics.J*E'*sqrt(gx.^2+gy.^2);
dg =  sqrt(gx.^2+gy.^2);
nrm.n1  = gx./dg;
nrm.n2  = gy./dg;
nrm.ds  = ds;
nrm.N1  = spdiags(nrm.n1,0,n,n);
nrm.N2  = spdiags(nrm.n2,0,n,n);
nrm.DS  = spdiags(nrm.ds,0,n,n);


end
