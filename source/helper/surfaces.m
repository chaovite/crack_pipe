function surf = surfaces(metrics,op)

  % Surface terms
  nrm                    = normal(metrics,op.restrictions,'north');
  surf.normal.north      = nrm;
  surf.quadrature.north  = nrm.DS*op.Px;

  nrm                    = normal(metrics,op.restrictions,'south');
  surf.normal.south      = nrm;
  surf.quadrature.south  = nrm.DS*op.Px;
  
  nrm                    = normal(metrics,op.restrictions,'west');
  surf.normal.west       = nrm;
  surf.quadrature.west   = nrm.DS*op.Py;
  
  nrm                    = normal(metrics,op.restrictions,'east');
  surf.normal.east       = nrm;
  surf.quadrature.east   = nrm.DS*op.Py;

end
