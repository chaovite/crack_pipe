function g = curved_grids(grid_type,file,nx,ny,operator_type)


  if nargin < 4
    operator_type = 'weak';
  end

  grid_types = select_grid(grid_type);


  cad_file     = iges2matlab(file);
  nurbs        = get_nurbs(cad_file);
  connectivity = nurbs_to_mesh(nurbs);
  
  % Temporary hack to resolve the rotation of blocks
  connectivity = fix_rotations(connectivity);
  
  % Apply transfinite interpolation to each element in the mesh
  msh = nurbs_tfi(nurbs,connectivity.E2e,2*nx+1,2*ny+1);
  % TODO: Add support for multi-block
  msh = msh{1};
  
  x = linspace(-10,10,2*nx+1);
  y = linspace(4,24,2*ny+1);
  x = linspace(0,1,2*nx+1);
  y = linspace(0,1,2*ny+1);
  %[msh.X,msh.Y] = meshgrid(x,y);
  
  pid = @(n) 1:2:(2*n+1);
  mid = @(n) [1 (2:2:2*n) 2*n+1];
  
  for i=1:length(msh)
    g.nx = nx;
    g.ny = ny;
    g.hx = 1/nx;
    g.hy = 1/ny;
    %g.hx = x(2) - x(1);
    %g.hy = y(2) - y(1);

    msh.X = msh.X;
    msh.Y = msh.Y;

  if grid_types.is_pp
    g.X       = msh.X(pid(ny),pid(nx));
    g.Y       = msh.Y(pid(ny),pid(nx));
  end
  
  if grid_types.is_mm
    g.X       = msh.X(mid(ny),mid(nx));
    g.Y       = msh.Y(mid(ny),mid(nx));
  end
  
  if grid_types.is_pm
    g.X       = msh.X(mid(ny),pid(nx));
    g.Y       = msh.Y(mid(ny),pid(nx));
  end
  
  if grid_types.is_mp
    g.X       = msh.X(pid(ny),mid(nx));
    g.Y       = msh.Y(pid(ny),mid(nx));
  end

  g.n1 = size(g.X,2);
  g.n2 = size(g.X,1);

  g.vec = @(u) reshape(u,g.n1*g.n2,1);
  g.grd = @(u) reshape(u,g.n2,g.n1);
  
  end
