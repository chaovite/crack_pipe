function grid_types = select_grid(grid_type)

  % grid types
  grid_types.is_pp = false;
  grid_types.is_mm = false;
  grid_types.is_pm = false;
  grid_types.is_mp = false;
  grid_types.is_p  = false;
  grid_types.is_m  = false;

  switch grid_type
    case 'pp'
      grid_types.is_pp = true;
    case 'mm'
      grid_types.is_mm = true;
    case 'pm'
      grid_types.is_pm = true;
    case 'mp'
      grid_types.is_mp = true;
    case 'p'
      grid_types.is_p  = true;
    case 'm'
      grid_types.is_m  = true;
  end
