function [data A] = load_transfer_function_data(file,root)
% data = load_transfer_function_data(file,root)
%
% Input:
% file: FDMAP project to load
% root: path to FDMAP project
%
% Output:
% data: struct that contains t (time vector), u (width averaged velocity at
% fracture mouth), and p
% (pressure at fracture mouth) as well as Z (fluid impedance, Z = rho*cp)

A = init(file,root,'l',8);

data.name = file;
u = loadfast(A,'I_u');
p = loadfast(A,'I_p',size(u.u,2),1);
data.t = u.t;
data.u = u.u(1,:);
data.p = p.p(1,:);  
data.Z = A.HF{6}.rho0*A.HF{6}.c0;

end
