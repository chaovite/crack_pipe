function [check, M]= pipe1d_precheck(M)
% Pre check that M valid structure for pipe1d model.
%

fprintf('------------------check pipe input----------------------\n');
fprintf('pipe name: ''%s'' \n', M.name);

check = true;

% check geometry and material properties have same dimension.
N  =  [length(M.L), length(M.R), length(M.rho), length(M.c),...
                               length(M.mu), length(M.S), length(length(M.K))];

% M.nz must have the same dimension as M.L, which is the grid points at
% each section.

if length(M.nz)~=length(M.L)
    check = false;
    warning('Must specify number of grid points in each section');
    return
end

% However, material properties and radius can be uniform, they will be
% updated in the grid.
%

for i = 2: length(N)
    if N(i) ~= length(M.L) && N(i) ~= 1
        disp(i)
        check = false;
        warning('Dimension of input doesn''t match, check number of pipe sections');
        return;
    end
end

nL = length(M.L);
fprintf('%d pipe sections defined. \n', nL);

% check number of grid points in each section.
switch M.order
    case 2
        MIN_nz = 4;
    case 4
        MIN_nz = 8;
    case 6
        MIN_nz = 12;
    case 8
        MIN_nz = 16;
    otherwise
        error('>8 higher order discretization not implemented!');
end

nz_min = min(M.nz);

if nz_min < MIN_nz
    check = false;
    warning('order %d requires minimum %d grid points in each section', M.order, MIN_nz);
    return
end

% check number interfaces.
n_int = length(M.interfaces);
if n_int ~= nL - 1
    check = false;
    warning('Number of interfaces must be 1 less than the number of pipe sections');
    return
end

fprintf('%d interfaces defined. \n', n_int);

M.crack_interfaces = {};% collect all the crack interfaces.
M.material_interfaces = {};% collect all the material interfaces.

% check each interface condition
for i = 1: n_int
    int_f    = M.interfaces{i};
    type = int_f.type;
    
    % if the interface doesn't have name give it to a default name.
    if ~isfield(int_f, 'name')
        fprintf('Interface %d miss a name, set name to ''interface %d'' \n', i, i);
        M.interfaces{i}.name = ['interface_',num2str(i)];
    end
    
    switch type
        case 'c'
            % crack coupling interface
            if ~isfield(int_f, 'link') || isempty(int_f.link)
                check = false;
                warning('No crack name for the link of crack interface %d. \n', i);
                return
            end
            M.crack_interfaces = [M.crack_interfaces; M.interfaces{i}];
            fprintf('Interface %d is linked to ''%s'' \n', i, int_f.link);
        case 'm'
            fprintf('Interface %d is a material property jump interface \n', i);
            M.material_interfaces = [M.material_interfaces; M.interfaces{i}];
        otherwise
            check = false;
            warning('Interface %d is not a valid interface type (either ''c'' or ''m'') \n');
            return
    end
end

% check boundary condition.
tp = M.bc.tp;  % top.
bt = M.bc.bt;  % bottom.

if ~isfield(M.bc.tp,'name')
    M.bc.tp.name = 'tp';
end

if ~isfield(M.bc.bt,'name')
    M.bc.bt.name = 'bt';
end

if ~isfield(tp, 'f') || ~isa(tp.f, 'function_handle')
    check = false;
    warning('Top boundary miss a function handle');
    return
end

% BT bc doesn't need a function handle, it's either p=0, v=0 or crack.
% if (~isfield(bt, 'f') || ~isa(bt.f, 'function_handle')) && ~strcmp(bt.type,'c')
%     check = false;
%     warning('Bottom boundary miss a function handle');
%     return
% end

% right now the top is only implemented with moving h boundary condition.
type = tp.type; 

switch type
    case 'h' % moving surface boundary condition.
        fprintf('Top b.c. : moving height \n');
    case 'v' % velocity bc
        fprintf('Top b.c. : velocity b.c.\n');
    case 'p' % pressure bc
        fprintf('Top b.c. : pressure b.c. \n');
    otherwise 
        check=false;
        warning('Invalid top b.c. type');
        return
end

% bottom boundary condition.
type = bt.type; 
switch type
    case 'c' % crack boundary condition.
            if ~isfield(bt, 'link') || isempty(bt.link)
                check = false;
                warning('No crack name for the link of crack bottom b.c. \n');
                return
            end
            M.crack_interfaces = [M.crack_interfaces; M.bc.bt];
            fprintf('Bottom b.c.: linked to ''%s'' \n', bt.link);
    case 'v' % velocity bc
        fprintf('Bottom b.c. : velocity b.c.\n');
    case 'p' % pressure bc
        fprintf('Bottom b.c. : pressure b.c. \n');
    otherwise
        check = false;
        warning('Invalid bottom b.c. type');
        return
end
%% print out all the crack interface
fprintf('%d material interfaces are defined:  ', length(M.material_interfaces));
for i = 1: length(M.material_interfaces)
    fprintf(' ''%s'' ', M.material_interfaces{i}.name);
end
fprintf('\n')

fprintf('%d crack interfaces are defined:  ', length(M.crack_interfaces));
for i = 1: length(M.crack_interfaces)
    fprintf(' ''%s'' ', M.crack_interfaces{i}.name);
end
fprintf('\n')
end

