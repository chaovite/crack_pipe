function [edge2el, el2edge, edgeBC, transT, edge_len, connDist, TR] = LoadGmsh2D(file, isplot)
%    [edge2el, el2edge, edgeBC, transT, edge_len, connDist, TR] = LoadGmsh2D(file, isplot)
%
%------------------------------------------------------------------------%
%------ Gmsh to Matlab script: Import 2d triangular mesh to matlab---------------------%
%------------------------------------------------------------------------%
%

if nargin<2
    isplot=false;
end

%-----------------------------------------------------------------------%
% dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range 
% bounded by row offsets R1 and R2 and column offsets C1 and C2.
%-----------------------------------------------------------------------%
%

% no of nodes is mentioned in 5th row and first column
N_n      = dlmread(file,'',[5-1 1-1 5-1 1-1]);
N_e      = dlmread(file,'',[7+N_n 0 7+N_n 0]);

nodes       = dlmread(file,'',[5 1 4+N_n 3]);
elements    = dlmread(file,'',[8+N_n 0 7+N_n+N_e 7]);

%-------- 2D Geometry
nodes2d = nodes(:,1:2);
elem_type   = elements(:,2);

%------- 2D Elements
elements2d = elements(elem_type==2, 6:8);
elements2d = sort(elements2d, 2);

% create triangulation object
TR = triangulation(elements2d, nodes2d);

% centers of triangles.
C = incenter(TR); 

% obtain all the edges
edges2d = edges(TR);
edges2d = sort(edges2d, 2);

n_edges = size(edges2d, 1);

xy1 = nodes2d(edges2d(:, 1), :);
xy2 = nodes2d(edges2d(:, 2), :);

% edge lengths
edge_len = sqrt(sum((xy2 - xy1).^2, 2));

% edge2el, connection list
edge2el = edgeAttachments(TR, edges2d);

edgeBC  = zeros(n_edges, 1);  % boundary edge markers
connDist = nan(n_edges, 1);     % connection distance

for i = 1: n_edges
    els = edge2el{i};
    
    if length(els)==1
        %boundary edge
        edgeBC(i) = 1;
        continue
    end
    
    c1 = C(els(1), :);
    c2 = C(els(2), :);
    
    connDist(i) =  norm(c1 - c2);
end

edgeBC = logical(edgeBC);

transT = edge_len./connDist;

%------- el2edge ------

%------- hash edges into integer number
hashBits   = 32;
edgemap  = cell(2^16, 1);

edgecode = mod(edges2d(:,1)*(N_n + 1) +  edges2d(:,2), 2^hashBits) + 1;
% check this edgecode is unique;
assert(length(edgecode)==length(unique(edgecode)));

edgemap(edgecode) = num2cell([1:n_edges]');

N_e2d = size(elements2d, 1);
el2edge = cell(N_e2d, 1);

for i  = 1: N_e2d 
    edge_index = zeros(1, 3);
    cnt = 0;
    
    % 3 edges for triangular elements
    % (1, 2), (1, 3), (2, 3).
    
    for k = 1: 2
        for m = k+1: 3
                cnt = cnt + 1;
                n1 = elements2d(i, k); 
                n2 = elements2d(i, m); 

                c12 =  mod(n1*(N_n + 1) +  n2, 2^hashBits) + 1;% hashcode
                edge_index(cnt) = edgemap{c12};
        end
        
    end
    
    el2edge{i} = edge_index;
    
end

if isplot
%---- visualize in matlab ---------------------
figure(1)
triplot(elements2d, nodes2d(:,1),nodes2d(:,2))
hold on;
plot(C(:,1), C(:,2),'r.'); % cell centers

% plot boundary edges
edges_bc= edges2d(edgeBC, :);
xs = nodes2d(edges_bc', 1);
ys = nodes2d(edges_bc', 2);
plot(xs, ys,'k*'); % cell centers
hold off;
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
legend({'elements','centers','boundary'});
title('GMsh to MATLAB import','fontsize',14);
fh = figure(1);
set(fh, 'color', 'white'); 
daspect([1,1,1]);
end

%-------------------------------------------------------------------------

end

