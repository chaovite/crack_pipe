function [edge2el, el2edge, edgeBC, transT, areaEl, Points, el2pt, Cs, dX, dY] = rectMesh2d(x, y, isplot)
%
%   [edge2el, el2edge, edgeBC, transT, areaEl] = rectMesh2d(x, y)
%
%   Generate 2D finite volume grid for a rectangular mesh.
%
%  x: x coordinates of vertices in x direction
%  y: y coordinates of vertices in y direction
%

tic;
disp('Creating mesh...');
% number of vertices
nvx = length(x);
nvy = length(y);

% number of cells.
ncx = nvx-1;
ncy = nvy-1;

N_edges = ncx*nvy + ncy*nvx;
N_el = ncx*ncy;

edge2el = cell(N_edges, 1);
el2edge = cell(N_el, 1);
edgeBC  = zeros(N_edges, 1);
transT    = nan(N_edges, 1);
areaEl    = zeros(N_el, 1);
el2pt      = zeros(N_el, 4);

[X, Y]     = meshgrid(x, y);
Points    = [X(:), Y(:)];

edge2pt = zeros(N_edges, 2);
Cs          = zeros(N_el, 2);
dX          = zeros(N_el, 1);
dY          = zeros(N_el, 1);

for i = 1: ncx
    for j = 1: ncy
        % element index
        ind_el = (i-1) * ncy + j;
        
        Cs(ind_el, :) = [(x(i)+x(i+1))/2, (y(j)+y(j+1))/2];
        
        dx_i = x(i+1) - x(i);
        dy_j = y(j+1) - y(j);
        dX(ind_el) = dx_i;
        dY(ind_el) = dy_j;
        
        % element area
        areaEl(ind_el) = dx_i*dy_j;
        
        % 4 bounding vertices
        pt_lb = (i-1)*nvy + j;%left-bottom
        pt_lt = pt_lb + 1;%left-top
        pt_rb = pt_lb + nvy;%right-bottom
        pt_rt = pt_rb + 1;%right-top
        
        % el2pt
        el2pt(ind_el, :) = [pt_lb, pt_rb, pt_rt, pt_lt];
        
        % el2edge
        edge_l   = ind_el; %left
        edge_r  = ind_el + ncy;%right
        edge_b = (i-1)*nvy + j + nvx*ncy;%bottom
        edge_t  = edge_b + 1;%top
        
        edge_index = [edge_l, edge_r, edge_b, edge_t];
        el2edge{ind_el} = edge_index;
        
        % left boundary cells
        if i==1
            % left edge is a boundary edge
            edge2el{edge_l} = [ind_el];
            transT(edge_l)   = 0;
            edgeBC(edge_l) = 1;
        else
            % left edge is not a boundary edge
            edge2el{edge_l} = [(i-2) * ncy + j, ind_el];
            transT(edge_l)   =   2*dy_j/(x(i+1) - x(i-1));
        end
        
        % bottom boundary cells
        if j==1
            % left edge is a boundary edge
            edge2el{edge_b} = [ind_el];
            transT(edge_b)   = 0;
            edgeBC(edge_b) = 1;
        else
            % left edge is not a boundary edge
            edge2el{edge_b} = [(i-1) * ncy + j-1, ind_el];
            transT(edge_b)   =  2* dx_i/(y(j+1) - y(j-1));
        end
        
        % add edge to edge2pts.
        edge2pt(edge_l, :) = [pt_lb, pt_lt]; % left
        edge2pt(edge_b, :) = [pt_lb, pt_rb]; % bottom
        
        % top boundary cells
        if j==ncy
            % top edge is a boundary edge
            edge2el{edge_t} = [ind_el];
            transT(edge_t)   = 0;
            edgeBC(edge_t) = 1;
            edge2pt(edge_t, :) = [pt_lt, pt_rt];
        else
            % top edge is not a boundary edge
            % do nothing
        end
        
        % right boundary cells
        if i==ncx
            % right edge is a boundary edge
            edge2el{edge_r} = [ind_el];
            transT(edge_r)   = 0;
            edgeBC(edge_r) = 1;
            edge2pt(edge_r, :) = [pt_rb, pt_rt];
        else
            % right edge is not a boundary edge
            % do nothing
        end
        
    end
    
end

toc;

edgeBC = logical(edgeBC);

% plot the mesh

if isplot
    figure(1111);
    Xs = reshape(Points(edge2pt',1), 2, N_edges);
    Ys = reshape(Points(edge2pt',2), 2, N_edges);
    plot(Xs, Ys, 'k');
    hold on;
    % plot boundary edges in red.
    plot(Xs(:, edgeBC), Ys(:, edgeBC), 'r','linew',2);
    plot(Cs(:, 1), Cs(:, 2), 'b.');
    xlabel('X');
    ylabel('Y');
    title('RectMesh for Finite Volume');
    daspect([1,1,1]);
    xlim(x([1, end]));
    ylim(y([1, end]));
    hold off;
end

end

