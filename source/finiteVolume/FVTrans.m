function T = FVTrans(cell2face, face2cell, faceBC, transT)
%
%   T = FVTrans(cell2face, face2cell, faceBC, transT)
%
%     Build the transmisibility matrix to compute flux from pressure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%       cell2face: a cell array Nc by 1 that contains index of faces of each cell
%
%       face2cell: a cell array Nf by 1 that contains index of cells for each interface
% 
%       faceBC  : a  Nf by 1 vector contains markers for boundary cells, 
%                      1 is boundary, 0 not boundary
%
%        transT   : a Nf by 1 vector contains  A_k/Delta_k 
%                 the ratio between surface area  A_k and the distance Delta_k
%                 between two cell centers
%
% output:
%    T: transmissibility matrix Nf by Nc to compute flux from pressure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nc = length(cell2face);
Nf  = length(face2cell);

T   = zeros(Nf, Nc);

for k = 1: Nf
    
    isboundary = faceBC(k);
    
    if isboundary
        continue
    end
    
    cells = face2cell{k};
    i = cells(1);
    j = cells(2);
    
    T(k,[i, j]) = [1, -1]*transT(k);
    
end

T = sparse(T);

end

