function [D] = FVDiv(cell2face, face2cell, faceBC)
% Construct finite vlume divergence operator
%
% Input:
%  cell2face: a cell array Nc by 1 that contains index of faces of each cell
%  face2cell: a cell array Nf by 1 that contains index of cells for each interface
%  faceBC  : a  Nf by 1 vector contains markers for boundary cells, 
%                      1 is boundary, 0 not boundary
%
% D: divergence matrix Nc by Nf that is used to sum all the flux from the
% non-boundary faces
%
% This implementation assumes zero-flow at the boundaries.
%

Nc = length(cell2face);
Nf  = length(face2cell);

D  = zeros(Nc, Nf);

for i = 1: Nc
    
    faces = cell2face{i};
    
    for j = 1: length(faces)
        k = faces(j); % index of face
        if faceBC(k)==1
            % boundary face skip
            continue
        end
        
        is_outward = face2cell{k}(1)==i;
        if is_outward
            D(i,k) = 1;
        else
            D(i,k) = -1;
        end
        
    end
    
end

D = sparse(D);

end

