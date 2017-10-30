%DISLOC3D  [U D S flag] = disloc3domp(m,x,mu,nu,nthreads)
%
% Returns the deformation at point 'x', given dislocation model 'm'.  'mu'
% specifies the shear modulus and 'nu' specifies Poisson's ratio.
%
% Both 'm' and 'x' can be matrices, holding different models and observation
% coordinates in the columns.  In this case, the function returns the
% deformation at the points specified in the columns of 'x' from the sum of
% all the models in the columns of 'm'.  'x' must be 3xi (i = number of
% observation coordinates) and 'm' must be 10xj (j = number of models).
%
% The coordinate system is as follows:
%   east = positive X, north = positive Y, up = positive Z.
% Observation coordinates with positive Z values return a warning.
%
% The outputs are 'U', the three displacement components: east, north, and up
% (on the columns); 'D', the nine spatial derivatives of the displacement: Uxx,
% Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, and Uzz (on the columns); and 'S', the 6
% independent stress components: Sxx, Sxy, Sxz, Syy, Syz, and Szz (on the
% columns). All these outputs have the same number of columns as 'x'.
%
% Output 'flag' is set for a singularity.
%
% The dislocation model is specified as:
%   length width depth dip strike east north strike-slip dip-slip opening.
% The coordinates (depth, east, north) specify a point at the middle of the
% bottom edge of the fault for positive dips and the middle of the top edge
% for negative dips. The coordinates (east, north) are specified just as for
% obs; however, note carefully that depth is a positive number in 'm', while
% in 'x' the equivalent coordinate is z and so is negative.
%
% nthreads is the number of OpenMP threads to use.
%
% Mind your units: for example, if lengths are given in km, then so should be
% slips.
