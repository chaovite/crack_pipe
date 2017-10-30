%DISLOC    U=disloc(m,x,nu)
%
%Returns the deformation at point 'x', given dislocation
%model 'm'. 'nu' specifies Poisson's ratio.
%
%Both 'm' and 'x' can be matrices, holding different models
%and observation coordinates in the columns.  In this case,
%the function returns the deformation at the points specified
%in the columns of 'x' from the sum of all the models in the
%columns of 'm'.  'x' must be 2xi (i = number of observation
%coordinates) and 'm' must be 10xj (j = number of models).
%
%The coordinate system is as follows: east = positive X,
%and north = positive Y. Depths should be given positive.
%
%The output, 'U', has the three displacement components:
%east, north, and up (in the rows). 
%
%The dislocation model is specified as: length, width,
%depth, dip, strike, east, north, strike-slip, dip-slip,
%and opening.  The coordinates (depth, east, and north)
%specify a point at the middle of the bottom edge of
%the fault for positive dips and the middle of the top
%edge for negative dips.
%
%The units of the displacements will be the units of the
%slip (opening); the units of the observation coordinates
%should be consistent with the units of the model.
%
%For more information, see: Okada, Y., Surface deformation
%due to shear and tensile faults in a half-space, Bull.
%Seismol. Soc. Am., 75, 1135–1154, 1985.
