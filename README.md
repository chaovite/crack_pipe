# pipe_crack
This repository contains the Matlab code that couples a multi-section pipes with multiple rectangular fluid filled cracks. The subroutines for SBP operators and imex time integrators are contributed by Ossian O'Reilly.

Currently, the 1d pipe model can be used to couple with multiple square cracks (assuming Poiseuille's flow across the width direction). However, the 2D radial symmetric pipe is only allow to be coupled to a crack at the bottom while the crack itself can be 3D tabular crack.

The coupled 2D pipe with 3D tabular crack model implements the governing equations in Liang et al., 2019, Part I, submitted JGR-Solid Earth. This model captures acoustic gravity waves in a 2D radially symmetric conduit and crack waves along 3D crack with proper treatment of viscous boundary layers in both the conduit and the crack.
