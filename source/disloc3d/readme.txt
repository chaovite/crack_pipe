disloc3d v0.1
Crustal Deformation and Fault Mechanics Group, Stanford University

This package provides wrappers to the code DC3.F by Y. Okada. Cite:

    Okada, Y., 1992. Internal deformation due to shear and tensile faults in
    a half space, Bull. seism. Soc. Am., 82, 1018--1040.

The wrappers differ only in software details; the math and usage are the same.

In Matlab, type
    help disloc3d
to see basic usage. Read
    FaultDiagram.pdf
to learn the geometry conventions. In Matlab, type
    ex
to run an example, and
    edit ex.m
to view the code for example usage.

---------

This package contains three versions of the same underlying routine. All three
are used in Matlab. In increasing order of difficult of (a) getting these
working and (b) speed, they are as follows:

1. The first version is disloc3dpm.m. It is pure Matlab code and so is much slower
than the other two. But it can be used for small problems if you can't get
either of the next two versions to build.

2a. The second version is a serial mex file that is a basic wrapper to DC3.F. To
build it, type
    mex disloc3d.F dc3d.f
If this does not work for you, then you need to set up mex and possibly install
compilers. See
    http://www.mathworks.com/help/matlab/ref/mex.html

2b. We have implemented a fix to a numerical problem that occurs in and nearly
in the plane of a rectangular element very close to the rays extending from the
edges. This numerical error is essentially never an issue for most users of this
code, so you should not worry about it. We provide the fix in a separate version
of the code in case you wish to try it. Type
    mex disloc3d.F dc3quadrant.f

3. The third version is also a mex file. It uses OpenMP to parallelize the
operations. It also compensates for the numerical error reported in 2b.
  For a variety of reasons, the code can be much harder to build than version 2.
On a UNIX system, attempt to build it as follows:
     mex -DHMMVP_MEXSVD CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" disloc3domp.c dc3omp.f

---------

Version History:
0.0. 100812. AMB. First release for general use.
0.1. 111612. AMB. New method based on xi,eta>0 quadrant to get rid of
  cancellation errors. Implemented in dc3quadrant.f and dc3omp.f but *not*
  dc3d.f.
