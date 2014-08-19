%------------------------------------------------------------------------------
README FAIR.2011
%------------------------------------------------------------------------------

W E L C O M E   T O   F A I R:

  FLEXIBLE ALGORITHMS FOR IMAGE REGISTRATION

The package is described in the SIAM book, 
  Jan Modersitzki: FAIR - flexible algorithms for image registration, SIAM 2009
  see http://www.ec-securehost.com/SIAM/FA06.html

This package comes with the following files and folders:
FAIRcopyright.m   your copyright agreement
FAIRstartup.m     a startup file, adding all paths
README.m          this file 
data              a collection of test cases
distances         a collection of distance measures
examples          a collection of examples
interpolation     the interpolation toolbox
kernel            the optimization kernel 
landmarks         landmark based registration related files
matrixfree        matrix-free related code 
regularization    the regularization toolbox
tests             a (draft) test suite 
transformations   the transformation toolbox
viewer            the image visualization toolbox% 

A good starting point is to enter the directory and type FAIRstartup at the MATLAB 
prompt to get all the paths set right. Explore then the various examples and tutorials.

%------------------------------------------------------------------------------

NEWS in 2011
We changed a few things and added a couple of new examples and features, 
which are summarized as follows.

GENERAL
  Minimal examples for various mfiles, just run the mfiles with no arguments.

DATA and SETUP files (new)
  setup2Ddisc2CData                      disc to C, see Christensens 1995
  setup2DEPIData                         2D slices of Echo Planar MRI (heavily distorted)
  setup2DGaussianData                    academic example for mass preservation, see F. Gigengack 
  setup2DGaussianNoisyData               as above plus noise
  setup2DMPsquare                        academic example for mass preservation, see F. Gigengack 
  setup3DmiceData                        3D cardiac PET images from mice, systoly and diastoly, see F. Gigengack
  setup3DphantomData                     3D  Echo Planar MRI's of a water bottle (heavily distorted)
  checkSetupDataFile                     debugging tool

  plus some images in term of *.jpg and *.mat.

INTERPOLATION
  We merge the schemes for dimensions 1,2, and 3 into one m-file (no flag specifying 
  the dimension is required), and added CPP versions and wrapper, see setup3DbrainData 
  for an example involving splineInterMex.

  nnInter.m                              next neighbour interpolation 
  linearInter.m                          linear interpolation
  linearInterMex.m                       wrapper for linear interpolation using CPP code 
  linearInterMexC.cpp                    CPP file for linear interpolation 
  linearInterSmooth.m                    a smooth C1 variant of linear interpolation
  linearInterSmoothMex.m                 analog to linear
  linearInterSmoothMexC.cpp              -"-
  linearInterMatlab.m                    linear interpolation based on MATLAB's scheme and 
                                         finite difference approximation to derivatives
  splineInter.m                          spline interpolation
  splineInterMex.m                       analog to linear
  splineInterMexC.cpp                    -"-
  cubicInter.m                           See:
                                         Catmull, E., and Rom,  R. A class of local interpolating splines.
                                         In Computer Aided Geo Design, R. E. Barnhill and R. F. Reisenfeld,
                                         Eds. Academic Press, New York, 1974, pp. 317-326.
  cubicInterMex.m                        analog to cubicInter
  cubicInterMexC.cpp                     -"-

  motherSpline: more efficient matlab implementation
  fixed BUG in splineIner: Valid = @(j) (-1<x(:,j) & x(:,j)<m(j)+2);

TRANSFORMATION
  We implemented memory efficient versions of affine linear and spline transformations. 

  affine3Dsparse                        affine linear transformation for 3D (memory efficient)
  splineTransformation2Dsparse          spline transformation for 2D (memory efficient)
  splineTransformation3Dsparse          spline transformation for 3D (memory efficient)
  tensorProd.c                          C implementation of kron(Q3,Q2,Q1) * w, 
                                        see splineTransformation3Dsparse.m
DISTANCES
  We added various C files and a weighted SSD measure.

KERNEL
  Substituted getCenteredGrid by getCellCenteredGrid
  GaussNewtonArmijo is no longer supported, please use GaussNewton instead.

TESTSUITES
  We also included are suite of tests for debugging the software.

EXAMPLES
  We included many more examples and aimed to make the examples more consistent.
