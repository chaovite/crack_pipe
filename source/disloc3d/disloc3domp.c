/* Original: May 2011. AMB is rewriting disloc3d.F in C so we can use
   OpenMP. OpenMP can be used in F77, but for a variety of reasons, including
   issues with non-thread-safe matlab calls, we want the flexibility of C.

   To build:
     Serial:
       mex -O disloc3domp.c dc3omp.f
     Parallel:
       mex -DHMMVP_MEXSVD CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" disloc3domp.c dc3omp.f
*/

#include <mex.h>
#include <math.h>
#include <omp.h>

#ifdef __cplusplus
extern "C"
#endif
void dc3d_
(double* ALPHA, double* X, double* Y, double* Z, double* DEPTH, double* DIP,
 double* AL1, double* AL2, double* AW1, double* AW2, double* DISL1,
 double* DISL2, double* DISL3, double* UX, double* UY, double* UZ, double* UXX,
 double* UYX, double* UZX, double* UXY, double* UYY, double* UZY, double* UXZ,
 double* UYZ, double* UZZ);

#define DEG2RAD (M_PI / 180.0)
/* disloc3d.F and related files use the following definition of DEG2RAD. Only 13
   digits are specified. That is a mistake. At least 3 more should be appended
   to get full double precision. The results of this file are completely
   identical to those of disloc3d.F if this line, rather than the previous, is
   used.  */
//#define DEG2RAD 0.01745329251994
#define cosd(a) (cos((a)*DEG2RAD))
#define sind(a) (sin((a)*DEG2RAD))

void disloc3domp
(double *mdl, int nmdl, double *obs, int nobs, double mu, double nu,
 double *U, double *D, double *S, double *flag, int nthreads,
 void (*errfn)(const char *))
{
  double lambda, alpha;

  omp_set_num_threads(nthreads);

  lambda = 2.0*mu*nu/(1.0 - 2.0*nu);
  alpha = (lambda + mu)/(lambda + 2.0*mu);

  for (int i = 0; i < 3*nobs; i++) U[i] = 0;
  for (int i = 0; i < 9*nobs; i++) D[i] = 0;
  for (int i = 0; i < nobs; i++) flag[i] = 0;
  
  for (int k = 0; k < nmdl; k++) {
    double *m, strike, cs, ss, dip, cd, sd, cs2, ss2, csss,
      disl1, disl2, disl3, al1, al2, aw1, aw2, depth;

    m = mdl + 10*k;
    
    strike = m[4] - 90.0;
    cs = cosd(strike);
    ss = sind(strike);
    cs2 = cs*cs;
    ss2 = ss*ss;
    csss = cs*ss;

    dip = m[3];
    cd = cosd(dip);
    sd = sind(dip);

    disl1 = m[7];
    disl2 = m[8];
    disl3 = m[9];
    al1 = al2 = 0.5*m[0];
    aw1 = aw2 = 0.5*m[1];

    depth = m[2] - 0.5*m[1]*sd;

    if (errfn &&
	((m[2] - sd*m[1] < 0.0 && fabs((m[2] - sd*m[1]) / depth) > 1.0e-14)
	 || m[0] <= 0.0 || m[1] <= 0.0 || m[2] < 0.0))
      errfn("Unphysical model");

    /* NB: I'm using schedule(static) here because the behavior of the common
       blocks protected by THREADPRIVATE in dc3omp.f is not defined for
       schedule(dynamic). But this thread is so deterministic that 'static' is
       probably just fine anyway.  */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < nobs; i++) {
      double *o, *u, *d, x, y, z, ux, uy, uz,
	uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz;

      o = obs + 3*i;
      if (k == 0 && errfn && o[2] > 0.0) errfn("obs has a positive z value");

      // disloc3d coords -> dc3d coords

      // Let
      //     R = [cs ss 0; -ss cs 0; 0 0 1].
      // Apply some translations and R'.
      x = cs*(-m[5] + o[0]) - ss*(-m[6] + o[1]);
      y = -0.5*cd*m[1] + ss*(-m[5] + o[0]) + cs*(-m[6] + o[1]);
      z = o[2];

      // Okada rectangular dislocation
      dc3d_(&alpha, &x, &y, &z, &depth, &dip,
            &al1, &al2, &aw1, &aw2, &disl1, &disl2, &disl3,
            &ux, &uy, &uz,
            &uxx, &uyx, &uzx, &uxy, &uyy, &uzy, &uxz, &uyz, &uzz);

      // dc3d coords -> disloc3d coords

      // U = U + R*[UX UY UZ]'
      u = U + 3*i;
      u[0] +=  cs*ux + ss*uy;
      u[1] += -ss*ux + cs*uy;
      u[2] +=  uz;
      // D = D + R*[UXX UXY UXZ; UYX UYY UYZ; UZX UZY UZZ]*R'
      d = D + 9*i;
      d[0] +=  cs2*uxx + csss*(uxy + uyx) + ss2*uyy;
      d[1] +=  cs2*uxy - ss2*uyx + csss*(-uxx + uyy);
      d[2] +=  cs*uxz + ss*uyz;
      d[3] += -ss*(cs*uxx + ss*uxy) + cs*(cs*uyx + ss*uyy);
      d[4] +=  ss2*uxx - csss*(uxy + uyx) + cs2*uyy;
      d[5] += -ss*uxz + cs*uyz;
      d[6] +=  cs*uzx + ss*uzy;
      d[7] += -ss*uzx + cs*uzy;
      d[8] +=  uzz;
    }
  }

  // Stress
  for (int i = 0; i < nobs; i++) {
    double *s, *d, theta;

    d = D + 9*i;
    s = S + 6*i;
    theta = d[0] + d[4] + d[8];
    s[0] = lambda*theta + 2.0*mu*d[0];
    s[1] = mu*(d[1] + d[3]);
    s[2] = mu*(d[2] + d[6]);
    s[3] = lambda*theta + 2.0*mu*d[4];
    s[4] = mu*(d[5] + d[7]);
    s[5] = lambda*theta + 2.0*mu*d[8];
  }
}

void threadsafe_mexWarnMsgTxt(const char* msg)
{
#pragma omp critical (warn)
  mexWarnMsgTxt(msg);
}

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  // Input
  double *mdl, *obs, mu, nu;
  int nmdl, nobs, nthreads;
  // Output
  double *U, *D, *S, *flag;

  // Arguments
  if (nrhs < 4 || nlhs != 4)
    mexErrMsgTxt("Usage: [U D S flag] = disloc3domp(m,x,mu,nu,[nthreads])");

  if (mxGetM(prhs[0]) != 10) mexErrMsgTxt("m is a 10x(nmdl) array");
  mdl = mxGetPr(prhs[0]);
  nmdl = mxGetN(prhs[0]);

  if (mxGetM(prhs[1]) != 3) mexErrMsgTxt("x is a 3x(nobs) array");
  obs = mxGetPr(prhs[1]);
  nobs = mxGetN(prhs[1]);

  if (mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1)
    mexErrMsgTxt("mu is a scalar");
  mu = mxGetScalar(prhs[2]);

  if (mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 1)
    mexErrMsgTxt("nu is a scalar");
  nu = mxGetScalar(prhs[3]);

  nthreads = 4;
  if (nrhs > 4) {
    if (mxGetM(prhs[4]) != 1 || mxGetN(prhs[4]) != 1)
      mexErrMsgTxt("nthreads is a scalar");
    nthreads = (int)mxGetScalar(prhs[4]);
    if (nthreads < 1) mexErrMsgTxt("nthreads >= 1");
  }

  // Output
  plhs[0] = mxCreateDoubleMatrix(3, nobs, mxREAL);
  U = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(9, nobs, mxREAL);
  D = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(6, nobs, mxREAL);
  S = mxGetPr(plhs[2]);
  plhs[3] = mxCreateDoubleMatrix(1, nobs, mxREAL);
  flag = mxGetPr(plhs[3]);

  disloc3domp(mdl, nmdl, obs, nobs, mu, nu, U, D, S, flag, nthreads,
		&threadsafe_mexWarnMsgTxt);
}
