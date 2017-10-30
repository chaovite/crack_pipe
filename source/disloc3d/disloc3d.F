C     AMB. This line seems necessary in Matlab 2009b on a 64-bit Linux
C     machine.
#define __x86_64__
#include "fintrf.h"
C     disloc3d

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     AMB. "implicit none" prevents undeclared variables, possibly
C     resulting from a typo. As one example of why it should be used,
C     THETA was not declared in the original code and so must have
C     frequently been compiled as a real*4 instead of a real*8!
      implicit none
 
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

      mwPointer mxCreateDoubleMatrix, mxGetPr
      mwSize mxGetM, mxGetN

C     Oct 2010. AMB. Fixed a bug in earlier versions of disloc3d.F,
C     added some indentation and comments, and changed mex variable
C     types to conform with Matlab's API. (Hopefully the latter will
C     help with compilation issues.) The bug was that the displacements
C     and strains were being accumulated before being rotated by the
C     STRIKE rotation. Then a spurious rotation was performed at the
C     end. By circumstance, this spurious rotation used the value of
C     STRIKE in the final column of MODEL. Hence this function produced
C     erroneous results if multiple models were passed and at least one
C     model had a value for STRIKE that differed from the others.

      mwSize m,n,nstat,nmod
      mwPointer pMu,pNu,pModel,pStat
      mwPointer UOUT,DOUT,SOUT,FLAGOUT
      integer I,J,FLAG
      real*8 MODEL(10),STAT(3),MU,LAMBDA,NU,THETA
      real*8 ALPHA,STRIKE,DIP,DEG2RAD,SD,CD,SS,CS,DEPTH
      real*8 X,Y,Z,DISL1,DISL2,DISL3,AL1,AL2,AW1,AW2
      real*8 U(3),D(9),S(6),UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ

C     For Windows only! This resets the floating point exception to
C     allow divide by zero, overflow and invalid numbers.
C
C     integer(2) CONTROL
C     call GETCONTROLFPQQ(CONTROL)
C     CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
C     CONTROL = CONTROL .OR. FPCW$INVALID
C     CONTROL = CONTROL .OR. FPCW$OVERFLOW
C     call SETCONTROLFPQQ(CONTROL)

C     CHECK FOR PROPER NUMBER OF ARGUMENTS

      if (nrhs .ne. 4) then
         call mexErrMsgTxt('Usage: [U,D,S,flag]=disloc3d(m,x,mu,nu)')
      endif

C     CHECK THE ARGUMENTS

      M = mxGetM(prhs(1))
      nmod = mxGetN(prhs(1))
      if (M .ne. 10) then
         call mexErrMsgTxt('m must be 10x1 model vector')
      endif
      
      M = mxGetM(prhs(2))
      nstat = mxGetN(prhs(2))
      if (M .ne. 3) then
         call mexErrMsgTxt('x must be 3xn.')
      endif
      
      M = mxGetM(prhs(3))
      N = mxGetN(prhs(3))
      if ((M .ne. 1) .OR. (N .ne. 1)) then
         call mexErrMsgTxt('mu must be a scalar.')
      endif
      
      M = mxGetM(prhs(4))
      N = mxGetN(prhs(4))
      if ((M .ne. 1) .OR. (N .ne. 1)) then
         call mexErrMsgTxt('nu must be a scalar.')
      endif


C     CREATE A MATRIX FOR RETURN ARGUMENT

      plhs(1) = mxCreateDoubleMatrix(3, nstat, 0)
      plhs(2) = mxCreateDoubleMatrix(9, nstat, 0)
      plhs(3) = mxCreateDoubleMatrix(6, nstat, 0)
C     AMB. Change this from nstatx1 to 1xnstat to conform with the other
C     outputs.
      plhs(4) = mxCreateDoubleMatrix(1, nstat, 0)

C ASSIGN POINTERS TO THE VARIOUS PARAMETERS

      UOUT = mxGetPr(plhs(1))
      DOUT = mxGetPr(plhs(2))
      SOUT = mxGetPr(plhs(3))
      FLAGOUT = mxGetPr(plhs(4))

      pMODEL = mxGetPr(prhs(1))
      pSTAT = mxGetPr(prhs(2))
      pMU = mxGetPr(prhs(3))
      pNU = mxGetPr(prhs(4))

C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES

      call mxCopyPtrToReal8(pMU, MU, 1)
      call mxCopyPtrToReal8(pNU, NU, 1)
      LAMBDA = 2.d0*MU*NU/(1.d0 - 2.d0*NU)

C LOOP OVER STATIONS

      do 111 I = 1,nstat 

         do 120 J = 1,3
            U(J) = 0.d0
 120     continue

         do 130 J = 1,9
            D(J) = 0.d0
 130     continue   
         
         call mxCopyPtrToReal8(pSTAT+(I-1)*24, STAT, 3)
         
C     GENERATE WARNINGS FOR POSITIVE DEPTHS
         
         if (STAT(3) .GT. 0) then
            call mexEvalString("warning('Positive depth given.')")
         ELSE
            
C     LOOP OVER MODELS
            
            do 222 J = 1,nmod
               
               call mxCopyPtrToReal8(pMODEL+(J-1)*80, MODEL, 10)
               
C     AMB May 2011. Write this out to full double precision; otherwise,
C     one can't get results more accurate than 13 digits and likely
C     fewer.
C               DEG2RAD = 0.01745329251994d0
               DEG2RAD = 0.017453292519943295d0

               ALPHA = (LAMBDA + MU)/(LAMBDA + 2.d0*MU)

               STRIKE = (MODEL(5) - 90.d0)*DEG2RAD
               CS = dcos(STRIKE)
               SS = dsin(STRIKE)

               DIP = MODEL(4)
               CD = dcos(DIP*DEG2RAD)
               SD = dsin(DIP*DEG2RAD)

               DISL1 = MODEL(8)
               DISL2 = MODEL(9)
               DISL3 = MODEL(10)

               DEPTH = MODEL(3) - 0.5d0*MODEL(2)*SD

               AL1 = MODEL(1)/2.d0
               AL2 = AL1
               AW1 = MODEL(2)/2.d0
               AW2 = AW1

C     Convert from disloc3d to dc3d coordinate system. Let the rotation
C     matrix
C       R = [cs ss 0; -ss cs 0; 0 0 1].
C     Apply some translations plus R'. See below for the inverse.
               X = CS*(-MODEL(6) + STAT(1)) - SS*(-MODEL(7) + STAT(2))
               Y = -0.5*CD*MODEL(2) + SS*(-MODEL(6) + STAT(1)) +
     *              CS*(-MODEL(7) + STAT(2))
               Z = STAT(3)

C     AMB. Reinstate the check.
               if (I .eq. 1 .AND. ((MODEL(3)-SD*MODEL(2) .lt. 0) .or.
     *              (MODEL(1) .le. 0) .or.
     *              (MODEL(2) .le. 0) .or.
     *              (MODEL(3) .lt. 0))) then
                  call mexEvalString
     *                 ("warning('Unphysical model.')")
               endif

               call DC3D(ALPHA,X,Y,Z,DEPTH,DIP,
     *              AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,
     *              UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,FLAG)
               
C     SUM CONTRIBUTIONS FROM ALL PHYSICAL MODELS. AMB. This is now done
C     after rotating back from dc3d to disloc3d coordinate system. In
C     earlier versions, it was incorrectly done on the accumulated
C     results based on only the final model. The translation part of the
C     earlier transformation does not need to be inverted because the
C     results -- change in position and the derivative of this -- are
C     translation invariant.

C     Compute U = U + R*[UX UY UZ]'.
               U(1) = U(1) + CS*UX + SS*UY
               U(2) = U(2) - SS*UX + CS*UY
               U(3) = U(3) + UZ
               
C     Compute D = D + R*[UXX UXY UXZ; UYX UYY UYZ; UZX UZY UZZ]*R'.
               D(1) = D(1) + CS**2.d0*UXX + CS*SS*(UXY + UYX) +
     *              SS**2.d0*UYY
               D(2) = D(2) + CS**2.d0*UXY - SS**2.d0*UYX +
     *              CS*SS*(-UXX + UYY)
               D(3) = D(3) + CS*UXZ + SS*UYZ
               D(4) = D(4) - SS*(CS*UXX + SS*UXY) +
     *              CS*(CS*UYX + SS*UYY)
               D(5) = D(5) + SS**2.d0*UXX - CS*SS*(UXY + UYX) +
     *              CS**2.d0*UYY
               D(6) = D(6) - SS*UXZ + CS*UYZ
               D(7) = D(7) + CS*UZX + SS*UZY
               D(8) = D(8) - SS*UZX + CS*UZY
               D(9) = D(9) + UZZ
               
 222        continue
            
C     CALCULATE STRESSES:
C       sigma = 2*mu*[(D + D')/2] + lambda*sum(diag(D))*eye(3).

            THETA = D(1) + D(5) + D(9)
            
            S(1) = LAMBDA*THETA + 2.d0*MU*D(1)
            S(2) = MU*(D(2) + D(4))
            S(3) = MU*(D(3) + D(7))
            S(4) = LAMBDA*THETA + 2.d0*MU*D(5)
            S(5) = MU*(D(6) + D(8))
            S(6) = LAMBDA*THETA + 2.d0*MU*D(9)
            
C     COPY TO MATLAB
            
            call mxCopyReal8ToPtr(U, UOUT+(I-1)*24, 3)
            call mxCopyReal8ToPtr(D, DOUT+(I-1)*72, 9)
            call mxCopyReal8ToPtr(S, SOUT+(I-1)*48, 6)
            call mxCopyInteger4ToPtr(FLAG, FLAGOUT+(I-1)*8, 1)
            
         endif
 111  continue  

      return
      end
