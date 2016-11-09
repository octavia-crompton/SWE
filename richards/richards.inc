************************************************************************

      parameter (nt = 11, nz=101)

************************************************************************

	    implicit real*8(a-h,o-z)

  	  common/reals/  stop_tol, alpha, theta_S, theta_R, 
     &               glambda, gn, gm, gKsat, dz, dt,
     &               htop, hbottom,
     &               gCbottom,gKbottom,thetabottom,
     &               gCtop,gKtop,thetatop
   	  common/ints/  kmin, kmax, itmin, itmax, 
     &              iBottomBoundCon
      common/param/  phi(6)
      common/array/  z(nz), t(nt), hinit(nz),
     &               thetainit(nz), gKinit(nz), gCinit(nz)
      common/array/  DeltaPlus(nz,nz), DeltaMinus(nz,nz),
     &               APlus(nz,nz), AMinus(nz,nz),
     &               Hsol(nz,nt), thetaSol(nz,nt)


