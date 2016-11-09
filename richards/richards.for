       program richards
C======================================================================
C     gC  (C)  =  specific moisture storage  
C     gK  (K)  =  hydraulic conductivity  
C     theta   =  vollumetric moisture content
C     Sulution arrays:
C       Hsol (H) ,   ThetaSol (THETA)
C     Hn  (hnp1m)
C     Bottom boundary condition:
C       iBottomBoundCon (BottomBoundCon)
C     Within the iteration solver, vanGenuchten ouput is:
C       dumtheta (thetanp1m), dumK (knp1m), dumC (cnp1m)
C     dMMinus (MMminus)
C     barKplus (kbarplus)
C     barplusK (Kbarplus)
C     force boundary conditions:
C       Hn1  (hnp1mp1)
C    dumThetaP1 (thetanp1mp1 )
************************************************************************
      include 'richards.inc'
      dimension R_MPFD(nz), dumPlus(nz, nz), dumMinus(nz, nz),
     &          dumPlusH(nz), dumMinusH(nz), test(nz, nz),
     &          Ainv(nz, nz), deltam(nz), Hn1(nz), dumThetaP1(nz)

C Read input, set up grid, initial conditions.
      call input   
!         write(*, 201) C(nz-10:nz,nz-10:nz)
!         print *, dMplus*kdum
      do i = 1,nt
        Hn  = Hsol(1:nz,i)   
        thetan = Thetasol(1:nz,i)
        istopFlag = 0
        niter = 0
        do while (istopFlag.eq.0)
          do k=1,nz
            call vanGenuchten(Hn(k), dumtheta(k), dumK(k), dumC(k))
          enddo
          call fdiag(hinit, C)
          barKplus = matmul(dMPlus, dumK)
          call fdiag(barKplus, barplusK)
          barKminus = matmul(dMMinus, dumK)
          call fdiag(barKminus, barminusK) 
          dumPlus = matmul(barplusK, DeltaPlus)
          dumPlusH = matmul(dumPlus, Hn)
          dumMinus = matmul(barminusK, DeltaMinus) 
          dumMinusH = matmul(dumMinus, Hn)
          A = (1.d0/dt)*C-1.d0/(dz**2.d0)*(dumPlus - dumMinus)
C         Compute the residual of MPFD (RHS)
          R_MPFD = (1.d0/(dz**2.d0))*(dumPlusH - dumMinusH) +
     &         (1.d0/dz)*(barKplus - barKminus) - 
     &          1.d0/dt*(dumtheta - thetan)
C         Compute deltam for iteration level m+1
          call inverse(A,Ainv,nz)
          deltam = matmul(Ainv, R_MPFD)
          niter = niter + 1
          if (maxval(abs(deltam(1:(nz-1)))).lt. stop_tol) then
              istopFlag = 1
              Hn1 = Hn + deltam
              Hn1(1) = htop
              if (BottomBoundCon.eq.0) then
                Hn1(nz) = hbottom
              elseif (BottomBoundCon.eq.1) then
                Hn1(nz) = Hn1(nz-1)
              endif
              do k=1,nz
                call vanGenuchten(Hn1(k), dumtheta(k), dumK(k), dumC(k))
                dumthetaP1(k) =dumtheta(k)
              enddo
          else
              Hn1 = Hn + deltam
              Hn = Hn1 ! Force boundary conditions
              Hn(1) = htop;
              if (BottomBoundCon.eq.0) then
                Hn1(nz) = hbottom
              elseif (BottomBoundCon.eq.1) then
                Hn1(nz) = Hn1(nz-1)
              endif
          endif
        enddo
        thetaSol(1:nz,i) = dumthetaP1 
        Hsol(1:nz,i) = Hn1 
      enddo
      
      stop

 201  format(' ', 11f10.1)
 200  format(' ', 3f15.10)
!  201  format(' ','time is ',f8.1,' maximum CFL number is ',f7.4 )
      end
      
************************************************************************
      subroutine input
      include 'richards.inc' 
      
      stop_tol = 0.01d0
      alpha = 0.0335d0
      theta_S = 0.368d0
      theta_R = 0.102d0
      glambda = 0.25d0  !  lambda --> glambda
      gKsat = 0.09d0    !  Ksat --> gKsat
      
      gn = glambda + 1.d0   !  n --> gn
      gm = glambda/gn       ! m -->  gm
      
      phi(1) = alpha
      phi(2) = theta_S
      phi(3) = theta_R
      phi(4) = gn
      phi(5) = gm
      phi(6) = gKsat
      
      dz = 1.d0
      kmin = 1
      kmax = nz
      do j=kmin,kmax
        z(j+1) = z(j) + dz
      enddo
      
      dt=1800
C      Define time variables
      t(1) = 0.d0
      do i=1,nt
        t(i+1) = t(i) + dt
      enddo
      htop = -30.d0
      hbottom = -1000.d0    
      hinit(1) = htop 
      hinit(2:nz) = hbottom 

      iBottomBoundCon = 1  ! BottomBoundCon  --> iBottomBoundCon
      
      DeltaPlus = 0.d0
      DeltaMinus = 0.d0
      do i = 2,nz-1
        DeltaPlus(i,i) = -1.d0
        DeltaMinus(i,i) = 1.d0
        DeltaMinus(i-1,i) = -1.d0
      enddo
      
      do i = 2,nz-1
        dMPlus(i,i) = 1.d0
        dMPlus(i,i-1) = 1.d0  
        dMMinus(i,i) = 1.d0
        dMMinus(i-1,i) = 1.d0        
      enddo
      dMPlus(1,1) = 2.d0
      dMPlus(nz,nz) = 2.d0    
      dMPlus(nz,nz-1) = 1.d0            
      dMPlus(2,1) = 0.d0

      dMMinus(1,1) = 2.d0 
      dMMinus(nz,nz) = 2.d0           
      
      call vanGenuchten(htop, thetatop, gKtop, gCtop)
      call vanGenuchten(hbottom, thetabottom, gKbottom, gCbottom)
      do k=1,nz
        call vanGenuchten(hinit(k), thetainit(k), gKinit(k), gCinit(k))
      enddo
      
      Hsol(1:nz,1) = hinit
      thetaSol(1:nz,1) = thetainit
            
      return
      end

************************************************************************
      subroutine vanGenuchten(h,theta, gK, gC)
      include 'richards.inc'
      theta = (theta_S - theta_R)/(1.d0+(alpha*abs(h))**gn)**gm +theta_R
      Se = ((theta - theta_R)/(theta_S - theta_R))
      gK = gKsat*Se**(1.d0/2.d0)*(1.d0-(1.d0-Se**(1.d0/gm))**gm)**2.d0
      if (h .le. 0.d0) then
        gC =  alpha*gn*(1.d0/gn - 1.d0)*(alpha*abs(h))**(gn - 1.d0)
     &        *(theta_R - theta_S)*((alpha*abs(h))**gn + 1)**(1/gn - 2)
      else
        write(*,*) 'should h be positive?'
        gC =  -alpha*gn*(1.d0/gn - 1.d0)*(alpha*abs(h))**(gn - 1.)
     &        *(theta_R - theta_S)*((alpha*abs(h))**gn + 1.d0)
     &        **(1.d0/gn - 2.d0)
      endif
      
      return
      end   

************************************************************************
      subroutine fdiag(h, hh)
      include 'richards.inc'
      real(kind=8), dimension(nz) :: h
      real(kind=8), dimension(nz, nz) :: hh
      do k=1,nz
        hh(k,k) = h(k)
      enddo
      
      return
      end
************************************************************************
      
      subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse