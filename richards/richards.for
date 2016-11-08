       program richards
C======================================================================
C     gC    =  specific moisture storage  (C)
C     gK     =  hydraulic conductivity  (K)
C     theta   =  vollumetric moisture content
************************************************************************
      include 'richards.inc'
      open(3,file ='richards.dat')

C C Read input, set up grid, initial conditions.
      call input

      
!       write(*,200) gCtop, gKtop, thetatop
!         write(*,*) gCbottom, gKbottom, thetabottom
!       write(*,201) thetaSOl(1:nz, 1)
!       write(*,201) dMminus(1:11, 1:11)
!       write(*,*) '  '
!       write(*, 201) dMminus(nz-10:nz,nz-10:nz)
!       write(*,*)

      do i = 1, nt
        hnp1m  = Hsol(1:nz,i)
        thetan = Thetasol(1:nz,i)
      enddo
      stop

 201  format(' ', 8f12.5)
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
      glambda = 0.25d0
      gKsat = 0.09d0
      
      gn = glambda + 1.d0
      gm = glambda/gn
      
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

      iBottomBoundCon = 1
      
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
      do i=1,nz
        call vanGenuchten(hinit(i), thetainit(i), gKinit(i), gCinit(i))
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
      if (h .lt. 0.d0) then
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