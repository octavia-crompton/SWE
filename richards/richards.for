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
      call vanGenuchten(0.1d0, theta, gK, gC)
      write(*,200) gC, gK, theta
      stop

 200  format(' ',3f15.8)
!  201  format(' ','time is ',f8.1,' maximum CFL number is ',f7.4 )
      end
      
************************************************************************
      subroutine input
      include 'richards.inc' 
!       character*72 dum, filename(100)
C     Read input from file 'richards.dat'.
      read(3,'(a72)') dum
      read(3,*) stop_tol
      read(3,'(a72)') dum
      read(3,*)  alpha, theta_S, theta_R, glambda
      read(3,'(a72)') dum
      read(3,*)  gKsat
      
      gn = glambda + 1.d0
      gm = glambda/gn
      
      phi(1) = alpha
      phi(2) = theta_S
      phi(3) = theta_R
      phi(4) = gn
      phi(5) = gm
      phi(6) = gKsat
      
      return
      end

************************************************************************
      subroutine vanGenuchten(h,theta, gK, gC)
      include 'richards.inc'
      real*8 theta
      
      theta = (theta_S - theta_R)/(1.d0+(alpha*abs(h))**gn)**gm +theta_R
      Se = ((theta - theta_R)/(theta_S - theta_R))
      gK = gKsat*Se**(1.d0/2.d0)*(1.d0-(1.d0-Se**(1.d0/gm))**gm)**2.d0
      if (h .lt. 0.d0) then
          gC =  alpha*gn*(1.d0/gn - 1.d0)*(alpha*abs(h))**(n - 1.d0)
     &      *(theta_R - theta_S)*((alpha*abs(h))**gn + 1)**(1/gn - 2)
      else
         write(*,*) 'should h be positive?'
          gC =  -alpha*gn*(1.d0/gn - 1.d0)*(alpha*abs(h))**(gn - 1.)
     &        *(theta_R - theta_S)*((alpha*abs(h))**gn + 1.d0)
     &        **(1.d0/gn - 2.d0)
      endif
      
      return
      end