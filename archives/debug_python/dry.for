	program swe2d
C======================================================================
C This program uses a finite volume discretization to solve the
C 2-D shallow water equations on an arbitrary mesh consisting of
C quadrilaterals.  Roe's Approximate Riemann Solver is used along
C with a MUSCL reconstruction of primitive variables to obtain
C 2nd order accuracy in space.  Hancock's	Predictor-Corrector method  
C is used to acheive 2nd order accuracy in time.  The model can handle
C a variably sloping bed as well as wetting and drying in the domain.
C Bed drag is parameterized by the Manning equation.  
C
C List of variables:
C
C     dt     = time increment.
C     epsh   = depth tolerence for dry bed problems.
C     grav   = acceleration due to gravity.
C     inum   = number of boundary faces in a boundary cell.
C     ipos   = orientation of boundary face:=1 bottom,=2 right,=3 top,
C              =4 left.
C     nbcell = number of boundary cells.
C     ndir   = number of bodes at which velocity, depth, volume flux
C              etc. are specified.
C     xsplit  = dividing line for dam-break problem.
C     tmax   = time of computation.
C     xn     = Manning bed roughness coefficient.
************************************************************************
      include 'dry.inc'
      real flux1, flux2, flux3, flux4
      real  vol, vol0, dvol, flux, infl
      open(2,file ='coords')
      open(3,file ='dryin.dat')
C     open(7,file ='object.dat')
C     open(8,file ='dryout.dat')
C     open(9,file ='movie.dat')
      open(11,file='restart.out')
c	    open(13,file='slice.out')
	    open(100,file='myout.out')
	    open(101,file='time.out')
      open(104,file='fluxes.out')
      open(120,file='dvol.out')                  
      open(1, file = 'test.out')
      call cpu_time(start)       
C Read input, set up grid, initial conditions.
      call input
      if(ifront .eq. 1 .or. imass .eq. 1) open(12,file='diag.out')
	    write(101, *) "time  |  print step  |  time step"            
      iprt = 0 
C Begin time loop.
      do while (t.lt.tmax)
	      it=it+1
        if(t.gt.60.)   dt=0.4
        if(t.gt.6960.) dt=0.4
C       nprt = tmax/100/dt        
c save every step for a movie     
C	      call movie
        t = t + dt
        amax = 0.D0
        write(*, *) t
        
C Compute predictor.
        do j=1,ncol
         do k=kbeg(j),kend(j)
          call bconds(j,k,h,u,v)
          call predict(j,k)
         enddo
        enddo
C Loop over cells to compute fluxes.
        flux1 = 0.D0
        flux2 = 0.D0
        flux3 = 0.D0
        flux4 = 0.D0
        do j=1,ncol
        do k=kbeg(j),kend(j)
          call bconds(j,k,hp,up,vp)
          call fluxes(j-1,j,k,k,1)          ! vertical faces.   
          call fluxes(j,j,k-1,k,2)          ! horizontal faces.
          do i=1,inum(j,k)
            if(ipos(j,k,i) .eq. 3) then
              call fluxes(j,j,k,k+1,2)       ! top boundaries. 
              flux3 =  flux3 - f(j,k+1,1, 2)*ds(j,k+1, 2)*dt
            elseif(ipos(j,k,i) .eq. 2) then
              call fluxes(j,j+1,k,k,1)       ! right boundaries.
              flux2 = flux2  - f(j+1,k,1, 1)*ds(j+1, k, 1)*dt
            elseif(ipos(j,k,i) .eq. 1) then                 ! lower boundary
              flux1 =  flux1 + f(j,k,1, 2)*ds(j, k, 2)*dt
            elseif(ipos(j,k,i) .eq. 4) then                 ! left boundary
              flux4 =  flux4 + f(j,k,1, 1)*ds(j, k, 1)*dt                
            endif
          enddo
        enddo
        enddo
       
C         do j=0,ncol+1
C          do k=0,kend(j)+1
C           write(1, *) j,k,f(j,k,1,2)
C          enddo
C         enddo
C Compute corrector solution.
        infl = 0.d0
        do j=1,ncol
        do k=kbeg(j),kend(j)
          call source(j,k,hp(j,k),up(j,k),vp(j,k))
          infl = infl + dt*qs(1)*area(j,k)
          do l=1,3
            q(j,k,l) = q(j,k,l) + dt/area(j,k)*(
     &        f(j,k,l,2)*ds(j,k,2) + f(j,k,l,1)*ds(j,k,1) -  
     &        f(j+1,k,l,1)*ds(j+1,k,1) - f(j,k+1,l,2)*ds(j,k+1,2)) 
     &        + dt*qs(l)
          enddo	
        enddo
        
        enddo

C Store solution. 
C Solve continuity equation even in all cells, but 
C solve momentum equations in wet cells only.
       do j=1,ncol
       do k=kbeg(j),kend(j)
C Check for negative depth.
        if(q(j,k,1) .gt. 0.D0) then
	        h(j,k) = q(j,k,1)
	      else
	        q(j,k,1) = 0.D0
          h(j,k) = 0.D0
        endif
C Neglect momentum in nearly dry cells.
        if(h(j,k) .lt. epsh) then  
           u(j,k) = 0.d0
           v(j,k) = 0.d0
	         t0(j,k) = t 
           do l=2,3
            q(j,k,l) = 0.d0
           enddo
        elseif(h(j,k) .ge. epsh) then  
          u(j,k) = q(j,k,2)/h(j,k)
          v(j,k) = q(j,k,3)/h(j,k)
        endif
       enddo
       enddo 
       write(*,201) t, amax*dt
       if (t .gt. tmax ) then
         write(*, *) '!!!'
         do j=0,ncol+1
          do k=0,kend(j)+1
           write(1, *) j,k, h(j,k)
          enddo
         enddo       
       endif
C     Record infiltration and volume change
        vol0 = vol
        vol = 0.D0
C        call myinterp
        do j=1,ncol 
          do k=kbeg(j),kend(j)
            call source(j,k,hp(j,k),up(j,k),vp(j,k))
            vol = vol + h(j,k)*area(j,k)
          enddo
        enddo 
        dvol = vol - vol0
        flux =  flux1 + flux2 + flux3 + flux4     
        write(120,200) t, dvol, flux, infl         
C   Below here only executed every nprt time steps        
      iprt = iprt + 1
      if(iprt .eq. nprt) then
	      iprt=0
        itp = itp +  1
        call myoutput
        
      endif 
      enddo

C       if (t .eq. 0.1D0) then
C         write(*, *) '!!!'
C         do j=0,ncol+1
C          do k=0,kend(j)+1
C           write(1, *) j,k, h(j,k)
C          enddo
C         enddo
C       endif
      
      call cpu_time(finish)
      print '("Time = ",f6.3," seconds.")',finish-start
  
C End of time loop.
C     call movie
c	    call object
c	    call output
      stop
 200  format(' ',f8.1, 5e19.8)
 201  format(' ','time is ',f8.1,' maximum CFL number is ',f7.4 )
 202  format(' ', 2i4, 5e23.15) 
      end
************************************************************************
      subroutine myoutput
      include 'dry.inc'

C 	  file 101 is 'time.out'  -  to keep track of the time stpes
      write(101, 203)   t, itp, it
C 	  file 100 is 'myout.out' 
	    write(100, 202)  itp
      do j=1,ncol
        do k=kbeg(j),kend(j)
          write(100, *) j, k, h(j,k)
        enddo
      enddo

C    Loop over cells to compute fluxes.  file 104 is 'fluxes.out'
      do j=1,ncol
        do k=kbeg(j),kend(j)
          do i=1,inum(j,k)
            if( ipos(j,k,i) .eq. 3) then   ! horizontal boundaries
C                   write(104, *) j,k, f(j, k,1, 2)
                write(104, *) j,k,f(j,k+1,1, 2)*ds(j,k+1, 2)
            endif
          enddo
        enddo
      enddo
      write(104, 202)  itp
      
      return
 200  format(' ', i6, i6, 5e15.7, 5e15.7)
 202  format(' ', 4i6)	 
 203  format(' ', f7.1 , 4i6, 4i6 ) 
      end			
      
************************************************************************
      subroutine source(j,k,hdum,udum,vdum)
      include 'dry.inc'

      if(hdum .ge. epsh) then
	      tnew = dmax1(0.d0, t - t0(j,k))
        if(tnew .le. 0.d0) then
           znew = 0.d0
          elseif(tnew .gt. 0.d0 .and. tnew .le. tclip) then
           znew = zslope*tnew
          elseif(tnew .gt. tclip .and. tnew .le. tc) then
           znew = xk*tnew**ainflt 
          elseif(tnew .gt. tc) then
           znew = xk*tc**ainflt + binflt*(tnew - tc)
        endif
        told = dmax1(0.d0, t - t0(j,k) - dt)
	      if(told .le. 0.d0) then
          zold = 0.d0
        elseif(told .gt. 0.d0 .and. told .le. tclip) then
          zold = zslope*told
        elseif(told .gt. tclip .and. told .le. tc) then
          zold = xk*told**ainflt
        elseif(told .gt. tc) then
           zold = xk*tc**ainflt + binflt*(told - tc)
	      endif
        winflt = (zold - znew)/dt
        vmag = dsqrt(udum*udum + vdum*vdum)
        fricx = grav*xn*xn*udum*vmag/hdum**(1.D0/3.D0)
        fricy = grav*xn*xn*vdum*vmag/hdum**(1.D0/3.D0)
        qs(1) = winflt
	      qs(2) = 0.5D0*udum*winflt - fricx - grav*hdum*sx(j,k)
        qs(3) = 0.5D0*vdum*winflt - fricy - grav*hdum*sy(j,k)
      else
        qs(1) = 0.d0
        qs(2) = 0.d0
        qs(3) = 0.d0
      endif

      return
      end
************************************************************************
      subroutine fluxes(jl,jr,kl,kr,i1)
      include 'dry.inc'
C MUSCL extrapolation at cell interface.
      hl = hp(jl,kl) + 0.5D0*dh(jl,kl,i1)
      ul = up(jl,kl) + 0.5D0*du(jl,kl,i1)
      vl = vp(jl,kl) + 0.5D0*dv(jl,kl,i1)
      hr = hp(jr,kr) - 0.5D0*dh(jr,kr,i1)
      ur = up(jr,kr) - 0.5D0*du(jr,kr,i1)
      vr = vp(jr,kr) - 0.5D0*dv(jr,kr,i1)
      snn = sn(jr,kr,i1)
      cnn = cn(jr,kr,i1)
	    if(i1 .eq. 1) then
        dx =  deta(jr,kr,2)*area(jr,kr) 
	      dy = -deta(jr,kr,1)*area(jr,kr)
      else
        dx = -dxi(jr,kr,2)*area(jr,kr) 
	      dy =  dxi(jr,kr,1)*area(jr,kr)
      endif
C Needed for dry bed problems.
	    if(hl .lt. 0.D0) hl = 0.D0
      if(hr .lt. 0.D0) hr = 0.D0
C Compute arithmatic averages for source terms.
      havg = 0.5D0*(hl + hr)
      uavg = 0.5D0*(ul + ur)
      vavg = 0.5D0*(vl + vr)
C Prevent leakage into cells with higher bed elevation.
      etal = hp(jl,kl) + zc(jl,kl)
      etar = hp(jr,kr) + zc(jr,kr)
C Fluxes and source terms.
      if(havg .le. 0.D0 ) then
	      do i=1,3  
          f(jr,kr,i,i1) = 0.D0
C           qsource(jr,kr,i,i1) = 0.d0
        enddo
      else 
        call solver(hl,hr,ul,ur,vl,vr,fdum,snn,cnn,dx,dy,sdumx,sdumy,sw,
     &              winflt)
	      do i=1,3
          f(jr,kr,i,i1) = fdum(i)
        enddo
      endif

	    return
      end
      
************************************************************************
      subroutine solver(hl,hr,ul,ur,vl,vr,f,sn,cn,dx,dy,sx,sy,sw,winflt)
      implicit real*8(a-h,o-z)
      dimension ws(3), e(3,3), a(3), astar(3), da(3), f(3), dum(3), z(3)
      common/m/ grav, amax, xn, epsh
C Compute Roe averages at cell face.
      duml  = dsqrt(hl)
      dumr  = dsqrt(hr)
      hhat  = duml*dumr
      uhat  = (duml*ul + dumr*ur)/(duml + dumr)
      vhat  = (duml*vl + dumr*vr)/(duml + dumr)
      chat  = dsqrt(0.5D0*grav*(hl + hr))
      uperp = uhat*cn + vhat*sn
C Compute eigenvalues.
      a(1) = uperp - chat
      a(2) = uperp
      a(3) = uperp + chat
C Compute approximate wave strengths.
      dh    = hr - hl
      du    = ur - ul
      dv    = vr - vl
      dupar = -du*sn + dv*cn
      duperp=  du*cn + dv*sn
      ws(1) = 0.5D0*(dh - hhat*duperp/chat)
      ws(2) = hhat*dupar     
      ws(3) = 0.5D0*(dh + hhat*duperp/chat)
C Compute right eigenvectors.
      e(1,1) = 1.D0
      e(2,1) = uhat - chat*cn
      e(3,1) = vhat - chat*sn
      e(1,2) = 0.D0
      e(2,2) = -sn
      e(3,2) =  cn
      e(1,3) = 1.D0
      e(2,3) = uhat + chat*cn
      e(3,3) = vhat + chat*sn
C Entropy fix.
	    dl = dsqrt(dx*dx + dy*dy)
      cl = dsqrt(grav*hl)
      cr = dsqrt(grav*hr)
      uperpl = ul*cn + vl*sn
      uperpr = ur*cn + vr*sn
      al1 = uperpl - cl
      al3 = uperpl + cl
      ar1 = uperpr - cr
      ar3 = uperpr + cr
      da(1) = dmax1(0.D0, 4.D0*(ar1 - al1))
      da(2) = 0.d0
      da(3) = dmax1(0.D0, 4.D0*(ar3 - al3))
      do i=1,3   
         if(dabs(a(i)) .lt. 0.5D0*da(i)) then
           astar(i) = a(i)*a(i)/da(i) + 0.25D0*da(i)
          else
           astar(i) = dabs(a(i))
         endif
      if(astar(i)/dl .gt. amax) amax = astar(i)/dl
      enddo
C Compute flux increments.
      do i=1,3
        dum(i) = 0.D0
        do l=1,3
          dum(i) = dum(i) + (astar(l)*ws(l))*e(i,l)   ! - sw*z(l)
        enddo
      enddo
C Add flux to appropriate cell faces.
      f(1) = 0.5D0*(f1(hl,uperpl) + f1(hr,uperpr) - dum(1))
      f(2) = 0.5D0*(f2(hl,ul,uperpl,cn) 
     &            + f2(hr,ur,uperpr,cn) - dum(2))
      f(3) = 0.5D0*(f3(hl,vl,uperpl,sn)   
     &            + f3(hr,vr,uperpr,sn) - dum(3))

      return
      end
************************************************************************
      subroutine predict(j,k)
      include 'dry.inc'

      do kk=1,2                  ! loop over coord. directons.
        if(kk .eq. 1)then
          jr = j + 1
          jl = j - 1
          kr = k
          kl = k
        else
          jr = j
          jl = j
          kr = k + 1
          kl = k - 1
        endif
C Compute gradients only in wet cells.
	      if(h(j,k) .ge. epsh) then
C Limit free surface elevation to reduce dissipation..
          dh1 = h(j,k) + zc(j,k) - h(jl,kl) - zc(jl,kl)
          dh2 = h(jr,kr) + zc(jr,kr) - h(j,k) - zc(j,k)
C Needed to minimize dissipation at wet/dry interfaces.
          if(h(jl,kl) .lt. epsh) dh1 = 2.d0*dh1
          if(h(jr,kr) .lt. epsh) dh2 = 2.d0*dh2
          call limitr(ilim,beta,dh1,dh2,dhh)
          dh(j,k,kk) = dhh - dz(j,k,kk)	    
C U velocity.
          du1 = u(j,k) - u(jl,kl)
          du2 = u(jr,kr) - u(j,k)
          call limitr(ilim,beta,du1,du2,duu)
          du(j,k,kk) = duu
C V velocity.
          dv1 = v(j,k) - v(jl,kl)
          dv2 = v(jr,kr) - v(j,k)
          call limitr(ilim,beta,dv1,dv2,dvv)
          dv(j,k,kk) = dvv
	      else
          dh(j,k,kk) = 0.D0
          du(j,k,kk) = 0.D0
          dv(j,k,kk) = 0.D0
        endif
      enddo
C Generalized velocities.
	    uxi = u(j,k)*dxi(j,k,1) + v(j,k)*dxi(j,k,2)
	    ueta = u(j,k)*deta(j,k,1) + v(j,k)*deta(j,k,2)
C Predictor.
      call source(j, k, h(j,k), u(j,k), v(j,k))
      if(h(j,k) .ge. epsh**(0.75d0)) then
        qs(2) = qs(2)/h(j,k)
        qs(3) = qs(3)/h(j,k)
      else
        qs(1) = 0.d0
        qs(2) = 0.d0
        qs(3) = 0.d0
	    endif
      hp(j,k) = h(j,k) - 0.5D0*dt*(
     &    uxi*dh(j,k,1) + h(j,k)*(dxi(j,k,1)*du(j,k,1) + 
     &                            dxi(j,k,2)*dv(j,k,1)) +   
     &    ueta*dh(j,k,2) + h(j,k)*(deta(j,k,1)*du(j,k,2) +
     &                             deta(j,k,2)*dv(j,k,2)) + qs(1))   
      up(j,k) = u(j,k)	- 0.5D0*dt*(
     &    grav*dxi(j,k,1)*dh(j,k,1) + uxi*du(j,k,1) + 
     &    grav*deta(j,k,1)*dh(j,k,2) + ueta*du(j,k,2) + qs(2)) 
      vp(j,k) = v(j,k) - 0.5D0*dt*(
     &    grav*dxi(j,k,2)*dh(j,k,1) + uxi*dv(j,k,1) +
     &    grav*deta(j,k,2)*dh(j,k,2) + ueta*dv(j,k,2) + qs(3))
C Correct any negative depths.
	    if(hp(j,k) .lt. 0.d0) then
        hp(j,k) = 0.d0
        dh(j,k,1) = 0.d0
        dh(j,k,2) = 0.d0
      endif
C Neglect momentum in nearly dry cells.
	    if(hp(j,k) .le. 0.d0) then
        up(j,k) = 0.D0
        vp(j,k) = 0.D0
        do i=1,2
          du(j,k,i) = 0.d0
          dv(j,k,i) = 0.d0
        enddo
      endif

      return
      end
************************************************************************
      subroutine bconds(j,k,hdum,udum,vdum)
      include 'dry.inc'
      dimension hdum(0:nx,0:ny), udum(0:nx,0:ny), vdum(0:nx,0:ny)
C Loop over all boundary faces in the cell. 
      do i=1, inum(j,k)
        if(ipos(j,k,i) .eq. 1) then  	 ! front face.
            jj = j
            kk = k-1
            jl = j
            kl = k 
            j2 = j
            k2 = k+1
            io = 2 
        elseif(ipos(j,k,i) .eq. 2) then	! right face.
            jj = j+1
            kk = k 
            jl = j+1 
            kl = k  
            j2 = j-1
            k2 = k
            io = 1
        elseif(ipos(j,k,i) .eq. 3) then   ! back face.
            jj = j
            kk = k+1
            jl = j
            kl = k+1
            j2 = j
            k2 = k-1 
            io = 2 
        elseif(ipos(j,k,i) .eq. 4) then  ! left face.
            jj = j-1
            kk = k 
            jl = j
            kl = k 
            j2 = j+1
            k2 = k  
            io = 1  
        endif
	      t0(jj,kk) = t0(j,k)
C Open boundary.
        if(itype(j,k,i) .eq. 0)then	
          dh(jj,kk,io) = dh(j,k,io)
          du(jj,kk,io) = du(j,k,io)
          dv(jj,kk,io) = dv(j,k,io)	                               
          hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
          udum(jj,kk) = 2.D0*udum(j,k) - udum(j2,k2)
          vdum(jj,kk) = 2.D0*vdum(j,k) - vdum(j2,k2)
C Wall boundary.
        elseif(itype(j,k,i) .eq. 1)then
          dh(jj,kk,io) = dh(j,k,io) 
          du(jj,kk,io) = du(j,k,io)
          dv(jj,kk,io) = dv(j,k,io)	                           
          hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
          udum(jj,kk) = udum(j,k)*(sn(jl,kl,io)*sn(jl,kl,io) - 
     &                     cn(jl,kl,io)*cn(jl,kl,io)) -
     &                      2.D0*vdum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
          vdum(jj,kk) = vdum(j,k)*(cn(jl,kl,io)*cn(jl,kl,io) - 
     &                      sn(jl,kl,io)*sn(jl,kl,io)) -
     &                      2.D0*udum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
C Specified depth and velocity (supercritical).
        elseif(itype(j,k,i) .eq. 2)then
          dh(jj,kk,io) = 0.D0
          du(jj,kk,io) = 0.D0
          dv(jj,kk,io) = 0.D0                 
          hdum(jj,kk) = fix(j,k,1) 
          udum(jj,kk) = fix(j,k,2)
          vdum(jj,kk) = fix(j,k,3)
	        if(isurf .eq. 1) hdum(jj,kk) = hdum(jj,kk) - zc(j,k)
C Tidal (subcritical).
        elseif(itype(j,k,i) .eq. 3)then
          dh(jj,kk,io) = 0.d0
          du(jj,kk,io) = 0.d0
          dv(jj,kk,io) = 0.d0	                               
          hdum(jj,kk) = fix(j,k,1) + amp*dcos(-6.283185D0*t/period)
          udum(jj,kk) = udum(j,k) 
          vdum(jj,kk) = vdum(j,k) 
	        if(isurf .eq. 1) hdum(jj,kk) = hdum(jj,kk) - zc(j,k)
C Specified flow rate (subcritical).
        elseif(itype(j,k,i) .eq. 4)then 
          du(jj,kk,io) =  0.d0
          dv(jj,kk,io) =  0.d0 
          if(hdum(j,k) .ge. epsh) then
            hdum(jj,kk) = hdum(j,k) 
            dh(jj,kk,io) =  0.d0
          elseif(hdum(j,k)*hdum(j2,k2) .ge. epsh) then
            hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
            dh(jj,kk,io) = dh(j,k,io)
          endif
          if(hdum(jj,kk) .ge. epsh) then
           udum(jj,kk) = fix(j,k,2)/hdum(jj,kk)
           vdum(jj,kk) = fix(j,k,3)/hdum(jj,kk)
          else
           udum(jj,kk) = 0.D0
           vdum(jj,kk) = 0.D0
	        endif 
C Solitary wave profile (subcritical).
        elseif(itype(j,k,i) .eq. 5)then 
             dh(jj,kk,io) = 0.D0
             du(jj,kk,io) = 0.d0
             dv(jj,kk,io) = 0.d0 
    	       dum = period - t
    	       if(dum .lt. -36.D0) dum = -36.D0	 
    	       hdum(jj,kk) = fix(j,k,1) + 2.D0*amp/(dcosh(dum)*dcosh(dum))
    	       udum(jj,kk) = u(j,k) 
    	       vdum(jj,kk) = v(j,k) 
    	       if(isurf .eq. 1)  hdum(jj,kk) = hdum(jj,kk) - zc(j,k)
        endif
    	  if(hdum(jj,kk) .lt. 0.D0) hdum(jj,kk) = 0.D0
    	  if(hdum(jj,kk) .lt. epsh) then
    	    udum(jj,kk) = 0.D0
    	    vdum(jj,kk) = 0.D0
	    endif
      enddo

      return
      end
***********************************************************************
      subroutine grid
      include 'dry.inc'
C Read grid data from file 'coords'.
      read(2,*) np, ne
      do i=1,np
        read(2,*) x(i), y(i), z(i)
C         write(1,*) x(i), y(i), z(i)
      enddo
      do j=1,ncol
        do k=kbeg(j),kend(j)          
          read(2,*) (nop(j,k,i), i=1,4)
        enddo
      enddo
C Compute grid metrics.
      do j=1,ncol
      do k=kbeg(j),kend(j)
        n1 = nop(j,k,1)
        n2 = nop(j,k,2)
        n3 = nop(j,k,3)
        n4 = nop(j,k,4)
        xc(j,k) = 0.25D0*(x(n1) + x(n2) + x(n3) + x(n4))
        yc(j,k) = 0.25D0*(y(n1) + y(n2) + y(n3) + y(n4))
        dxdxi =  0.5D0*(-x(n1) + x(n2) + x(n3) - x(n4))
        dxdeta = 0.5D0*(-x(n1) - x(n2) + x(n3) + x(n4))
        dydxi =  0.5D0*(-y(n1) + y(n2) + y(n3) - y(n4))
        dydeta = 0.5D0*(-y(n1) - y(n2) + y(n3) + y(n4))
	      area(j,k) = dxdxi*dydeta - dxdeta*dydxi
        if(area(j,k) .le. 0.D0)then
          write(*,*) 'area error in cell ',j,k
          stop
        endif
	      dxi(j,k,1) = dydeta/area(j,k)
	      deta(j,k,1) = -dydxi/area(j,k)
        dxi(j,k,2) = -dxdeta/area(j,k)
        deta(j,k,2) = dxdxi/area(j,k)
       sx(j,k)=((z(n2)-z(n4))*(y(n3)-y(n1))-(z(n3)-z(n1))*(y(n2)-y(n4)))
     &   /(2.D0*area(j,k))
       sy(j,k)=((z(n3)-z(n1))*(x(n2)-x(n4))-(z(n2)-z(n4))*(x(n3)-x(n1)))
     &   /(2.D0*area(j,k))
        zc(j,k) = 0.25D0*(z(n1) + z(n2) + z(n3) + z(n4))
	      dz(j,k,1) = sx(j,k)*dxdxi + sy(j,k)*dydxi
        dz(j,k,2) = sx(j,k)*dxdeta + sy(j,k)*dydeta
      enddo
      enddo
C Compute cell face angles.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       ddx = x(nop(j,k,2)) - x(nop(j,k,1))
       ddy = y(nop(j,k,2)) - y(nop(j,k,1))
       ds(j,k,2) = dsqrt(ddx*ddx + ddy*ddy)
       sn(j,k,2) = ddx/ds(j,k,2)
       cn(j,k,2) = -ddy/ds(j,k,2)          ! Horizontal face.
       ddx = x(nop(j,k,4)) - x(nop(j,k,1))
       ddy = y(nop(j,k,4)) - y(nop(j,k,1))
       ds(j,k,1) = dsqrt(ddx*ddx + ddy*ddy)
       sn(j,k,1) = -ddx/ds(j,k,1)            ! Vertical face.
       cn(j,k,1) = ddy/ds(j,k,1)
       do i=1,inum(j,k)
        if(ipos(j,k,i) .eq. 3) then
           ddx = x(nop(j,k,3)) - x(nop(j,k,4))
           ddy = y(nop(j,k,3)) - y(nop(j,k,4))
           ds(j,k+1,2) = dsqrt(ddx*ddx + ddy*ddy)
           sn(j,k+1,2) = ddx/ds(j,k+1,2)     ! Top (boundary) faces.
           cn(j,k+1,2) = -ddy/ds(j,k+1,2)
         elseif(ipos(j,k,i) .eq. 2) then
           ddx = x(nop(j,k,3)) - x(nop(j,k,2))
           ddy = y(nop(j,k,3)) - y(nop(j,k,2))
           ds(j+1,k,1) = dsqrt(ddx*ddx + ddy*ddy)
           sn(j+1,k,1) = -ddx/ds(j+1,k,1)     ! Right (boundary) faces.
           cn(j+1,k,1) = ddy/ds(j+1,k,1)
        endif
       enddo
      enddo
      enddo
C Set some things in ghost cells.
      do j=1,ncol
      do k=kbeg(j),kend(j) 
        do i=1, inum(j,k)
          call findbc(i,j,k,jj,kk,j2,k2)
          area(jj,kk) = area(j,k)
          sx(jj,kk) = sx(j,k)
          sy(jj,kk) = sy(j,k)
          dxi(jj,kk,1) =dxi(j,k,1)
          deta(jj,kk,1) = deta(j,k,1)
          dxi(jj,kk,2) = dxi(j,k,2)
          deta(jj,kk,2) = deta(j,k,2)
          xc(jj,kk) = 2.D0*xc(j,k) - xc(j2,k2)
          yc(jj,kk) = 2.D0*yc(j,k) - yc(j2,k2)
          zc(jj,kk) = 2.D0*zc(j,k) - zc(j2,k2)
         enddo
      enddo
      enddo

      return
      end 
************************************************************************
      subroutine findbc(i,j,k,jj,kk,j2,k2)
      include 'dry.inc'
      if(ipos(j,k,i) .eq. 1) then  
         jj = j
         kk = k-1
	       j2 = j
         k2 = k+1
       elseif(ipos(j,k,i) .eq. 2) then
         jj = j+1
         kk = k 
         j2 = j-1
         k2 = k
       elseif(ipos(j,k,i) .eq. 3) then 
         jj = j
         kk = k+1
         j2 = j
         k2 = k-1
       elseif(ipos(j,k,i) .eq. 4) then 
         jj = j-1
         kk = k
         j2 = j+1
         k2 = k 
       endif
      return
      end   
************************************************************************
      subroutine myinterp
      include 'dry.inc'
    	do j=1,ncol
    	do k=kbeg(j), kend(j)
C Extrapolate to corner cells.
    	 ibc1 = 0
    	 ibc2 = 0
    	 ibc3 = 0
    	 ibc4 = 0
    	 do i=1, inum(j,k)
    	  if(ipos(j,k,i) .eq. 1) ibc1 = 1
    	  if(ipos(j,k,i) .eq. 2) ibc2 = 1
    	  if(ipos(j,k,i) .eq. 3) ibc3 = 1
    	  if(ipos(j,k,i) .eq. 4) ibc4 = 1	  
    	 enddo
       if(ibc1 .eq. 1 .and. ibc4 .eq. 1)then
         h(j-1,k-1) = 2.D0*h(j,k) - h(j+1,k+1)
         u(j-1,k-1) = 2.D0*u(j,k) - u(j+1,k+1)
         v(j-1,k-1) = 2.D0*v(j,k) - v(j+1,k+1)
       elseif(ibc1 .eq. 1 .and. ibc2 .eq. 1)then
         h(j+1,k-1) = 2.D0*h(j,k) - h(j-1,k+1)
         u(j+1,k-1) = 2.D0*u(j,k) - u(j-1,k+1)
         v(j+1,k-1) = 2.D0*v(j,k) - v(j-1,k+1)
       elseif(ibc2 .eq. 1 .and. ibc3 .eq. 1)then
         h(j+1,k+1) = 2.D0*h(j,k) - h(j-1,k-1)
         u(j+1,k+1) = 2.D0*u(j,k) - u(j-1,k-1)
         v(j+1,k+1) = 2.D0*v(j,k) - v(j-1,k-1)
       elseif(ibc3 .eq. 1 .and. ibc4 .eq. 1)then
         h(j-1,k+1) = 2.D0*h(j,k) - h(j+1,k-1)
         u(j-1,k+1) = 2.D0*u(j,k) - u(j+1,k-1)
         v(j-1,k+1) = 2.D0*v(j,k) - v(j+1,k-1)
       endif
C           call bconds(j,k,h,u,v)
        enddo
        enddo 
C Interpolate points for plotting.

      return
      end
************************************************************************
      subroutine interp
      include 'dry.inc'

	do j=1,ncol
	do k=kbeg(j), kend(j)
C Extrapolate to corner cells.
	 ibc1 = 0
	 ibc2 = 0
	 ibc3 = 0
	 ibc4 = 0
	 do i=1, inum(j,k)
	  if(ipos(j,k,i) .eq. 1) ibc1 = 1
	  if(ipos(j,k,i) .eq. 2) ibc2 = 1
	  if(ipos(j,k,i) .eq. 3) ibc3 = 1
	  if(ipos(j,k,i) .eq. 4) ibc4 = 1	  
	 enddo
       if(ibc1 .eq. 1 .and. ibc4 .eq. 1)then
         xc(j-1,k-1) = 2.D0*xc(j,k) - xc(j+1,k+1)
         yc(j-1,k-1) = 2.D0*yc(j,k) - yc(j+1,k+1)
         h(j-1,k-1) = 2.D0*h(j,k) - h(j+1,k+1)
         u(j-1,k-1) = 2.D0*u(j,k) - u(j+1,k+1)
         v(j-1,k-1) = 2.D0*v(j,k) - v(j+1,k+1)
       elseif(ibc1 .eq. 1 .and. ibc2 .eq. 1)then
         xc(j+1,k-1) = 2.D0*xc(j,k) - xc(j-1,k+1)
         yc(j+1,k-1) = 2.D0*yc(j,k) - yc(j-1,k+1)
         h(j+1,k-1) = 2.D0*h(j,k) - h(j-1,k+1)
         u(j+1,k-1) = 2.D0*u(j,k) - u(j-1,k+1)
         v(j+1,k-1) = 2.D0*v(j,k) - v(j-1,k+1)
       elseif(ibc2 .eq. 1 .and. ibc3 .eq. 1)then
         xc(j+1,k+1) = 2.D0*xc(j,k) - xc(j-1,k-1)
         yc(j+1,k+1) = 2.D0*yc(j,k) - yc(j-1,k-1)
         h(j+1,k+1) = 2.D0*h(j,k) - h(j-1,k-1)
         u(j+1,k+1) = 2.D0*u(j,k) - u(j-1,k-1)
         v(j+1,k+1) = 2.D0*v(j,k) - v(j-1,k-1)
       elseif(ibc3 .eq. 1 .and. ibc4 .eq. 1)then
         xc(j-1,k+1) = 2.D0*xc(j,k) - xc(j+1,k-1)
         yc(j-1,k+1) = 2.D0*yc(j,k) - yc(j+1,k-1)
         h(j-1,k+1) = 2.D0*h(j,k) - h(j+1,k-1)
         u(j-1,k+1) = 2.D0*u(j,k) - u(j+1,k-1)
         v(j-1,k+1) = 2.D0*v(j,k) - v(j+1,k-1)
       endif
       call bconds(j,k,h,u,v)
      enddo
      enddo
C Interpolate points for plotting.
      itest = 0
      do i=1,np
       call sample(1,i,x(i),y(i),hn(i),un(i),vn(i))
      enddo
	                 
      return
      end
************************************************************************
      subroutine input
      include 'dry.inc' 
      character*72 dum, filename(100)        
C Read input from file 'dry.dat'.
      read(3,'(a72)') dum
      read(3,*) grav, dt, tmax, xsplit, xn
	if(grav .ge. 32.1D0) xn = xn/(3.28D0)**(1.D0/3.D0)
      read(3,'(a72)') dum
      read(3,*) epsh, beta
      read(3,'(a72)') dum
      read(3,*) xk, ainflt, binflt, tc, cappa
      read(3,'(a72)') dum
      read(3,*) istart, imass, ifront, nprt
      read(3,'(a72)') dum
      read(3,*) nbcell
      read(3,'(a72)') dum
      do ii=1,nbcell
         read(3,*)j,k,inum(j,k),
     &   (itype(j,k,i),i=1,inum(j,k)),(ipos(j,k,i),i=1,inum(j,k))
      enddo
      read(3,'(a72)') dum
      read(3,*) ncol
      read(3,'(a72)') dum
      do j=1,ncol       
         read(3,*) idum, kbeg(j), kend(j)
      enddo
      kbeg(ncol+1) = kbeg(ncol)
      kend(ncol+1) = kend(ncol)
	write(*,*) '  grid setup is complete'
      read(3,'(a72)') dum
      read(3,*) h0l, u0l, v0l
      read(3,'(a72)') dum
      read(3,*) h0r, u0r, v0r
      read(3,'(a72)') dum
      read(3,*) ndir
      read(3,'(a72)') dum
      do i=1,ndir
       read(3,*) j,k,fix(j,k,1),fix(j,k,2),fix(j,k,3),period,amp
	    if(period .eq. 0.D0) period = 1.D0
      enddo
      isurf = 2
C       write(*,*) '  choose a scheme'
C       write(*,*) '      0 = Upwind'
C       write(*,*) '      1 = Lax-Wendroff'
C       write(*,*) '      2 = Beam-Warming'
C       write(*,*) '      3 = Fromm'
C       write(*,*) '      4 = Double Minmod'
C       write(*,*) '      5 = Beta Family'
C       read (*,*) ilim
      ilim = 5
      write(*,*) ' '  

C Set up computational grid.
      call grid
C Set initial conditions.
      t = 0.D0 
c     dt = tmax/dfloat(nt)
      do j=1,ncol
        do k=kbeg(j),kend(j)
          n1 = nop(j,k,1)
          n2 = nop(j,k,2)
          n3 = nop(j,k,3)
          n4 = nop(j,k,4)
          if(xc(j,k) .le. xsplit)then
            if(isurf .eq. 1) then
  		        h(j,k) = h0l - zc(j,k)
  	        else
              h(j,k) = h0l
  	        endif
              u(j,k) = u0l
              v(j,k) = v0l
          else
            if(isurf .eq. 1) then
  		        h(j,k) = h0r - zc(j,k)
            else
              h(j,k) = h0r
            endif
             u(j,k) = u0r
             v(j,k) = v0r
          endif
  	      if(xk .gt. 0.D0 .and. h(j,k) .eq. 0.D0) t0(j,k) = t 
        enddo
      enddo     
      do j=1,ncol
        do k=kbeg(j),kend(j)
          if(h(j,k) .lt. 0.D0) h(j,k) = 0.D0
          q(j,k,1) = h(j,k)
          q(j,k,2) = h(j,k)*u(j,k)
          q(j,k,3) = h(j,k)*v(j,k)
C For fixed flux BCs, must initialize depth.          
          do i=1,inum(j,k)
            if(itype(j,k,i) .eq. 4) then
              call findbc(i,j,k,jj,kk,j2,k2)
	            if (ipos(j,k,i) .eq. 1 .or. ipos(j,k,i) .eq. 3) then
	              qflux = fix(j,k,2)*cn(j,k,2) + fix(j,k,3)*sn(j,k,2)
                dx = -dxi(j,k,2)*area(j,k)
                dy = dxi(j,k,1)*area(j,k)
	              dss = dsqrt(dx*dx + dy*dy)
                if(h(j,k) .lt. epsh) then
                  h(jj,kk) = 0.4D0
                  hp(jj,kk) = h(jj,kk)     
                endif
                qflux = dabs(qflux)
                tclip2 = (xk*dss/(qflux*(1.d0-cappa)))**
     &                      (1.d0/(1.d0-ainflt))
                if(tclip2 .gt. tclip) tclip = tclip2
              endif
            endif
C For specified free surface BCs.
	        enddo
        enddo
      enddo
C For problems with bed seepage.
      if(xk .gt. 0.d0) then
       if(tclip .gt. 0.d0) then
          zslope = xk*tclip**(ainflt - 1.d0)
          write(*,*) ' '
          write(*,*) '*************************************************'
          write(*,*) 'Kostiakov infiltration formula modified'
          write(*,*) 'fot time less than ', tclip
          write(*,*) 'infiltration slope is ', zslope
          write(*,*) '*************************************************'
         else
          tclip = 0.d0
       endif
      endif
C       itest = 0
      write(*,*) ' '
	    write(*,*) 'initial conditions are set'
C       pause 'hit return to continue'

      return
      end
************************************************************************
      subroutine sample(iflag,i,xdum,ydum,hdum,udum,vdum)
      include 'dry.inc' 
C This subroutine interpolates h, u, and v at 
C user selected points and prints them to output files.
C Compute interpolation weights first time through only.
      if(itest .eq. 0) then
C Bracket sampling location.
       do j=1,ncol
       do k=kbeg(j),kend(j)+1
	   if(xdum .ge. xc(j-1,k) .and. xdum .lt. xc(j,k)) then
	     if(ydum .ge. yc(j,k-1) .and. ydum .le. yc(j,k)) then	
	      j1(i) = j
	      k1(i) = k
	      goto 10
	     endif
	    elseif(xdum .ge. xc(j,k) .and. xdum .le. xc(j+1,k)) then
	     if(ydum .ge. yc(j,k-1) .and. ydum .le. yc(j,k)) then	
	      j1(i) = j+1
	      k1(i) = k
	      goto 10
	     endif
         endif
	 enddo
	 enddo
C Compute interpolation distances.
10     d1 = dsqrt((xdum - xc(j1(i),k1(i)))**2 + 
     &            (ydum - yc(j1(i),k1(i)))**2)
       d2 = dsqrt((xdum - xc(j1(i)-1,k1(i)))**2 + 
     &            (ydum - yc(j1(i)-1,k1(i)))**2)
       d3 = dsqrt((xdum - xc(j1(i)-1,k1(i)-1))**2 + 
     &            (ydum - yc(j1(i)-1,k1(i)-1))**2)
       d4 = dsqrt((xdum - xc(j1(i),k1(i)-1))**2 + 
     &            (ydum - yc(j1(i),k1(i)-1))**2)
C Compute weights as distance inverses.
C Tecplot states that the 3.5 exponent yields the smoothest results.
	  if(d1 .gt. 0.D0) then
	    w1(i) = d1**(-3.5D0)
	   else
	    w1(i) = 1.D0
	    w2(i) = 0.D0
	    w3(i) = 0.D0
	    w4(i) = 0.D0
	    goto 4
	  endif
	  if(d2 .gt. 0.D0) then
	    w2(i) = d2**(-3.5D0)
	   else
	    w1(i) = 0.D0
	    w2(i) = 1.D0
	    w3(i) = 0.D0
	    w4(i) = 0.D0
	    goto 4
	  endif
	  if(d3 .gt. 0.D0) then
	    w3(i) = d3**(-3.5D0)
	   else
	    w1(i) = 0.D0
	    w2(i) = 0.D0
	    w3(i) = 1.D0
	    w4(i) = 0.D0
	    goto 4
	  endif
	  if(d4 .gt. 0.D0) then
	    w4(i) = d4**(-3.5D0)
	   else
	    w1(i) = 0.D0
	    w2(i) = 0.D0
	    w3(i) = 0.D0
	    w4(i) = 1.D0
	    goto 4
	  endif
4	  sumw = w1(i) + w2(i) + w3(i) + w4(i)
      endif
C Interpolate data.
  	  sumw = w1(i) + w2(i) + w3(i) + w4(i)
	    hdum = (w1(i)*h(j1(i),k1(i)) + w2(i)*h(j1(i)-1,k1(i)) + 
     &        w3(i)*h(j1(i)-1,k1(i)-1) + w4(i)*h(j1(i),k1(i)-1))/sumw
	    udum = (w1(i)*u(j1(i),k1(i)) + w2(i)*u(j1(i)-1,k1(i)) + 
     &        w3(i)*u(j1(i)-1,k1(i)-1) + w4(i)*u(j1(i),k1(i)-1))/sumw
	    vdum = (w1(i)*v(j1(i),k1(i)) + w2(i)*v(j1(i)-1,k1(i)) + 
     &        w3(i)*v(j1(i)-1,k1(i)-1) + w4(i)*v(j1(i),k1(i)-1))/sumw
	    zdum = (w1(i)*zc(j1(i),k1(i)) + w2(i)*zc(j1(i)-1,k1(i)) + 
     &        w3(i)*zc(j1(i)-1,k1(i)-1) + w4(i)*zc(j1(i),k1(i)-1))/sumw
C     if(iflag .lt. 1) then
C     idev = 13 + i
c	    write(idev,99) i, t, hdum, udum, vdum, hdum + zdum
C       endif

      return
 99   format(' ',i5,5e15.7)
	end
************************************************************************
      subroutine output
      include 'dry.inc' 
C Write to file 'restart.out'.
      do j=1,ncol
      do k=kbeg(j),kend(j)
        write(11,199) j, k, h(j,k), u(j,k), v(j,k), t0(j,k)
      enddo
      enddo
C Push cell-centered values to corner nodes for plotting.
	    call interp
C Write to file 'dry.plt'.
      write(*,*) ' '
      write(*,*) 'OUTPUT OPTIONS:'
      write(*,*) 'enter 1 to plot depth'
	    write(*,*) '      2 to plot bed'
	    write(*,*) '      3 to plot free surface'
	read(*,*) iout
	if(iout .eq. 1) then
	write(8,*) 'VARIABLES = "X", "Y", "H", "U", "V", "Z", "Fr"'
      write(8,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT', ' ET=QUADRILATERAL'
      do n=1,np
	 if(hn(n) .gt. 0.D0) then
	  fr = dsqrt((un(n)*un(n) + vn(n)*vn(n))/(grav*hn(n)))
	  else
	  fr = 0.D0
	 endif
       write(8,200) x(n), y(n), hn(n), un(n), vn(n), z(n), fr
      enddo
	elseif(iout .eq. 2) then
	write(8,*) 'VARIABLES = "X", "Y", "Z", "U", "V", "H", "Fr"'
      write(8,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT', ' ET=QUADRILATERAL'
      do n=1,np
	 if(hn(n) .gt. 0.D0) then
	  fr = dsqrt((un(n)*un(n) + vn(n)*vn(n))/(grav*hn(n)))
	  else
	  fr = 0.D0
	 endif
       write(8,200) x(n), y(n), z(n), un(n), vn(n), hn(n), fr
      enddo
	elseif(iout .eq. 3) then
	write(8,*) 'VARIABLES = "X", "Y", "H+Z", "U", "V", "H", "Z", "Fr"'
      write(8,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT', ' ET=QUADRILATERAL'
      do n=1,np
	 if(hn(n) .gt. 0.D0) then
	  fr = dsqrt((un(n)*un(n) + vn(n)*vn(n))/(grav*hn(n)))
	  else
	  fr = 0.D0
	 endif
       write(8,200) x(n),y(n),hn(n)+z(n),un(n),vn(n),hn(n),z(n),fr
      enddo
	endif
      do j=1,ncol
      do k=kbeg(j),kend(j)
	 write(8,201) (nop(j,k,i), i=1,4)
	enddo
	enddo
	write(*,*) 'enter a 1 to extract a slice of output'
	read(*,*) iext
	if(iext .eq. 1) then
	 write(*,*) 'enter 0 for constant x slice'
	 write(*,*) 'else enter 1 for constant y slice'
	 read(*,*) islice
	 if(islice .eq. 0) then
	  write(*,*) 'enter x value of slice'
	  read(*,*) xslice
        do n=1,np
	   if(x(n) .eq. xslice) then
	    pdum = un(n)*hn(n)
	    qdum = vn(n)*hn(n)
	    if(hn(n) .gt. 0.D0) 
     &     fr = dsqrt(un(n)**2+vn(n)**2)/dsqrt(grav*hn(n))
          write(13,200) y(n), hn(n), un(n), vn(n), pdum, qdum, fr, z(n)
	   endif
        enddo
	 else
	  write(*,*) 'enter y value of slice'
	  read(*,*) yslice
        do n=1,np
	   if(y(n) .eq. yslice) then
	    pdum = un(n)*hn(n)
	    qdum = vn(n)*hn(n)
	    eta = hn(n) + z(n)
	    if(hn(n) .gt. 0.D0) 
     &     fr = dsqrt(un(n)**2+vn(n)**2)/dsqrt(grav*hn(n))
          write(13,200) x(n),hn(n),un(n),vn(n),pdum,qdum,fr,eta,z(n)
	   endif
        enddo
	 endif
	endif

      return
 199  format(' ',2i6,4e15.7)
 200  format(' ',9e15.5)
 201  format(' ',4i6)
      end
************************************************************************
      subroutine movie
      include 'dry.inc'
      do i=1,np
	xdum=x(i)
	ydum=y(i)
C Bracket sampling location.
       do j=1,ncol
       do k=kbeg(j),kend(j)+1
	   if(xdum .ge. xc(j-1,k) .and. xdum .lt. xc(j,k)) then
	     if(ydum .ge. yc(j,k-1) .and. ydum .le. yc(j,k)) then	
	      j1(i) = j
	      k1(i) = k
	      goto 10
	     endif
	    elseif(xdum .ge. xc(j,k) .and. xdum .le. xc(j+1,k)) then
	     if(ydum .ge. yc(j,k-1) .and. ydum .le. yc(j,k)) then	
	      j1(i) = j+1
	      k1(i) = k
	      goto 10
	     endif
         endif
	 enddo
	 enddo
C Compute interpolation distances.
10     d1 = dsqrt((xdum - xc(j1(i),k1(i)))**2 + 
     &            (ydum - yc(j1(i),k1(i)))**2)
       d2 = dsqrt((xdum - xc(j1(i)-1,k1(i)))**2 + 
     &            (ydum - yc(j1(i)-1,k1(i)))**2)
       d3 = dsqrt((xdum - xc(j1(i)-1,k1(i)-1))**2 + 
     &            (ydum - yc(j1(i)-1,k1(i)-1))**2)
       d4 = dsqrt((xdum - xc(j1(i),k1(i)-1))**2 + 
     &            (ydum - yc(j1(i),k1(i)-1))**2)
C Compute weights as distance inverses.
C Tecplot states that the 3.5 exponent yields the smoothest results.
	  if(d1 .gt. 0.D0) then
	    w1(i) = d1**(-3.5D0)
	   else
	    w1(i) = 1.D0
	    w2(i) = 0.D0
	    w3(i) = 0.D0
	    w4(i) = 0.D0
	    goto 4
	  endif
	  if(d2 .gt. 0.D0) then
	    w2(i) = d2**(-3.5D0)
	   else
	    w1(i) = 0.D0
	    w2(i) = 1.D0
	    w3(i) = 0.D0
	    w4(i) = 0.D0
	    goto 4
	  endif
	  if(d3 .gt. 0.D0) then
	    w3(i) = d3**(-3.5D0)
	   else
	    w1(i) = 0.D0
	    w2(i) = 0.D0
	    w3(i) = 1.D0
	    w4(i) = 0.D0
	    goto 4
	  endif
	  if(d4 .gt. 0.D0) then
	    w4(i) = d4**(-3.5D0)
	   else
	    w1(i) = 0.D0
	    w2(i) = 0.D0
	    w3(i) = 0.D0
	    w4(i) = 1.D0
	    goto 4
      endif
4	  sumw = w1(i) + w2(i) + w3(i) + w4(i)
      
C Interpolate data.
  	sumw = w1(i) + w2(i) + w3(i) + w4(i)
	hn(i) = (w1(i)*h(j1(i),k1(i)) + w2(i)*h(j1(i)-1,k1(i)) + 
     &        w3(i)*h(j1(i)-1,k1(i)-1) + w4(i)*h(j1(i),k1(i)-1))/sumw
c	z(i) = (w1(i)*zc(j1(i),k1(i)) + w2(i)*zc(j1(i)-1,k1(i)) + 
c     &        w3(i)*zc(j1(i)-1,k1(i)-1) + w4(i)*zc(j1(i),k1(i)-1))/sumw
	enddo

C Write to file 'movie.dat'.
	if(t.eq.0.)then
	write(9,*) 'VARIABLES = "X", "Y", "Z", "H", "D"'
      write(9,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT',' ET=QUADRILATERAL'
        do n=1,np
          write(9,200) x(n),y(n),z(n),z(n),z(n)
	  enddo
	  do j=1,ncol
          do k=kbeg(j),kend(j)
	      write(9,201) (nop(j,k,i), i=1,4)
	    enddo
	  enddo
	elseif(mod(it,30).eq.0)then
      write(9,*) 'VARIABLES = "X", "Y", "Z", "H", "D"'
      write(9,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT',' ET=QUADRILATERAL'
     &  ,' D=(FECONNECT,1,2,3)'
	  do n=1,np
          write(9,200) hn(n)+z(n),hn(n)
	  enddo
	endif
      return
 200  format(' ',5e15.5)
 201  format(' ',4i6)
      end
************************************************************************
      subroutine limitr(i,beta,dq1,dq2,dq)
      implicit real*8(a-h,o-z)
C Lax-Wendroff.
      if(i .eq. 1)then
        dq = dq2 
C Beam-Warming.
      elseif(i .eq. 2)then
        dq = dq1 
C Fromm
      elseif(i .eq. 3)then
        dq = 0.5D0*(dq1 + dq2) 
C Double minmod.
      elseif(i .eq. 4)then
        a = 0.5D0*(dq1 + dq2)
        b = 2.D0*dq1
        c = 2.D0*dq2
        if(a*b .gt. 0.D0 .and. b*c .gt. 0.D0) then
           dq = fmin1(a, b, c)  
         else
           dq = 0.D0
        endif
C Beta Family.
      else
        if(dq1*dq2 .le. 0.D0) then
           dq = 0.D0
         else
           dq = fmin2(fmax2(dq1, dq2), beta*fmin2(dq1, dq2))
        endif
      endif

      return     
      end
************************************************************************
      real*8 function f1(h,up)
      implicit real*8(a-h,o-z)
      f1 = h*up
      return
      end
************************************************************************
      real*8 function f2(h,u,up,cn)
      implicit real*8(a-h,o-z)
      common/m/ grav, amax, epsh, xn
      f2 = h*u*up + 0.5D0*grav*h*h*cn
      return
      end
************************************************************************
      real*8 function f3(h,v,up,sn)
      implicit real*8(a-h,o-z)
      common/m/ grav, amax, epsh, xn
      f3 = h*v*up + 0.5D0*grav*h*h*sn
      return
      end
************************************************************************
      real*8 function fmin1(a,b,c)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
         fmin1 = -1.D0*dmin1(dabs(a),dabs(b),dabs(c)) 
        else
         fmin1 = dmin1(a,b,c) 
      endif 
      return
      end
************************************************************************
      real*8 function fmin2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
          fmin2 = -dmin1(dabs(a),dabs(b))
         else
          fmin2 = dmin1(a,b)
      endif 
      return
      end
************************************************************************
      real*8 function fmax2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
          fmax2 = -dmax1(dabs(a),dabs(b))
         else
          fmax2 = dmax1(a,b) 
      endif
      return
      end
