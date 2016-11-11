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

C   Added :  flux tracking
C             meta.out -  print dt and hnorm to file for python use
C             code to turn of influx at tc in bconds subroutine
C   Removed:  movie subroutine, ilim = 0, volume and front tracking 
C             unused boundary conditions  (2, 3 and 5)
C             unused subroutines interp, movie, output, sample
************************************************************************
      include 'dry.inc'  
      open(2,file ='coords')
      open(3,file ='dryin.dat')
C 4 output files: i    
	    open(100,file= 'output/h.out')
	    open(101,file='output/time.out')
      open(104,file='output/fluxes.out')
      open(105,file='output/allfluxes.out')
      open(120,file='output/dvol.out')                  
      open(102,file='output/meta.out')                  
      call cpu_time(start)       
C Read input, set up grid, initial conditions.
      call input
      if(ifront .eq. 1 .or. imass .eq. 1) open(12,file='diag.out')
	    write(101, *) "time  |  print step  |  time step"            
      iprt = 0 
C Begin time loop.
      do while (t.lt.tmax)
	      it=it+1
        if(t.gt.60.)   dt=0.2d0
        if (t.gt.60. .and. t.lt.60. + dt) write(102, *) 'dt = ', dt
        if(t.gt.6960.) dt=0.4d0
c       nprt = tmax/100/dt                
        t = t + dt
        amax = 0.D0
C Compute predictor.
        do j=1,ncol
        do k=kbeg(j),kend(j)
         call bconds(j,k,h,u,v)
         call predict(j,k)
        enddo
        enddo
C Loop over cells to compute fluxes.
        flux1 = 0.d0
        flux2 = 0.d0
        flux3 = 0.d0
        flux4 = 0.d0
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
              flux2 = flux2  - f(j+1,k,1, 1)*ds(j+1,k,1)*dt
            elseif(ipos(j,k,i) .eq. 1) then                 ! lower boundary
              flux1 = flux1 + f(j,k,1, 2)*ds(j,k,2)*dt
            elseif(ipos(j,k,i) .eq. 4) then                 ! left boundary
              flux4 = flux4 + f(j,k,1,1)*ds(j,k,1)*dt                
            endif
          enddo
        enddo
        enddo
C         write(105, *) flux1, flux2, flux3, flux4
C Compute corrector solution.
        zinfl = 0.d0
        do j=1,ncol
        do k=kbeg(j),kend(j)
          call source(j,k,hp(j,k),up(j,k),vp(j,k))
          zinfl = zinfl + dt*qs(1)*area(j,k)
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
C     Record infiltration and volume change
        vol0 = vol
        vol = 0.D0
        do j=1,ncol 
          do k=kbeg(j),kend(j)
            call source(j,k,hp(j,k),up(j,k),vp(j,k))
            vol = vol + h(j,k)*area(j,k)
          enddo
        enddo 
        dvol = vol - vol0
        write(105,200) flux1, flux2, flux3, flux4
        flux =  flux1 + flux2 + flux3 + flux4     
        write(120,200) t, dvol, flux, zinfl        
        
C   Below here only executed every nprt time steps        
       iprt = iprt + 1
       if(iprt .eq. nprt) then
	      iprt = 0
        itp = itp +  1
         call myoutput         
        write(*,201) t, amax*dt
	     endif    

      enddo
      
      call cpu_time(finish)
      print '("Time = ",f7.3," seconds.")',finish-start
  
      stop
 200  format(' ',f8.1, 5e19.8)
 201  format(' ','time is ',f8.1,' maximum CFL number is ',f7.4 )
      end
************************************************************************
      subroutine myoutput
      include 'dry.inc'
  
C 	  file 101 is 'time.out'  -  to keep track of the time stpes
      write(101, 203)   t, itp, it
C 	  file 100 is 'myout.out' 
      do j=1,ncol
        do k=kbeg(j),kend(j)
          write(100, *) j, k, h(j,k)
        enddo
      enddo
	    write(100, 202)  itp

C    Loop over cells to compute fluxes.  file 104 is 'fluxes.out'
      do j=1,ncol
        do k=kbeg(j),kend(j)
          do i=1,inum(j,k)
            if( ipos(j,k,i) .eq. 3 .and. itype(j,k,i) .eq. 4 ) then   ! horizontal boundaries
C                   write(104, *) j,k, f(j,k,1, 2)
                write(104, *) j,k,f(j,k+1,1, 2)*ds(j,k+1, 2)
            endif
          enddo
        enddo
      enddo
      write(104, 202)  itp, t
C       write similar for all boundaries
      return
 202  format(' ', i8, f9.2)
 203  format(' ', f7.2 , 4i6, 4i6 )
      end
************************************************************************
      subroutine source(j,k,hdum,udum,vdum)
      include 'dry.inc'

      if(hdum .ge. epsh) then
        call picard(j,k,hdum, winflt)
!         tnew = dmax1(0.d0, t - t0(j,k))
!         if(tnew .le. 0.d0) then
!            znew = 0.d0
!           elseif(tnew .gt. 0.d0 .and. tnew .le. tclip) then
!            znew = zslope*tnew
!           elseif(tnew .gt. tclip .and. tnew .le. tc) then
!            znew = xk*tnew**ainflt
!           elseif(tnew .gt. tc) then
!            znew = xk*tc**ainflt + binflt*(tnew - tc)
!         endif
!         told = dmax1(0.d0, t - t0(j,k) - dt)
!         if(told .le. 0.d0) then
!           zold = 0.d0
!         elseif(told .gt. 0.d0 .and. told .le. tclip) then
!           zold = zslope*told
!         elseif(told .gt. tclip .and. told .le. tc) then
!           zold = xk*told**ainflt
!         elseif(told .gt. tc) then
!            zold = xk*tc**ainflt + binflt*(told - tc)
!         endif
!         winflt = (zold - znew)/dt
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
      subroutine picard(j,k,hdum, winflt)
      include 'dry.inc'
      
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
        enddo
      else 
        call solver(hl,hr,ul,ur,vl,vr,fdum,snn,cnn,dx,dy)
	      do i=1,3
          f(jr,kr,i,i1) = fdum(i)
        enddo
      endif

	    return
      end 
************************************************************************
      subroutine solver(hl,hr,ul,ur,vl,vr,f,sn,cn,dx,dy)
      implicit real*8(a-h,o-z)
      dimension ws(3), e(3,3), a(3), astar(3), da(3), f(3), dum(3)
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
          dum(i) = dum(i) + (astar(l)*ws(l))*e(i,l) 
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
      up(j,k) = u(j,k) - 0.5D0*dt*(
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
C   Loop over all boundary faces in the cell. 
      do i=1, inum(j,k)
        if(ipos(j,k,i) .eq. 1) then! front face.
            jj = j
            kk = k-1
            jl = j
            kl = k 
            j2 = j
            k2 = k+1
            io = 2 
        elseif(ipos(j,k,i) .eq. 2) then! right face.
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
     &                                 cn(jl,kl,io)*cn(jl,kl,io)) -
     &                   2.D0*vdum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
          vdum(jj,kk) = vdum(j,k)*(cn(jl,kl,io)*cn(jl,kl,io) - 
     &                      sn(jl,kl,io)*sn(jl,kl,io)) -
     &                      2.D0*udum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
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
          if (t .gt. tc)  itype(j,k,i) = 1
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
C     Read grid data from file 'coords'.
      read(2,*) np, ne
    	if(np .gt. nn) then
    	  write(*,*) 'ERROR: parameter nn in file dry.inc is too small'
    	  write(*,*) '       for desired number of grid points'
    	  stop
    	endif
      do i=1,np
        read(2,*) x(i), y(i), z(i)
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
        dxdxi = 0.5D0*(-x(n1) + x(n2) + x(n3) - x(n4))
        dxdeta = 0.5D0*(-x(n1) - x(n2) + x(n3) + x(n4))
        dydxi = 0.5D0*(-y(n1) + y(n2) + y(n3) - y(n4))
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
        sx(j,k)=((z(n2)-z(n4))*(y(n3)-y(n1))-(z(n3)-z(n1))
     &     *(y(n2)-y(n4)))/(2.D0*area(j,k))
        sy(j,k)=((z(n3)-z(n1))*(x(n2)-x(n4))-(z(n2)-z(n4))
     &     *(x(n3)-x(n1)))/(2.D0*area(j,k))
        zc(j,k) = 0.25D0*(z(n1) + z(n2) + z(n3) + z(n4))
    	  dz(j,k,1) = sx(j,k)*dxdxi + sy(j,k)*dydxi
    	  dz(j,k,2) = sx(j,k)*dxdeta + sy(j,k)*dydeta
      enddo
      enddo
      write(102,*) "max_slope = " , maxval(sqrt(sx**2+sy**2)*100)/2
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
      subroutine input
      include 'dry.inc' 
!       character*72 dum, filename(100)
C     Read 199rom file 'dryin.dat'.
      read(3,'(a72)') dum
      read(3,*) grav, dt, tmax, xsplit, xn
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
      enddo
C   write(*,*) 'enter a 1 to set free surface elevation'
C   write(*,*) 'enter a 2 to set flow depth'
C   read(*,*) isurf
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
      if(istart .ne. 1)then
        t = 0.D0 
c        dt = tmax/dfloat(nt)
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
      endif 
      do j=1,ncol
      do k=kbeg(j),kend(j)
        if(h(j,k) .lt. 0.D0) h(j,k) = 0.D0
        q(j,k,1) = h(j,k)
        q(j,k,2) = h(j,k)*u(j,k)
        q(j,k,3) = h(j,k)*v(j,k)
        do i=1,inum(j,k)
C For fixed flux BCs, must initialize depth.
          if(itype(j,k,i) .eq. 4) then
            write(*,*) ' '
            write(*,*) 'fixed flux BC is requested in cell ',j,k
            call findbc(i,j,k,jj,kk,j2,k2)
            if(ipos(j,k,i) .eq. 1 .or. ipos(j,k,i) .eq. 3) then
              qflux = fix(j,k,2)*cn(j,k,2) + fix(j,k,3)*sn(j,k,2)
              dx = -dxi(j,k,2)*area(j,k)
              dy = dxi(j,k,1)*area(j,k)
              dss = dsqrt(dx*dx + dy*dy)
              if(h(j,k) .lt. epsh) then
C                 write(*,*) ' *** bed adjacent to the boundary is dry ***'
                if(qflux*dz(j,k,2).lt.0.D0.and.xn.gt.0.d0) then
                  qflux = dabs(qflux)
                  hnorm=(qflux*xn/dsqrt(dabs(dz(j,k,2)/dss)))
     &                     **(3.D0/5.D0)                            
C           write(*,*) ' *** normal depth = ', hnorm,' is specified'
                    h(jj,kk) = hnorm 
                    hp(jj,kk) = hnorm    
	              else
C           write(*,*) 'adverse slope or zero Manning n in cell'
C           write(*,*) 'enter initial flow depth at specified'
C           write(*,*) 'flux boundary ', j,k,ipos(j,k,i)
C           read(*,*) h(jj,kk)
                    h(jj,kk) = (- qflux*xn/dsqrt(dabs(dz(j,k,2)/dss)))
     &                     **(3.D0/5.D0)                      
                    hp(jj,kk) = h(jj,kk)     
	              endif
                write(102, *)  j,k, h(jj,kk)           
              endif
              qflux = dabs(qflux)
              tclip2 = (xk*dss/(qflux*(1.d0-cappa)))
     &            **(1.d0/(1.d0-ainflt))
              if(tclip2 .gt. tclip) tclip = tclip2
            endif
          endif
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
      itest = 0
      write(*,*) ' '
	    write(*,*) 'initial conditions are set'
C       pause 'hit return to continue'
      return
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
