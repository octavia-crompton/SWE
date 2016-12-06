************************************************************************

      parameter (nn=7000,nx=410,ny=110)

************************************************************************

      implicit real*8(a-h,o-z)

      common/cords/ x(nn), y(nn), z(nn), xc(0:nx,0:ny), yc(0:nx,0:ny), 
     &              zc(0:nx,0:ny), dz(0:nx,0:ny,2)
      common/flux/  f(0:nx,0:ny,3,2), q(0:nx,0:ny,3), fdum(3), qs(3) 
      common/grads/ dh(0:nx,0:ny,2), du(0:nx,0:ny,2), dv(0:nx,0:ny,2)
      common/gridi/ itype(nx,ny,4), ipos(nx,ny,4), inum(nx,ny),
     &              nop(nx,ny,4), kbeg(nx), kend(nx)
      common/gridr/ sn(0:nx,0:ny,2), cn(0:nx,0:ny,2), 
     &              ds(0:nx,0:ny,2), dxi(0:nx,0:ny,2),
     &              deta(0:nx,0:ny,2), area(0:nx,0:ny) 
      common/data/  w1(nn), w2(nn), w3(nn), w4(nn), xmon(100), ymon(100)
      common/data1/ j1(nn), k1(nn)
      common/srctrm/ qsource(0:nx,0:ny,3,2), 
     &               fix(nx,ny,3), sx(0:nx,0:ny), sy(0:nx,0:ny), 
     &               t0(0:nx,0:ny)
      common/m/    grav, amax, xn, epsh
      common/ints/ ncol, np, nt, ne, nbcell, ilim, imass, ndum, isource,
     &             ifront, nprt, isurf, nmon, itest, ientropy, it, itp
      common/reals/ tmax, t, dt, tc, xsplit, sw, beta, period,  
     &              amp, xk, ainflt, binflt, tclip, zslope, cappa,
     &              flux1, flux2, flux3, flux4,
     &              vol, vol0, dvol, flux, zinfl, prate
      common/vars/ h(0:nx,0:ny), u(0:nx,0:ny), v(0:nx,0:ny),
     &             hn(nn), un(nn), vn(nn),
     &             hp(0:nx,0:ny), up(0:nx,0:ny), vp(0:nx,0:ny)
    