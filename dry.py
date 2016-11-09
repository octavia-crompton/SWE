import matplotlib.pylab as plt
import numpy as np

#%load_ext autoreload
#%autoreload 2
# from dry import *


plt.ion()
plt.style.use('ggplot')
cmap = 'rainbow'


# noinspection PyAttributeOutsideInit
class dry:
    def __init__(self):
        """
        Define global and class variables
        """
        # <editor-fold desc=" define global constants.">
        global grav, dt, xsplit, xn, epsh, beta, xk
        global ainflt, binflt, tc, cappa, nprt, nbcell, amp, period
        grav = 9.806
        tmax = .15
        xsplit = 100.  # dividing line for dam break problem
        xn = 0.05  # Manning n
        epsh = 0.0025  # depth tolerence for dry bed problems.
        beta = 2.0  #
        xk = 3.9217e-4  # m/s^a where a 0.5
        ainflt = 0.5
        binflt = 2.65e-7
        tc = 6840.
        cappa = 0.99
        amp = 0
        nprt = 10  # print frequency
        nbcell = 386  # number of boundary cells
        period = 1
        # </editor-fold>

        self.input()  # input calls grid
        try:
            self.get_weights()
        except IOError:
            self.transform()
        self.iprt = 0
        self.it = 0
        self.dt = 0.15

        # begin time loop:
        while self.t < tmax:
            self.step()
            self.plot()

    def step(self):
        self.it +=  1
        self.t += self.dt
        if self.t > 60.:
            self.dt = 0.2
        elif self.t > 6960.:
            self.dt = 1.0
        self.movie()
        self.amax = 0.0

        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[nrow + 1] + 1):
                self.h, self.u, self.v = self.bconds(j, k, self.h, self.u, self.v)
                self.predict(j, k)

        # Loop over cells to compute fluxes
        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[nrow + 1] + 1):
                self.hp, self.up, self.vp = self.bconds(j, k, self.hp, self.up, self.vp)
                self.fluxes(j - 1, j, k, k, 0)  # vertical faces
                self.fluxes(j, j, k - 1, k, 1)  # horizontal faces.
                for i in range(inum[j, k]):
                    if ipos[j, k, i] == 3:
                        self.fluxes(j, j, k, k + 1, 1)  # top boundaries
                    elif ipos[j, k, i] == 2:
                        self.fluxes(j, j + 1, k, k, 0)  # right boundaries

        # Compute corrector solution.
        if ilim == 0:
            for j in range(1, ncol + 1):
                for k in range(kbeg[j], kend[j] + 1):
                    for l in range(3):
                        self.q[j, k, l] += self.dt / area[j, k] * (
                            (self.qsource[j, k, l, 1] + self.f[j, k, l, 1]) * ds[j, k, 1] +
                            (self.qsource[j, k, l, 0] + self.f[j, k, l, 0]) * ds[j, k, 0] +
                            (self.qsource[j + 1, k, l, 0] - self.f[j + 1, k, l, 0]) * ds[j + 1, k, 0] +
                            (self.qsource[j, k + 1, l, 1] - self.f[j, k + 1, l, 1]) * ds[j, k + 1, 1])
        else:
            for j in range(1, ncol + 1):
                for k in range(kbeg[j], kend[j] + 1):
                    qs = self.source(j, k, self.hp[j, k], self.up[j, k], self.vp[j, k])
                    for l in range(3):
                        self.q[j, k, l] = self.q[j, k, l] + self.dt / area[j, k] * (
                            self.f[j, k, l, 1] * ds[j, k, 1] + self.f[j, k, l, 0] * ds[j, k, 0] -
                            self.f[j + 1, k, l, 0] * ds[j + 1, k, 0] - self.f[j, k + 1, l, 1] * ds[j, k + 1, 1]) + \
                                          self.dt * qs[l]

        # Store solution.
        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[j] + 1):
                # Check for negative depth. Solve continuity equation in all cells
                if self.q[j, k, 0] > 0.:
                    self.h[j, k] = self.q[j, k, 0]
                else:
                    self.q[j, k, 0] = 0.
                    self.h[j, k] = 0.
                # Neglect momentum in nearly dry cells / solve momentum equations in wet cells only.
                if self.h[j, k] < epsh:
                    self.u[j, k] = 0.
                    self.v[j, k] = 0.
                    self.t0[j, k] = self.t
                    self.q[j, k, 1:] = 0.
                elif self.h[j, k] >= epsh:
                    self.u[j, k] = self.q[j, k, 1] / self.h[j, k]
                    self.v[j, k] = self.q[j, k, 2] / self.h[j, k]

        # If imass or ifront = 1 we'lll open a file called 'diag.out'
        imass = 1
        ifront = 1
        self.iprt += 1
        if self.iprt == nprt:
            self.iprt = 0
            if imass == 1:  # Check volume error.
                vol = 0.
                for j in range(1, ncol + 1):
                    for k in range(kbeg[j], kend[j] + 1):
                        vol = vol + self.h[j, k] * area[j, k]
            if ifront == 1:  # Find location of wetting front for plotting.
                xf = 0.
                for j in range(1, ncol + 1):
                    for k in range(kbeg[j], kend[j] + 1):
                        n2 = nop[j, k, 1]
                        n3 = nop[j, k, 2]
                        xdum = max(x[n2], x[n3])
                        if (self.h[j, k] >= epsh) and (xdum > xf):
                            xf = xdum
            if (ifront == 1) or (imass == 1):
                file = open('diag_py.out', 'a')
                file.write("%f %f %f\n" % (self.t, vol, xf))
                file.close()

        print "time is ", self.t, ' and maximum CFL number is ', self.amax * self.dt

    def input(self):
        """
        Global variables:  constants
                inum : number of boundary faces in a boundary cell.
                inum = 1 for boundary cells with 1 face
                inum = 2 for corner cells
            isurf :  enter a 1 to set free surface elevation
                     enter a 2 to set flow depth
            ilim  : limit scheme
                    0 = Upwind, 1 = Lax-Wendroff,   2 = Beam-Warming
                    3 = Fromm,  4 = Double Minmod,  5 = Beta Family
        Class variables:  not constant
            tclip = 0 :  used to determine dt - timestep interval
        """

        self.istart = 0
        # <editor-fold desc="Read in dryin.dat:">
        dryin = []
        for\
                line in open("dryin.dat", 'r'):
            dryin.append(line)
        # </editor-fold>
        # <editor-fold desc="Get the number of columns and rows">
        global ncol, nrow
        nrow = []
        ncol = []
        for line in range(11, 11 + 386):
            a = dryin[line].strip().split(" ")
            a = [int(i) for i in a if i]
            ncol.append(a[0])
            nrow.append(a[1])
        nrow = len(np.unique(nrow))
        ncol = len(np.unique(ncol))
        # </editor-fold>
        # <editor-fold desc="Read in inum, ipos, itype">
        global inum, ipos, itype
        inum = np.zeros([ncol+1, nrow+1], dtype=int)
        itype = np.zeros([ncol+1, nrow+1, 2])
        ipos = np.zeros([ncol+1, nrow+1, 2])
        for line in range(11, 11 + 386):
            a = dryin[line].strip().split(" ")
            a = [int(x) for x in a if x]
            j = a[0]
            k = a[1]
            inum[j, k] = a[2]
            itype[j, k, 0] = a[3]
            if a[2] == 1:
                ipos[j, k, 0] = a[4]
            elif a[2] > 1:
                itype[j, k, 1] = a[4]
                ipos[j, k, 0] = a[5]
                ipos[j, k, 1] = a[6]
        # </editor-fold>
        # <editor-fold desc="read in kbed, kend">
        global kbeg, kend
        kbeg = [0]
        kend = [0]
        for line in range(400, 399 + ncol):
            # print dryin[line]
            a = dryin[line].strip().split(" ")
            a = [int(x) for x in a if x]
            kbeg.append(a[1])
            kend.append(a[2])
        kbeg.append(kbeg[-1])
        kend.append(kend[-1])
        # </editor-fold>
        # <editor-fold desc="initial conditions for dam break problem">
        global h0l, u0l, v0l, h0r, u0r, v0r
        h0l = 0.
        u0l = 0.
        v0l = 0.
        h0r = 0.
        u0r = 0.
        v0r = 0.
        # </editor-fold>
        # <editor-fold desc="Read in fixed BC cells and monitoring points">
        global fix, nmon, xmon, ymon
        ndir = 3  # number of fixed bc cells
        fix = np.zeros([ncol + 1, nrow + 1, 3])
        for i, line in enumerate(range(501, 501 + ndir)):
            a = dryin[line].strip().split(" ")
            a = [float(b.strip().replace("d", "")) for b in a if b]
            j = int(a[0])
            k = int(a[1])
            fix[j, k, 0] = a[2]
            fix[j, k, 1] = a[3]
            fix[j, k, 2] = a[4]

        nmon = 3  # number of monitoring points
        xmon = np.zeros(nmon)
        ymon = np.zeros(nmon)
        for i, line in enumerate(range(507, 507 + nmon)):
            a = dryin[line].strip().split(" ")
            a = [float(b.strip().replace("d", "")) for b in a if b]
            i = int(a[0] - 1)
            xmon[i] = a[1]
            ymon[i] = a[2]
        # </editor-fold>

        # <editor-fold desc="User sets isurf and ilim">
        global isurf, ilim, sw
        isurf = 2  # float(raw_input("isurf = "))
        ilim = 5  # float(raw_input("ilim = "))
        sw = 0.  # Upwind the source terms if upwind differencing is requested.
        if ilim == 0: sw = 1.0
        if self.istart == 1:  # incomplete code here
            self.t = float(raw_input("enter restart time: "))
        # </editor-fold>

        self.grid()  # call grid function

        # <editor-fold desc="create empty place holder arrays: hp, ">
        self.hp = np.zeros([ncol + 2, nrow + 2])
        self.up = np.zeros([ncol + 2, nrow + 2])
        self.vp = np.zeros([ncol + 2, nrow + 2])
        self.t0 = np.zeros([ncol + 2, nrow + 2])

        self.dh = np.zeros([ncol + 2, nrow + 2, 2])  # modified by conds, predict, and fluxes functions
        self.du = np.zeros([ncol + 2, nrow + 2, 2])  # same
        self.dv = np.zeros([ncol + 2, nrow + 2, 2])  # same
        self.f = np.zeros([ncol + 2, nrow + 2, 3, 2])
        self.qsource = np.zeros([ncol + 2, nrow + 2, 3, 2])  # source term if ilim = 0
        # </editor-fold>
        # <editor-fold desc="Compute initial h, u, v at cell centers if no restart file">
        self.h = np.zeros([ncol + 2, nrow + 2])
        self.u = np.zeros([ncol + 2, nrow + 2])
        self.v = np.zeros([ncol + 2, nrow + 2])
        if self.istart != 1:  # Set the initial conditions
            self.t = 0.
            for j in range(1, ncol + 1):
                for k in range(kbeg[j], kend[j] + 1):
                    if self.xc[j, k] < xsplit:  #  this could be more python friendly
                        if isurf == 1:  # if isurf = 1,  subtract z from h0l and h0r to get h
                            self.h[j, k] = h0l - self.zc[j, k]
                        else:  #
                            self.h[j, k] = h0l
                        self.u[j, k] = u0l
                        self.v[j, k] = v0l
                    else:
                        if isurf == 1:
                            self.h[j, k] = h0r - self.zc[j, k]
                        else:
                            self.h[j, k] = h0r
                        self.u[j, k] = u0r
                        self.v[j, k] = v0r

                    if xk > 0. and self.h[j, k] == 0.:  # set t0 to t in all the cells where h[j,k] = 0
                        self.t0[j, k] = self.t
        # </editor-fold>

        # <editor-fold desc="Compute initial q">
        self.q = np.zeros([ncol + 2, nrow + 2, 3])
        self.tclip = 0.
        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[j] + 1):
                if self.h[j, k] < 0.0:  # correct depth < 0
                    self.h[j, k] = 0.
                self.q[j, k, 0] = self.h[j, k]
                self.q[j, k, 1] = self.h[j, k] * self.u[j, k]
                self.q[j, k, 2] = self.h[j, k] * self.v[j, k]
                for i in range(inum[j, k]):
                    if (itype[j, k, i] == 5) & (fix[j, k, 0] != 0.):
                        # If solitary wave profile is subcritical and h > 0
                        if isurf == 1:
                            hdum = fix[j, k, 0] - self.zc[j, k]
                        else:
                            hdum = fix[j, k, 0]
                        if hdum <= 0.0:
                            print 'ERROR: depth specification for solitary wave \n BC is invalid'
                    if (itype[j, k, i] == 4):  # Specified flow rate (subcritical).
                        # For fixed flux BCs, must initialize depth.
                        print 'fixed flux BC is requested in cell ', j, k
                        jj, kk, j2, k2 = self.findbc(i, j, k)
                        if ipos[j, k, i] == 2 or ipos[j, k, i] == 4:  # If vertical faces
                            qflux = fix[j, k, 1] * cn[j, k, 0] + fix[j, k, 2] * sn[j, k, 0]  #
                            dx = deta[j, k, 1] * area[j, k]
                            dy = - deta[j, k, 0] * area[j, k]
                            dss = np.sqrt(dx * dx + dy * dy)
                            if self.h[j, k] < epsh:
                                print '  bed adjacent to the boundary is dry'
                                # if qflux_x *dz/dx < 0 and manning's N > 0
                                if (qflux * dz[j, k, 0] < 0.) & (xn > 0.):
                                    qflux = np.abs(qflux)  # qflux = abs(qflux)
                                    hnorm = (qflux * xn / np.sqrt(np.abs(dz[j, k, 0] / dss))) ** (3. / 5.)
                                    print '\t normal depth = ', hnorm, ' is specified'
                                    self.h[jj, kk] = hnorm
                                    self.hp[jj, kk] = hnorm
                                else:
                                    print '\t adverse slope or zero Manning n in cell'
                                    print '\t enter initial flow depth at specified'
                                    print '\t flux boundary ', j, k, ipos[j, k, i]
                                    self.h[jj, kk] = .1  # float(raw_input(""))
                                    self.hp[jj, kk] = self.h[jj, kk]
                            qflux = np.abs(qflux)
                            tclip1 = (xk * dss / (qflux * (1.0 - cappa))) ** (1.0 / (1.0 - ainflt))
                            if tclip1 > self.tclip:
                                self.tclip = tclip1

                        elif ipos[j, k, i] == 1 or ipos[j, k, i] == 3:  # If horizontal areas
                            qflux = fix[j, k, 1] * cn[j, k, 1] + fix[j, k, 2] * sn[j, k, 1]  #
                            dx = - dxi[j, k, 1] * area[j, k]
                            dy = dxi[j, k, 0] * area[j, k]
                            dss = np.sqrt(dx * dx + dy * dy)
                            if self.h[j, k] < epsh:
                                print 'bed adjacent to the boundary is dry'
                                if (qflux * dz[j, k, 1] < 0.) and (xn > 0.0):
                                    qflux = np.abs(qflux)
                                    hnorm = (qflux * xn / np.sqrt(np.abs(dz[j, k, 0] / dss))) ** (3. / 5.)
                                    print '\t normal depth =', hnorm, ' is specified'
                                    self.h[jj, kk] = hnorm
                                    self.hp[jj, kk] = hnorm
                                else:
                                    print '\t adverse slope or zero Manning n in cell'
                                    print '\t enter initial flow depth at specified'
                                    print '\t flux boundary =', j, k, ipos[j, k, i]
                                    self.h[jj, kk] = float(raw_input(""))
                                    self.hp[jj, kk] = self.h[jj, kk]  # set hp outside boundary equal to h inside
                                qflux = np.abs(qflux)
                                tclip2 = (xk * dss / (qflux * (1.0 - cappa))) ** (1.0 / (1.0 - ainflt))
                                if tclip2 > self.tclip:
                                    self.tclip = tclip2
                    if (itype[j, k, i] == 3) and (isurf == 1):
                        # For specified free surface BCs.
                        amp = 0
                        check = fix[j, k, 0] - self.zc[j, k] - amp
                        if check < 0.0:
                            print 'WARNING: specified tidal amplitude will'
                            print 'result in negative specified depth'
                            print 'specified depth will be set to zero'
        # </editor-fold>
        # <editor-fold desc="Compute zslope">

        if xk > 0.0:
            if self.tclip > 0.0:
                self.zslope = xk * self.tclip ** (ainflt - 1.0)
                print ' '
                print 'Kostiakov infiltration formula modified'
                print 'fot time less than ', self.tclip
                print 'infiltration slope is ', self.zslope
            else:
                self.tclip = 0.0
        # </editor-fold>

        self.itest = 0
        self.h, self.u, self.v = self.interp(self.h, self.u, self.v)
        for j in range(0, ncol + 1):
            for k in range(0, nrow + 1):
                if self.xc[j, k] == 0:
                    self.xc[j, k] = np.mean([self.xc[j - 1, k], self.xc[j + 1, k]])
                if self.yc[j, k] == 0:
                    self.yc[j, k] = np.mean([self.yc[j - 1, k], self.yc[j + 1, k]])
                if self.zc[j, k] == 0:
                    self.zc[j, k] = np.mean([self.zc[j - 1, k], self.zc[j + 1, k]])

        # <editor-fold desc="create the ouptut files">
        file = open("diag_py.out", "w")
        file.write('t\t vol \t xf  \n')
        file.close()
        file = open("dryout_py.dat", "w")
        file.close()
        # </editor-fold>

    def grid(self):
        """
        Called from the input function only
        Read grid data from file 'coords'.
        Creates global variables: x, y, z, nop
        nop:
            nop dimensions =  ncol+1,nrow+1,4
            each row of nop contains the node numbers of each cell face, counterclockwise
        """

        # <editor-fold desc="Read in coords">
        coords = []
        for line in open("coords", 'r'):
            item = line.rstrip()  # strip off newline and any other trailing whitespace
            coords.append(line)
        # </editor-fold>

        # <editor-fold desc="From 'coords': get npt, ne, x, y, z">
        global npt, ne, x, y, z
        npt, ne = [int(j) for j in coords[0].split()]
        x = np.zeros(npt + 1)
        y = np.zeros(npt + 1)
        z = np.zeros(npt + 1)
        # x, y, z are 1 indexed (the x[0] contains a dummy 0)
        # column vectors
        for i in range(1, npt + 1):
            x[i], y[i], z[i] = [float(l) for l in coords[i].split()]
        # </editor-fold>

        # <editor-fold desc="From 'coords: get nop">
        global nop
        nop = np.zeros([ncol + 1, nrow + 1, 4], dtype=int)
        count = 0
        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[j] + 1):
                nop[j, k, :] = [int(l) for l in coords[count + 1 + npt].split()]
                count += 1
        # </editor-fold>

        # <editor-fold desc="Compute grid metrics">
        global dxi, deta, dz, area, sx, sy  # grid transformation metric
        self.xc = np.zeros([ncol + 2, nrow + 2])  # coordinates of cell centers
        self.yc = np.zeros([ncol + 2, nrow + 2])
        self.zc = np.zeros([ncol + 2, nrow + 2])
        dxi = np.zeros([ncol + 2, nrow + 2, 2])
        deta = np.zeros([ncol + 2, nrow + 2, 2])
        dz = np.zeros([ncol + 2, nrow + 2, 2])
        area = np.zeros([ncol + 2, nrow + 2])
        sx = np.zeros([ncol + 2, nrow + 2])  # dz/dx
        sy = np.zeros([ncol + 2, nrow + 2])  # dz/sy
        # </editor-fold>
        for j in range(1, ncol + 1):  # Compute grid metrics
            for k in range(kbeg[j], kend[j] + 1):
                n1 = nop[j, k, 0]
                n2 = nop[j, k, 1]
                n3 = nop[j, k, 2]
                n4 = nop[j, k, 3]
                self.xc[j, k] = 0.25 * (x[n1] + x[n2] + x[n3] + x[n4])
                self.yc[j, k] = 0.25 * (y[n1] + y[n2] + y[n3] + y[n4])
                self.zc[j, k] = 0.25 * (z[n1] + z[n2] + z[n3] + z[n4])
                dxdxi = 0.5 * (-x[n1] + x[n2] + x[n3] - x[n4])  # 2
                dxdeta = 0.5 * (-x[n1] - x[n2] + x[n3] + x[n4])  # 0
                dydxi = 0.5 * (-y[n1] + y[n2] + y[n3] - y[n4])  # 0
                dydeta = 0.5 * (-y[n1] - y[n2] + y[n3] + y[n4])  #
                area[j, k] = dxdxi * dydeta - dxdeta * dydxi
                if area[j, k] < 0:
                    print "error"
                dxi[j, k, 0] = dydeta / area[j, k]
                deta[j, k, 0] = -dydxi / area[j, k]
                dxi[j, k, 1] = -dxdeta / area[j, k]
                deta[j, k, 1] = dxdxi / area[j, k]
                sx[j, k] = ((z[n2] - z[n4]) * (y[n3] - y[n1]) - (z[n3] - z[n1]) * (y[n2] - y[n4])) / (2. * area[j, k])
                sy[j, k] = ((z[n3] - z[n1]) * (x[n2] - x[n4]) - (z[n2] - z[n4]) * (x[n3] - x[n1])) / (2. * area[j, k])
                dz[j, k, 0] = sx[j, k] * dxdxi + sy[j, k] * dydxi
                dz[j, k, 1] = sx[j, k] * dxdeta + sy[j, k] * dydeta

        # <editor-fold desc="Compute cell face angles">
        global ds, sn, cn
        ds = np.zeros([ncol + 2, nrow + 2, 2])
        sn = np.zeros([ncol + 2, nrow + 2, 2])
        cn = np.zeros([ncol + 2, nrow + 2, 2])
        # </editor-fold>
        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[j] + 1):
                # Horizontal face.
                ddx = x[nop[j, k, 1]] - x[nop[j, k, 0]]  # x[n2] - x[n1]  = 1
                ddy = y[nop[j, k, 1]] - y[nop[j, k, 0]]  # y[n2] - y[n1]  = 0
                ds[j, k, 1] = np.sqrt(ddx * ddx + ddy * ddy)  # = 1 = area
                sn[j, k, 1] = ddx / ds[j, k, 1]  # 1
                cn[j, k, 1] = -ddy / ds[j, k, 1]  # 0
                # Vertical face
                ddx = x[nop[j, k, 3]] - x[nop[j, k, 0]]  # x[n4] - x[n1] = 0
                ddy = y[nop[j, k, 3]] - y[nop[j, k, 0]]  # y[n4] - y[n1]  = 1
                ds[j, k, 0] = np.sqrt(ddx * ddx + ddy * ddy)  # = 1
                sn[j, k, 0] = -ddx / ds[j, k, 0]  # 0
                cn[j, k, 0] = ddy / ds[j, k, 0]  # 1
                for i in range(inum[j, k]):
                    if ipos[j, k, i] == 3:  # Top (boundary) faces.
                        ddx = x[nop[j, k, 2]] - x[nop[j, k, 3]]
                        ddy = y[nop[j, k, 2]] - y[nop[j, k, 3]]
                        ds[j, k + 1, 1] = np.sqrt(ddx * ddx + ddy * ddy)
                        sn[j, k + 1, 1] = ddx / ds[j, k + 1, 1]
                        cn[j, k + 1, 1] = -ddy / ds[j, k + 1, 1]
                    elif ipos[j, k, i] == 2:  # Right (boundary) faces.
                        ddx = x[nop[j, k, 2]] - x[nop[j, k, 1]]
                        ddy = y[nop[j, k, 2]] - y[nop[j, k, 1]]
                        ds[j + 1, k, 0] = np.sqrt(ddx * ddx + ddy * ddy)
                        sn[j + 1, k, 0] = -ddx / ds[j + 1, k, 0]
                        cn[j + 1, k, 0] = ddy / ds[j + 1, k, 0]

        for j in range(1, ncol + 1):  # Compute some things in ghost cells
            for k in range(kbeg[j], kend[j] + 1):
                for i in range(inum[j, k]):  # loop over cell faces
                    #  jj, kk is the cell immediately across the boundary face i
                    #  j2, k2 is the cell on the opposing size of j,k from jj,kk
                    jj, kk, j2, k2 = self.findbc(i, j, k)
                    if (jj < 0) or (kk < 0):
                        print jj, kk
                    area[jj, kk] = area[j, k]
                    sx[jj, kk] = sx[j, k]
                    sy[jj, kk] = sy[j, k]
                    dxi[jj, kk, 0] = dxi[j, k, 0]
                    deta[jj, kk, 0] = deta[j, k, 0]
                    dxi[jj, kk, 1] = dxi[j, k, 1]
                    deta[jj, kk, 1] = deta[j, k, 1]
                    self.xc[jj, kk] = 2. * self.xc[j, k] - self.xc[j2, k2]
                    self.yc[jj, kk] = 2. * self.yc[j, k] - self.yc[j2, k2]
                    self.zc[jj, kk] = 2. * self.zc[j, k] - self.zc[j2, k2]

    def findbc(self, i, j, k):
        """
        findbc gets called from the grid function:
        :param i:  dummy index for the face number in ipos
        :param j:  dummy index for column number (eventually row)
        :param k:  dummy index for row number (eventually column)
        :return:  (jj,kk): the indices of the ghost cell
                   (j2, k2):  indices of the cell opposite the ghost cell
        """
        if ipos[j, k, i] == 1.:
            jj = j
            kk = k - 1
            j2 = j
            k2 = k + 1
        elif ipos[j, k, i] == 2.:
            jj = j + 1
            kk = k
            j2 = j - 1
            k2 = k
        elif ipos[j, k, i] == 3.:
            jj = j
            kk = k + 1
            j2 = j
            k2 = k - 1
        elif ipos[j, k, i] == 4.:
            jj = j - 1
            kk = k
            j2 = j + 1
            k2 = k
        else:
            print "error -  this is not a boundary cell!"
            jj = j
            kk = k
            j2 = j
            k2 = k
        return jj, kk, j2, k2

    def interp(self, h, u, v):
        """
        corrects  xc, yc, h, u, v in the corner cells
        :param h: h or hp at current time step
        :param u: u or up
        :param v: v or vp
        :return:  h, u, v
        """
        for j in range(1, ncol + 1):
            for k in range(kbeg[j], kend[j] + 1):
                # Extrapolate to corner cells
                ibc1 = 0
                ibc2 = 0
                ibc3 = 0
                ibc4 = 0
                for i in range(inum[j, k]):
                    if ipos[j, k, i] == 1:
                        ibc1 = 1
                    if ipos[j, k, i] == 2:
                        ibc2 = 1
                    if ipos[j, k, i] == 3:
                        ibc3 = 1
                    if ipos[j, k, i] == 4:
                        ibc4 = 1
                if (ibc1 == 1) and (ibc4 == 1):
                    self.xc[j - 1, k - 1] = 2. * self.xc[j, k] - self.xc[j + 1, k + 1]
                    self.yc[j - 1, k - 1] = 2. * self.yc[j, k] - self.yc[j + 1, k + 1]
                    self.zc[j - 1, k - 1] = 2. * self.zc[j, k] - self.zc[j + 1, k + 1]
                    h[j - 1, k - 1] = 2. * h[j, k] - h[j + 1, k + 1]
                    u[j - 1, k - 1] = 2. * u[j, k] - u[j + 1, k + 1]
                    v[j - 1, k - 1] = 2. * v[j, k] - v[j + 1, k + 1]
                elif (ibc1 == 1) and (ibc2 == 1):
                    self.xc[j + 1, k - 1] = 2. * self.xc[j, k] - self.xc[j - 1, k + 1]
                    self.yc[j + 1, k - 1] = 2. * self.yc[j, k] - self.yc[j - 1, k + 1]
                    self.zc[j + 1, k - 1] = 2. * self.zc[j, k] - self.zc[j - 1, k + 1]
                    h[j + 1, k - 1] = 2. * h[j, k] - h[j - 1, k + 1]
                    u[j + 1, k - 1] = 2. * u[j, k] - u[j - 1, k + 1]
                    v[j + 1, k - 1] = 2. * v[j, k] - v[j - 1, k + 1]
                elif (ibc2 == 1) and (ibc3 == 1):
                    self.xc[j + 1, k + 1] = 2. * self.xc[j, k] - self.xc[j - 1, k - 1]
                    self.yc[j + 1, k + 1] = 2. * self.yc[j, k] - self.yc[j - 1, k - 1]
                    self.zc[j + 1, k + 1] = 2. * self.zc[j, k] - self.zc[j - 1, k - 1]
                    h[j + 1, k + 1] = 2. * h[j, k] - h[j - 1, k - 1]
                    u[j + 1, k + 1] = 2. * u[j, k] - u[j - 1, k - 1]
                    v[j + 1, k + 1] = 2. * v[j, k] - v[j - 1, k - 1]
                elif (ibc3 == 1) and (ibc4 == 1):
                    self.xc[j - 1, k + 1] = 2. * self.xc[j, k] - self.xc[j + 1, k - 1]
                    self.yc[j - 1, k + 1] = 2. * self.yc[j, k] - self.yc[j + 1, k - 1]
                    self.zc[j - 1, k + 1] = 2. * self.zc[j, k] - self.zc[j + 1, k - 1]
                    h[j - 1, k + 1] = 2. * h[j, k] - h[j + 1, k - 1]
                    u[j - 1, k + 1] = 2. * u[j, k] - u[j + 1, k - 1]
                    v[j - 1, k + 1] = 2. * v[j, k] - v[j + 1, k - 1]
                    h, u, v = self.bconds(j, k, h, u, v)
        return h, u, v

    def bconds(self, j, k, hdum, udum, vdum):
        """
        :param j: dummmy column variable.  because this is python, we'll have to reverse this
        :param k: dummy row variable
        :param hdum:  h (ncol+2, nrow+2)
        :param udum:  u (ncol+2, nrow+2)
        :param vdum:  v (ncol+2, nrow+2)
        :return:  returns to hdum, udum, vdum
            hdum, udum, vdum may refer to h,u,v (in the input and corrector steps),
            or hp, up, vp (in the predictor and flux calculation steps)
        :update:  updates variables self.dh, self.du, self.dv
        """
        for i in range(inum[j, k]):
            if ipos[j, k, i] == 1:  # front face
                jj = j
                kk = k - 1
                jl = j
                kl = k
                j2 = j
                k2 = k + 1
                io = 1
            elif ipos[j, k, i] == 2:  # right face
                jj = j + 1
                kk = k
                jl = j + 1
                kl = k
                j2 = j - 1
                k2 = k
                io = 0
            elif ipos[j, k, i] == 3:  # back face
                jj = j
                kk = k + 1
                jl = j
                kl = k + 1
                j2 = j
                k2 = k - 1
                io = 1
            elif ipos[j, k, i] == 4:  # left face
                jj = j - 1
                kk = k
                jl = j
                kl = k
                j2 = j + 1
                k2 = k
                io = 0
            self.t0[jj, kk] = self.t0[j, k]
            if itype[j, k, i] == 0:  # This is an open boundary
                self.dh[jj, kk, io] = self.dh[j, k, io]
                self.du[jj, kk, io] = self.du[j, k, io]
                self.dv[jj, kk, io] = self.dv[j, k, io]
                hdum[jj, kk] = 2.0 * hdum[j, k] - hdum[j2, k2]
                udum[jj, kk] = 2.0 * udum[j, k] - udum[j2, k2]
                vdum[jj, kk] = 2.0 * vdum[j, k] - vdum[j2, k2]

            elif itype[j, k, i] == 1:  # Wall boundary
                self.dh[jj, kk, io] = self.dh[j, k, io]
                self.du[jj, kk, io] = self.du[j, k, io]
                self.dv[jj, kk, io] = self.dv[j, k, io]
                hdum[jj, kk] = 2.0 * hdum[j, k] - hdum[j2, k2]
                udum[jj, kk] = udum[j, k] * (sn[jl, kl, io] * sn[jl, kl, io] -
                                             cn[jl, kl, io] * cn[jl, kl, io]) - \
                               2.0 * vdum[j, k] * sn[jl, kl, io] * cn[jl, kl, io]
                vdum[jj, kk] = vdum[j, k] * (cn[jl, kl, io] * cn[jl, kl, io] -
                                             sn[jl, kl, io] * sn[jl, kl, io]) - \
                               2.0 * udum[j, k] * sn[jl, kl, io] * cn[jl, kl, io]

            elif itype[j, k, i] == 2:  # Specified depth and velocity (supercritical)
                self.dh[jj, kk, io] = 0.0
                self.du[jj, kk, io] = 0.0
                self.dv[jj, kk, io] = 0.0
                hdum[jj, kk] = fix[j, k, 0]
                udum[jj, kk] = fix[j, k, 1]
                vdum[jj, kk] = fix[j, k, 2]
                if isurf == 1:
                    hdum[jj, kk] = hdum[jj, kk] - self.zc[j, k]

            elif itype[j, k, i] == 3:  # Tidal (subcritical)
                self.dh[jj, kk, io] = 0.0
                self.du[jj, kk, io] = 0.0
                self.dv[jj, kk, io] = 0.0
                hdum[jj, kk] = fix[j, k, 0] + amp * np.cos(-6.283185 * self.t / period)
                udum[jj, kk] = udum[j, k]
                vdum[jj, kk] = vdum[j, k]
                if isurf == 1:
                    hdum[jj, kk] = hdum[jj, kk] - self.zc[j, k]

            elif itype[j, k, i] == 4:  # Specified flow rate (subcritical)
                self.du[jj, kk, io] = 0.0
                self.dv[jj, kk, io] = 0.0
                if hdum[j, k] > epsh:
                    hdum[jj, kk] = hdum[j, k]
                    self.dh[jj, kk, io] = 0.0
                elif hdum[j, k] * hdum[j2, k2] > epsh:
                    hdum[jj, kk] = 2.0 * hdum[j, k] - hdum[j2, k2]
                    self.dh[jj, kk, io] = self.dh[j, k, io]

                if hdum[jj, kk] > epsh:
                    udum[jj, kk] = fix[j, k, 1] / hdum[jj, kk]
                    vdum[jj, kk] = fix[j, k, 2] / hdum[jj, kk]
                else:
                    udum[jj, kk] = 0.0
                    vdum[jj, kk] = 0.0

            elif itype[j, k, i] == 5:  # Solitary wave profile (subcritical).
                self.dh[jj, kk, io] = 0.0
                self.du[jj, kk, io] = 0.0
                self.dv[jj, kk, io] = 0.0
                dum = period - self.t
                if dum < -36.0:
                    dum = -36.0
                hdum[jj, kk] = fix[j, k, 1] + 2.0 * amp / (np.cosh(dum) * np.cosh(dum))
                udum[jj, kk] = udum[j, k]
                vdum[jj, kk] = vdum[j, k]
                if isurf == 1:
                    hdum[jj, kk] = hdum[jj, kk] - self.zc[j, k]
                if hdum[jj, kk] < 0.0:
                    hdum[jj, kk] = 0.0
                if hdum[jj, kk] < epsh:
                    udum[jj, kk] = 0.0
                    vdum[jj, kk] = 0.0
            #if hdum[jj, kk] > 10:
            #    print 'error'
        return hdum, udum, vdum

    def transform(self):
        global j1, k1, w1, w2, w3, w4
        file = open("weights.out", "w")
        j1 = np.zeros(npt + 1, dtype=int)
        k1 = np.zeros(npt + 1, dtype=int)
        w1 = np.zeros(npt + 1, dtype=float)
        w2 = np.zeros(npt + 1, dtype=float)
        w3 = np.zeros(npt + 1, dtype=float)
        w4 = np.zeros(npt + 1, dtype=float)
        # bracket sampling location
        for i in range(1, npt + 1):
            xdum = x[i]
            ydum = y[i]
            for j in range(1, ncol + 1):
                for k in range(kbeg[j], kend[j] + 2):
                    if (xdum >= self.xc[j - 1, k]) and (xdum < self.xc[j, k]):
                        if (ydum >= self.yc[j, k - 1]) and (ydum < self.yc[j, k]):
                            j1[i] = j
                            k1[i] = k
                            # print xdum, j, xc[j-1,k]
                            break
                    elif (xdum >= self.xc[j, k]) and (xdum <= self.xc[j + 1, k]):
                        if (ydum >= self.yc[j, k - 1]) and (ydum <= self.yc[j, k]):
                            j1[i] = j + 1
                            k1[i] = k
                            break
            # Compute interpolation distances.
            d1 = np.sqrt((xdum - self.xc[j1[i], k1[i]]) ** 2 + (ydum - self.yc[j1[i], k1[i]]) ** 2)
            d2 = np.sqrt((xdum - self.xc[j1[i] - 1, k1[i]]) ** 2 + (ydum - self.yc[j1[i] - 1, k1[i]]) ** 2)
            d3 = np.sqrt((xdum - self.xc[j1[i] - 1, k1[i] - 1]) ** 2 + (ydum - self.yc[j1[i] - 1, k1[i] - 1]) ** 2)
            d4 = np.sqrt((xdum - self.xc[j1[i], k1[i] - 1]) ** 2 + (ydum - self.yc[j1[i], k1[i] - 1]) ** 2)
            # Compute weights as distance inverses.
            # Tecplot states that the 3.5 exponent yields the smoothest results.
            if (d1 > 10) or (d2 > 10) or (d3 > 10) or (d4 > 10):
                print d1, d2, d3, d4
            if d1 == 0.:
                w1[i] = 1.
                w2[i] = 0.
                w3[i] = 0.
                w4[i] = 0.
            elif d2 == 0.:
                w1[i] = 0.
                w2[i] = 1.
                w3[i] = 0.
                w4[i] = 0.
            elif d3 == 0.:
                w1[i] = 0.
                w2[i] = 0.
                w3[i] = 1.
                w4[i] = 0.
            elif d4 == 0.:
                w1[i] = 0.
                w2[i] = 0.
                w3[i] = 0.
                w4[i] = 1.
            else:
                w1[i] = d1 ** (-3.5)
                w2[i] = d2 ** (-3.5)
                w3[i] = d3 ** (-3.5)
                w4[i] = d4 ** (-3.5)
            file.write('{0}, {1}, {2}, {3}, {4}, {5}, {6}\n'.format(i, j1[i], k1[i], w1[i], w2[i], w3[i], w4[i]))
        file.close()

    def get_weights(self):
        global j1, k1, w1, w2, w3, w4
        weights = [[0, 0, 0, 0, 0, 0, 0]]
        for line in open("weights.out", 'r'):
            weights.append([float(x) for x in line.strip().split(", ")])
        weights = np.asarray(weights)
        j1 = weights[:, 1].astype(int)
        k1 = weights[:, 2].astype(int)
        w1 = weights[:, 3].astype(float)
        w2 = weights[:, 4].astype(float)
        w3 = weights[:, 5].astype(float)
        w4 = weights[:, 6].astype(float)

    def movie(self):
        self.hn = np.zeros(npt + 1)
        self.un = np.zeros(npt + 1)
        self.vn = np.zeros(npt + 1)

        for i in range(1, npt + 1):
            # Interpolate data
            sumw = w1[i] + w2[i] + w3[i] + w4[i]
            self.hn[i] = (w1[i] * self.h[j1[i], k1[i]] + w2[i] * self.h[j1[i] - 1, k1[i]] +
                          w3[i] * self.h[j1[i] - 1, k1[i] - 1] + w4[i] * self.h[j1[i], k1[i] - 1]) / sumw
            self.un[i] = (w1[i] * self.u[j1[i], k1[i]] + w2[i] * self.u[j1[i] - 1, k1[i]] +
                          w3[i] * self.u[j1[i] - 1, k1[i] - 1] + w4[i] * self.u[j1[i], k1[i] - 1]) / sumw
            self.vn[i] = (w1[i] * self.v[j1[i], k1[i]] + w2[i] * self.v[j1[i] - 1, k1[i]] +
                          w3[i] * self.v[j1[i] - 1, k1[i] - 1] + w4[i] * self.v[j1[i], k1[i] - 1]) / sumw

            # z[i] = (w1[i]*zc[j1[i],k1[i]] + w2[i]*zc[j1[i]-1,k1[i]] +
            #         w3[i]*zc[j1[i]-1,k1[i]-1] + w4[i]*zc[j1[i],k1[i]-1])/sumw

        if self.t == 0.:
            file = open("movie_py.dat", "w")
            file.write('VARIABLES = "X", "Y", "Z", "H", "D"\n')
            file.write('ZONE N=%i, E=%i, F=FEPOINT,  ET=QUADRILATERAL \n' % (npt, ne))
            for n in range(1, npt + 1):
                file.write('%e %e %e %e %e \n' % (x[n], y[n], z[n], z[n], z[n]))
            for j in range(1, ncol + 1):
                for k in range(kbeg[j], kend[j] + 1):
                    file.write('  '.join(['%i' % nop[j, k, i] for i in range(4)]))
                    file.write('\n')
            file.close()
        elif np.mod(self.it, 30) == 0:
            file = open("movie_py.dat", "a")
            file.write('VARIABLES = "X", "Y", "Z", "H", "D"\n')
            file.write('ZONE N=   %i E=   %i F=FEPOINT,  ET=QUADRILATERAL  D=(FECONNECT,1,2,3)\n' % (npt, ne))
            for n in range(1, npt + 1):
                file.write(''.join(['%f' % (self.hn[n] + z[n]), ', %f' % self.hn[n]]))
                file.write('\n')
            file.close()

    def predict(self, j, k):
        """
        :param j: column index
        :param k: row
        :return: computes gradients in wet cells :
            calls self.limitr --> modifes self.dh, self.du, self.dv
            Sets gradients in dry cells to zero.
        """
        if ilim > 0:
            # loop over coordinate directions.
            for kk in range(2):
                if kk == 1:
                    jr = j + 1
                    jl = j - 1
                    kr = k
                    kl = k
                else:
                    jr = j
                    jl = j
                    kr = k + 1
                    kl = k - 1
                if self.h[j, k] >= epsh:
                    # Compute gradients only in wet cells.
                    # Limit free surface elevation to reduce dissipation..
                    dh1 = self.h[j, k] + self.zc[j, k] - self.h[jl, kl] - self.zc[jl, kl]
                    dh2 = self.h[jr, kr] + self.zc[jr, kr] - self.h[j, k] - self.zc[j, k]
                    if self.h[jl, kl] < epsh:
                        dh1 = 2. * dh1
                    if self.h[jr, kr] < epsh:
                        dh2 = 2. * dh2
                    dhh = self.limitr(dh1, dh2)
                    self.dh[j, k, kk] = dhh - dz[j, k, kk]

                    # U velocity.
                    du1 = self.u[j, k] - self.u[jl, kl]
                    du2 = self.u[jr, kr] - self.u[j, k]
                    self.du[j, k, kk] = self.limitr(du1, du2)

                    # V velocity.
                    dv1 = self.v[j, k] - self.v[jl, kl]
                    dv2 = self.v[jr, kr] - self.v[j, k]
                    self.dv[j, k, kk] = self.limitr(dv1, dv2)
                else:
                    self.dh[j, k, kk] = 0.
                    self.du[j, k, kk] = 0.
                    self.dv[j, k, kk] = 0.

            # Generalized velocities.
            uxi = self.u[j, k] * dxi[j, k, 0] + self.v[j, k] * dxi[j, k, 1]  # if no transform: self.u
            ueta = self.u[j, k] * deta[j, k, 0] + self.v[j, k] * deta[j, k, 1]  # if no transform: self.v
            qs = self.source(j, k, self.h[j, k], self.u[j, k], self.v[j, k])
            if self.h[j, k] >= epsh ** 0.75:
                qs[1] = qs[1] / self.h[j, k]
                qs[2] = qs[2] / self.h[j, k]
            else:
                qs[0] = 0.
                qs[1] = 0.
                qs[2] = 0.
            self.hp[j, k] = self.h[j, k] - 0.50 * self.dt * (
                uxi * self.dh[j, k, 0] + self.h[j, k] * (dxi[j, k, 0] * self.du[j, k, 0] +
                                                         dxi[j, k, 1] * self.dv[j, k, 0]) +
                ueta * self.dh[j, k, 1] + self.h[j, k] * (deta[j, k, 0] * self.du[j, k, 1] +
                                                          deta[j, k, 1] * self.dv[j, k, 1]) + qs[0])

            self.up[j, k] = self.u[j, k] - 0.5 * self.dt * (
                grav * dxi[j, k, 0] * self.dh[j, k, 0] + uxi * self.du[j, k, 0] +
                grav * deta[j, k, 0] * self.dh[j, k, 1] + ueta * self.du[j, k, 1] + qs[1])

            self.vp[j, k] = self.v[j, k] - 0.5 * self.dt * (
                grav * dxi[j, k, 1] * self.dh[j, k, 0] + uxi * self.dv[j, k, 0] +
                grav * deta[j, k, 1] * self.dh[j, k, 1] + ueta * self.dv[j, k, 1] + qs[2])

        else:
            self.hp[j, k] = self.h[j, k]
            self.up[j, k] = self.u[j, k]
            self.vp[j, k] = self.v[j, k]
        if self.hp[j, k] < 0.0:  # Correct any negative depths.
            self.hp[j, k] = 0.0
            self.dh[j, k, 0] = 0.0
            self.dh[j, k, 1] = 0.0
        if self.hp[j, k] <= 0.:  # Neglect momentum in nearly dry cells.
            self.up[j, k] = 0.
            self.vp[j, k] = 0.
            for i in range(2):
                self.du[j, k, i] = 0.
                self.dv[j, k, i] = 0.

    def limitr(self, dq1, dq2):
        if (ilim == 1):  # Lax-Wendroff.
            dq = dq2
        elif ilim == 2.:  # Beam-Warming.
            dq = dq1
        elif (ilim == 3):  # Fromm
            dq = 0.5 * (dq1 + dq2)
        elif ilim == 4:  # Double minmod.
            a = 0.5 * (dq1 + dq2)
            b = 2.0 * dq1
            c = 2.0 * dq2
            if (a * b > 0.0) and (b * c > 0.0):
                dq = fmin1(a, b, c)
            else:
                dq = 0
        else:
            if dq1 * dq2 < 0.0:
                dq = 0.0
            else:
                dq = fmin2(fmax2(dq1, dq2), beta * fmin2(dq1, dq2))
        return dq

    ## global variables here are: t, t0, tc, ainflt, binflt, xk
    def source(self, j, k, hdum, udum, vdum):
        qs = np.zeros(3)
        if (hdum >= epsh):
            tnew = max(0.0, self.t - self.t0[j, k])
            if (tnew <= 0.0):
                znew = 0.
            elif (tnew > 0.) and (tnew <= self.tclip):
                znew = self.zslope * tnew
            elif (tnew > self.tclip) and (tnew <= tc):
                znew = xk * tnew ** ainflt
            elif tnew > tc:
                znew = xk * tc ** ainflt + binflt * (tnew - tc)
            told = max(0.0, self.t - self.t0[j, k] - self.dt)
            if told <= 0.0:
                zold = 0.0
            elif (told > 0.0) and (told <= self.tclip):
                zold = self.zslope * told
            elif (told > self.tclip) and (told <= tc):
                zold = xk * told ** ainflt
            elif told > tc:
                zold = xk * tc ** ainflt + binflt * (told - tc)
            self.winflt = (zold - znew) / self.dt
            vmag = np.sqrt(udum * udum + vdum * vdum)
            fricx = grav * xn * xn * udum * vmag / hdum ** (1. / 3.)
            fricy = grav * xn * xn * vdum * vmag / hdum ** (1. / 3.)
            qs[0] = self.winflt
            qs[1] = 0.5 * udum * self.winflt - fricx - grav * hdum * sx[j, k]
            qs[2] = 0.5 * vdum * self.winflt - fricy - grav * hdum * sy[j, k]
        else:
            qs[0] = 0.
            qs[1] = 0.
            qs[2] = 0.
            self.winflt = 0
        return qs

    def fluxes(self, jl, jr, kl, kr, i1):
        # MUSCL extrapolation at cell interface.

        hl = self.hp[jl, kl] + 0.5 * self.dh[jl, kl, i1]
        ul = self.up[jl, kl] + 0.5 * self.du[jl, kl, i1]
        vl = self.vp[jl, kl] + 0.5 * self.dv[jl, kl, i1]
        hr = self.hp[jr, kr] - 0.5 * self.dh[jr, kr, i1]
        ur = self.up[jr, kr] - 0.5 * self.du[jr, kr, i1]
        vr = self.vp[jr, kr] - 0.5 * self.dv[jr, kr, i1]
        snn = sn[jr, kr, i1]
        cnn = cn[jr, kr, i1]
        if i1 == 0:
            dx = deta[jr, kr, 1] * area[jr, kr]
            dy = -deta[jr, kr, 0] * area[jr, kr]
        else:
            dx = -dxi[jr, kr, 1] * area[jr, kr]
            dy = dxi[jr, kr, 0] * area[jr, kr]
        # Needed for dry bed problems.
        if (hl < 0.):
            hl = 0.
        if (hr < 0.):
            hr = 0.
        # Compute arithmetic averages for source terms.
        havg = 0.5 * (hl + hr)
        uavg = 0.5 * (ul + ur)
        vavg = 0.5 * (vl + vr)

        iflag = 0.
        if (ilim == 0):
            # Prevent leakage into cells with higher bed elevation.
            etal = self.hp[jl, kl] + self.zc[jl, kl]
            etar = self.hp[jr, kr] + self.zc[jr, kr]
            sdumx = 0.5 * (sx[jl, kl] + sx[jr, kr])
            sdumy = 0.5 * (sy[jl, kl] + sy[jr, kr])
            if (self.hp[jl, kl] >= epsh) and (self.hp[jr, kr] < epsh):
                if etal <= self.zc[jr, kr]:
                    iflag = 1.
            if (self.hp[jr, kr] >= epsh) and (self.hp[jl, kl] < epsh):
                if etar <= self.zc[jl, kl]:
                    iflag = 1.
            # C Bed infiltration, bottom drag.
            if havg < epsh:
                self.qsource[jr, kr, 0, i1] = 0.
                self.qsource[jr, kr, 1, i1] = -0.5 * (grav * havg * sdumx) * dx * cnn
                self.qsource[jr, kr, 2, i1] = -0.5 * (grav * havg * sdumy) * dy * snn
            else:
                tavg = 0.5 * (self.t0[jl, kl] + self.t0[jr, kr])
                topp = self.t - tavg
                if topp > 0.:
                    if topp <= tc:
                        zinflt = xk * topp ** ainflt
                    else:
                        zinflt = binflt
                else:
                    zinflt = 0.0

                vmag = np.sqrt(uavg * uavg + vavg * vavg)
                fricx = grav * xn * xn * uavg * vmag / havg **(1./3.)
                fricy = grav * xn * xn * vavg * vmag / havg **(1./3.)

                self.qsource[jr, kr, 0, i1] = -0.25 * self.winflt * (cnn * dx + snn * dy)
                self.qsource[jr, kr, 1, i1] = 0.5 * (0.5 * uavg * self.winflt - fricx -
                                                     grav * havg * sdumx) * dx * cnn
                self.qsource[jr, kr, 2, i1] = 0.5 * (0.5 * vavg * self.winflt - fricy -
                                                     grav * havg * sdumy) * dy * snn
        if havg <= 0 or iflag == 1:
            for i in range(3):
                self.f[jr, kr, i, i1] = 0.
                self.qsource[jr, kr, i, i1] = 0.
        else:
            sdumx = 0
            sdumy = 0
            fdum = self.solver(hl, hr, ul, ur, vl, vr, snn, cnn, dx, dy, sdumx, sdumy)
            for i in range(3):
                self.f[jr, kr, i, i1] = fdum[i]

    def solver(self, hl, hr, ul, ur, vl, vr, sn, cn, dx, dy, sx, sy):

        # <editor-fold desc="Compute Roe averages at cell face: duml, dumr, hhat, uhat, vhat, chat">
        duml = np.sqrt(hl)
        dumr = np.sqrt(hr)
        hhat = duml * dumr
        uhat = (duml * ul + dumr * ur) / (duml + dumr)
        vhat = (duml * vl + dumr * vr) / (duml + dumr)
        chat = np.sqrt(0.5 * grav * (hl + hr))
        uperp = uhat * cn + vhat * sn
        # </editor-fold>
        # <editor-fold desc="Compute eigenvalues (a) and approximate wave strengths (dh, du, dv)">
        a = np.array([uperp - chat, uperp, uperp + chat])
        # Compute approximate wave strengths.
        dh = hr - hl
        du = ur - ul
        dv = vr - vl

        dupar = -du * sn + dv * cn
        duperp = du * cn + dv * sn
        ws = np.array([0.5 * (dh - hhat * duperp / chat),
                       hhat * dupar,
                       0.5 * (dh + hhat * duperp / chat)])
        # </editor-fold>
        # <editor-fold desc="Compute right eigenvectors: (e)">
        e = np.zeros([3, 3])
        e[0, 0] = 1.
        e[1, 0] = uhat - chat * cn
        e[2, 0] = vhat - chat * sn
        e[0, 1] = 0.
        e[1, 1] = -sn
        e[2, 1] = cn
        e[0, 2] = 1.
        e[1, 2] = uhat + chat * cn
        e[2, 2] = vhat + chat * sn
        # </editor-fold>

        # Compute approximate source terms for upwind method only.
        if sw == 1:
            zdum = np.zeros(3)
            if hhat > epsh:
                vmag = np.sqrt(uhat * uhat + vhat * vhat)
                fricx = grav * xn * xn * uhat * vmag / hhat ** (1. / 3.)
                fricy = grav * xn * xn * vhat * vmag / hhat ** (1. / 3.)
            else:
                fricx = 0.
                fricy = 0.
            q1 = - self.winflt
            q2 = - 0.5 * cn * (chat * chat * sx + fricx - 0.5 * uhat * self.winflt) * dx * cn
            q3 = - 0.5 * sn * (chat * chat * sy + fricy - 0.5 * vhat * self.winflt) * dy * sn
            zdum[0] = np.abs(a[0]) / a[0] * (a[2] * q1 - cn * q2 - sn * q3) / (2. * chat)
            if (np.abs(sn * a[1]) > 1e-8):
                zdum[1] = np.abs(a[1]) / a[1] * ((uhat - uperp * cn) * q1 +
                                                 sn * sn * q2 - sn * cn * q3) / sn
            else:
                zdum[1] = 0.
            zdum[2] = np.abs(a[2]) / a[2] * (-a[0] * q1 + cn * q2 + sn * q3) / (2. * chat)

        # <editor-fold desc="Entropy fix">
        dl = np.sqrt(dx * dx + dy * dy)
        cl = np.sqrt(grav * hl)
        cr = np.sqrt(grav * hr)
        uperpl = ul * cn + vl * sn
        uperpr = ur * cn + vr * sn
        al1 = uperpl - cl
        al3 = uperpl + cl
        ar1 = uperpr - cr
        ar3 = uperpr + cr
        da = np.zeros(3)
        da[0] = max(0., 4. * (ar1 - al1))
        da[1] = 0.
        da[2] = max(0., 4. * (ar3 - al3))
        astar = np.zeros(3)
        for i in range(3):
            if (np.abs(a[i]) < 0.5 * da[i]):
                astar[i] = a[i] * a[i] / da[i] + 0.25 * da[i]
            else:
                astar[i] = np.abs(a[i])
            if (astar[i] / dl > self.amax):
                self.amax = astar[i] / dl
        # </editor-fold>
        # <editor-fold desc="Compute flux increments">
        dum = np.zeros(3)
        for i in range(3):
            dum[i] = 0
            for l in range(3):
                dum[i] += (astar[l] * ws[l] - sw * z[l]) * e[i, l]
        # </editor-fold>

        f = np.zeros(3)
        f[0] = 0.5 * (f1(hl, uperpl) + f1(hr, uperpr) - dum[0])
        f[1] = 0.5 * (f2(hl, ul, uperpl, cn)
                      + f2(hr, ur, uperpr, cn) - dum[1])
        f[2] = 0.5 * (f3(hl, vl, uperpl, sn)
                      + f3(hr, vr, uperpr, sn) - dum[2])
        return f

    def plot(self):
        #plt.close()
        plt.pcolormesh(self.h)
        #plt.colorbar()


def f1(h, up):
    f1 = h * up
    return f1


def f2(h, u, up, cn):
    f2 = h * u * up + 0.5 * grav * h * h * cn
    return f2


def f3(h, v, up, sn):
    f3 = h * v * up + 0.5 * grav * h * h *   sn
    return f3


def fmin2(a, b):
    if a < 0.:
        return - min(np.abs(a), np.abs(b))
    else:
        return min(a, b)


def fmax2(a, b):
    if a < 0.:
        return - max(np.abs(a), np.abs(b))
    else:
        return max(a, b)


def fmin1(a, b, c):
    if (a < 0.0):
        return -1. * min(np.abs(a), np.abs(b), np.abs(c))
    else:
        return min(a, b, c)


d = dry()
# while d.t < tmax:
#     d.step()
#     d.plot()
print "done"
