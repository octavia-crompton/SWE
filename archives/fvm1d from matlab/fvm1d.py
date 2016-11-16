import numpy as np
import matplotlib.pyplot as plt


class fvm1d:
    def __init__(self, nc=100, L=1000., dt=0.5, nt=100., ntplot=10, beta=2.):
        L = np.float128(L)
        dx = np.float128(L / nc)
        x = np.arange(0, L + dx, dx)
        xc = np.arange(dx/2, L+dx-dx/2, dx)
        self.beta = beta
        self.nc = nc  # number of cells
        self.nf = nc + 1  # number of edges
        self.L = L  # length of channel
        self.dx = dx  # length of cell
        self.x = x  # array of edge coordinates
        self.xc = xc

        # Set up time marching and output interval
        self.dt = dt  # time step (s)
        self.nt = nt  # number of time steps
        self.ntplot = ntplot  # plot interval (number of time steps)

        # Define bed elevation at faces z=f(x)
        self.z = np.zeros(x.shape);  # flat bed - can enter own function here, z=f(x)
        self.zc = np.zeros(x.shape[0] - 1)
        self.dz = np.zeros(x.shape[0] - 1)
        # Compute bed slope
        for i in range(0, nc):
            self.dz[i] = self.z[i + 1] - self.z[i];  # dimensions of length
            self.zc[i] = 0.50 * (self.z[i + 1] + self.z[i]);  # elevation of cell center

        self.grav = 9.806

        #  initial conditions
        xo = L / 2
        ulo = 0
        uro = 0

        self.h = np.zeros_like(xc)
        self.u = np.zeros_like(xc)
        self.eta = 10 - np.arange( 0, 10.,  .10)
        for i in range(0, nc):
            if xc[i] < xo:
                self.h[i] = self.eta[i] - self.zc[i]
                self.u[i] = ulo
            else:
                self.h[i] = self.eta[i] - self.zc[i]
                self.u[i] = uro

        self.uh = self.h*self.u
        self.deta = np.zeros_like(self.eta)
        self.du = np.zeros_like(self.u)

        self.t = 0

        self.h_plot = np.zeros([nt/ntplot + 1, nc])

        for n in range(nt):
            self.deta = self.limiter(self.eta)
            self.du = self.limiter(self.u)
            self.etap, self.up = self.predictor()
            self.hp = self.etap - self.zc
            self.S = 9.806 * self.hp * self.dz / self.dx
            self.F, self.amax = self.fluxes(n)
            self.uh, self.h, self.u = self.corrector()
            self.eta = self.h + self.zc
            #e = self.eta + 0.5*self.u** 2
            self.t += self.dt
            cr = self.amax * dt/dx
            if cr > 1:
                return

            if (np.mod(n,ntplot) == 0):
                self.h_plot[n/ntplot]  = self.h

    def limiter(self, f):
        df = np.zeros_like(f)
        for i in range(1, f.shape[0] - 1):
            df1 = f[i+1] - f[i]
            df2 = f[i] - f[i-1]
            if (df1*df2) < 0:
                df[i] = 0
            else:
                s = np.sign(df1)
                a = abs(df1)
                b = abs(df2)
                df[i] = s * np.min([np.max([a, b]), self.beta * np.min([a, b])])
        return df

    def predictor(self):
        dh = self.deta - self.dz
        etap = np.zeros_like(self.eta)
        up = np.zeros_like(self.u)
        for i in range(self.nc):
            etap[i] = self.eta[i] - 0.5*self.dt/self.dx*(self.h[i]*self.du[i] + self.u[i]*dh[i])
            up[i] = self.u[i] - 0.5*self.dt/self.dx*(self.u[i] * self.du[i] + 9.806 * self.deta[i])
        return etap, up

    def solver(self, hl, hr, ul, ur, vl, vr, sn, cn):
        grav = 9.806
        duml = hl**0.5
        dumr = hr**0.5
        cl = (grav*hl)**0.5
        cr = (grav*hr)**0.5

        hhat = duml*dumr
        uhat = (duml*ul + dumr*ur)/(duml+dumr)
        vhat = (duml*vl + dumr*vr)/(duml+dumr)

        chat = (0.5*grav*(hl + hr))**0.5

        uperp = uhat*cn + vhat*sn
        dh = hr - hl
        du = ur - ul
        dv = vr - vl
        dupar = -du*sn + dv*cn
        duperp = du*cn + dv*sn
        dW = np.array([0.5 * (dh - hhat * duperp / chat), hhat * dupar, 0.5 * (dh + hhat * duperp / chat)])
        dW = np.expand_dims(dW, 1)
        uperpl = ul * cn + vl * sn
        uperpr = ur * cn + vr * sn

        al1 = uperpl - cl
        al3 = uperpl + cl
        ar1 = uperpr - cr
        ar3 = uperpr + cr

        R = np.array([1., 0., 1.,
                      uhat - chat*cn, -sn, uhat + chat*cn,
                      vhat - chat*sn, cn, vhat + chat*sn]).reshape(3, 3)

        da1 = np.max([0, 2 * (ar1 - al1)])
        da3 = np.max([0, 2 * (ar3 - al3)])
        a1 = np.abs(uperp - chat)
        a2 = np.abs(uperp)
        a3 = np.abs(uperp + chat)

        # Critical flow fix
        if a1 < da1:
            a1 = 0.5 * (a1*a1/da1 + da1)

        if a3 < da3:
            a3 = 0.5 * (a3*a3/da3 + da3)

        # Compute interface flux
        A = np.diag([a1, a2, a3])
        FL = np.array([uperpl*hl, ul * uperpl * hl + 0.5 * grav * hl * hl * cn, vl * uperpl * hl + 0.5 * grav * hl * hl * sn]);
        FR = np.array([uperpr*hr, ur * uperpr * hr + 0.5 * grav * hr * hr * cn, vr * uperpr * hr + 0.5 * grav * hr * hr * sn]);
        F = 0.5 * (FL + FR - np.dot(R, np.dot(A, dW)).ravel())
        amax = chat + abs(uperp)

        return F[:2], amax

    def fluxes(self, n):
        sn = 0
        cn = 1
        vl = 0
        vr = 0

        F = np.zeros([self.nf, 2])
        # Left Boundary : model as wall
        hr = self.etap[0] - 0.5 * self.deta[0] - self.z[0]
        ur = self.up[0] - 0.5 * self.du[0]
        F[0, :], amax = self.solver(hr, hr, -ur, ur, vl, vr, sn, cn)

        # Right Boundary : model as wall
        hl = self.etap[-1] + 0.5*self.deta[-1]-self.z[-1]
        ul = self.up[-1] + 0.5 * self.du[-1]
        F[-1, :], amax = self.solver(hl, hl, ul, -ul, vl, vr, sn, cn)

        amax0 = 0
        for j in range(1, self.nf - 1):
            # Variable reconstruction
            hl = self.etap[j-1] + 0.5*self.deta[j-1] - self.z[j-1] # Left h(j-1/2, n+1/2)=eta(j-1, n+1/2)+0.5 d_eta(j-1, n)-z(j-1/2, n+1/2)
            ul = self.up[j-1] + 0.5*self.du[j-1]  # u(j-1/2) = u(i-1) + 0.5*du(i-1)
            hr = self.etap[j] -0.5*self.deta[j] - self.z[j-1]  # Right h(j-1/2)=eta(j)-0.5 d_eta(j)-z(j-1/2)
            ur = self.up[j] - 0.5*self.du[j]
            # Call solver
            F[j, :], amax = self.solver(hl, hr, ul, ur, vl, vr, sn, cn)
            amax0 = max([amax0, amax])  # Keep track of max wave speed to check CFL

        return F, amax0

    def corrector(self):
        h = self.h.copy()
        uh = self.uh.copy()
        u = self.u.copy()

        for j in range(self.nc):
            h[j] = self.h[j] - self.dt/self.dx*(self.F[j+1, 0] - self.F[j, 0])
            uh[j] = self.uh[j] - self.dt/self.dx*(self.F[j+1,1] - self.F[j,1]) - self.dt*self.S[j]
            u[j] = uh[j]/h[j]  # update for dry bed cases
            # Add vfr for dry bed cases to compute eta from h
        return uh, h, u


f = fvm1d(nt = 2, ntplot=10)

import scipy.io as sio

sio.savemat("../fvm1d_wetbed/h_py.mat", {"h_py": f.h, "u_py": f.u, "eta_py": f.eta,  "deta_py" : f.deta, "du_py" : f.du,
                                         "F_py" : f.F, "etap_py" : f.etap, })

sio.savemat("../fvm1d_wetbed/inputs_py.mat", {"F_py": F})


