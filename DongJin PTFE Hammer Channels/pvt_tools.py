

import numpy as np
import matplotlib.pyplot as plt


class PVT(object):

    def __init__(self, fileName = 'a.txt'):
        self.fileName = fileName
        self.data = []
        self.P = []
        self.V = []
        self.A = []
        self.time = []


    @property
    def t(self):
        return self.data[:,0]

    @property
    def ct(self):
        '''cumulative time'''
        return self.t.cumsum()

    @property
    def x(self):
        return self.data[:,1].cumsum()

    @property
    def y(self):
        return self.data[:,3].cumsum()

    @property
    def z(self):
        return self.data[:,5].cumsum()

    @property
    def vx(self):
        return self.data[:,2]

    @property
    def vy(self):
        return self.data[:,4]

    @property
    def vz(self):
        return self.data[:,6]

    @property
    def s(self):
        return self.data[:,7].cumsum()

    @property
    def vs(self):
        return self.data[:,8]


    def read(self):
        """
        Reads the contents of the PVT file fileName
        :param fileName: name of the file to read
        :return: error code
        """
        dat = np.loadtxt(self.fileName)
        ini = np.zeros((1, dat.shape[1]))
        self.data = np.concatenate((ini, dat))
        return 0


    def segment(self, index):

        '''

        :param index:
        :return:
        '''

        #TODO: check inde


        npts = 100;
        p_new = self.data[index, :]
        p_old = self.data[index-1, :]
        DT = p_new[0]
        DX = p_new[[1, 3, 5, 7]]
        Vin = p_old[[2, 4, 6, 8]]
        Vout = p_new[[2, 4, 6, 8]]
        t = np.linspace(0, DT, npts).reshape(npts, 1)

        # Profile coefficients
        Jerk = (6. * (DT * (Vin + Vout) - 2. * DX)) / (DT** 3.)
        Jerk = Jerk.reshape(Jerk.shape[0],1).transpose()
        Gin = 2. * (3. * DX - DT * (2. * Vin + Vout)) / (DT** 2.)
        Gin = Gin.reshape(Gin.shape[0], 1).transpose()
        Gout = 2. * (DT * (Vin + 2. * Vout) - 3. * DX) / (DT** 2.)
        Gout = Gout.reshape(Gout.shape[0], 1).transpose()

        # Profile equations
        A = Gin + np.matmul(t, Jerk)
        V = Vin + t * Gin + np.matmul(t**2, Jerk) / 2.
        P = t * Vin + np.matmul(t**2, Gin) / 2. + np.matmul(t**3, Jerk) / 6.

        return (P, V, A, t)


    def plot(self):
        plt.plot(self.x, self.y, color='0.8')
        for c in range(1, self.data.shape[0]):
            cP, cV, cA, cT = self.segment(c)
            x = cP[:,0] + self.x[c-1]
            y = cP[:,1] + self.y[c-1]
            s = cP[:,3] + self.s[c-1]
            t = cT[:] + self.ct[c-1]

            ind = np.where(s < 1.)
            if ind[0].shape[0] > 0:
                plt.plot(x[ind], y[ind], color='y')

            ind = np.where(s >= 1.)
            if ind[0].shape[0] > 0:
                plt.plot(x[ind], y[ind], 'ko')

            if c == 1:
                self.P = cP
                self.V = cV
                self.A = cA
                self.time = t
            else:
                self.P = np.concatenate((self.P, cP))
                self.V = np.concatenate((self.V, cV))
                self.A = np.concatenate((self.A, cA))
                self.time = np.concatenate((self.time, t))

        plt.show()








if __name__ == "__main__":

    p = PVT('a.txt')
    p.read()
    p.plot()
    #plt.plot(p.time, p.A)
    plt.figure()
    plt.plot(p.time, p.A[:,0])

    plt.figure()
    plt.plot(p.A[:,2])
    plt.figure()
    plt.plot(p.time,p.A[:,1])


