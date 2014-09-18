import numpy as np
import matplotlib.pyplot as plt

def uBasPLF(n, x0, xe, N):

    x = np.linspace(x0, xe, N)

    if n == 1:
        uBasPLF = (xe - x) / (xe - x0)
        return uBasPLF

    if n == 2:
        uBasPLF = (x0 - x) / (x0 - xe)
        return uBasPLF

    if n == 3:
        uBasPLF = np.interp(x, [x0, (x0 + xe) * 0.5, xe], [0, 1, 0])
        return uBasPLF

    l2 = np.floor(np.log2(n - 2))
    absInt = np.linspace(-x0, xe, 2 ** (l2 + 1) + 1)
    ordInt = np.zeros(2 ** (l2 + 1) + 1)
    ordInt[1] = 1

    index = (n - 2 - 2 ** l2) * 2
    index = int(index)
    ordInt = np.roll(ordInt, index)

    uBasPLF = np.interp(x, absInt, ordInt)

    return uBasPLF


def Haar_helper(x, x0, xe):
    if x >= x0 and x <= (xe - x0) / 2:
        return 1.0
    elif x > (xe - x0) / 2 and x <= xe:
        return -1.0
    else:
        return 0.0


def HaarWavelet(n, x0, xe, N):
    x = np.linspace(x0, xe, N)
    haar = np.ones((N, 1))
    
    if n == 1:
        return haar
    elif n == 2:
        haar[0:np.floor(N/2)] = 1
        haar[np.floor(N/2):N] = -1
        return haar
    else:
        l2 = np.floor(np.log2(n - 1))
        l3 = n - (2**l2) - 1
        
        for ii in range(0,N):
	    haar[ii] = np.sqrt(2.0**l2) * Haar_helper((2.0**l2)*x[ii] - (l3*(xe-x0)), x0, xe) #TODO: correct factor ?
        return haar


def massPL(Y, poddim, tmesh, haar = False):

    N = Y.shape[0]
    Nts = len(tmesh)

    Yg = np.zeros((N, poddim))
    My = np.zeros((poddim, poddim))
    NU = np.zeros((Nts, poddim))

    for j in range(0, poddim):
        if haar:
	    NU[:, j] = HaarWavelet(j + 1, tmesh[0], tmesh[-1], Nts).T
	else:
	    NU[:, j] = uBasPLF(j + 1, tmesh[0], tmesh[-1], Nts).T
        #plt.plot(tmesh, NU[:, j])
        #plt.show()

    for i in range(0, N):
        for j in range(0, poddim):
            Yg[i, j] = time_norm(tmesh, Y[i, :] * NU[:, j])

    for i in range(0, poddim):
        for j in range(0, i + 1):
            My[i, j] = time_norm(tmesh, NU[:, i] * NU[:, j])
            My[j, i] = My[i, j]
    
    return Yg, My


def time_norm(tmesh, Y):

    #dtc = tmesh[1] - tmesh[0]
    #err_old = 0.5 * dtc * Y[0]
    #err = err_old

    #for k in range(1, len(tmesh)):
        #dtc = tmesh[k] - tmesh[k - 1]
        #err_new = 0.5 * dtc * Y[k]
        #err += err_new
        #err_old = err_new
    dtvec = tmesh[1:] - tmesh[:-1]
    trapvec = 0.5*(Y[:-1] + Y[1:])
    err = (dtvec*trapvec).sum()
    
    return err

# if __name__ == '__main__':
#     get_mayf()
