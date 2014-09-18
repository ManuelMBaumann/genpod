import numpy as np
from scipy.integrate import odeint
import pod_utils as pdu
# import matplotlib.pyplot as plt
import matlibplots.conv_plot_utils as cpu
import json


def solve_podapprx(M, A, Y, f, B, poddim, k, tmesh, genpod=True):
    y0 = Y[:, 0]
    (A_red, y_red, rhs_red,
     Uk) = pdu.get_podred_model(M, A, Y, y0, f, poddim, k, tmesh=tmesh,
                                genpod=genpod, plotsvs=False, verbose=True)

    def f_pod(y, t):
        return (-np.dot(A_red, y).flatten() + rhs_red.flatten())

    def Df_pod(y, t):
        return -A_red

    # Solve reduced model:
    Yred = odeint(f_pod, y_red, tmesh, Dfun=Df_pod)
    Yappr = np.dot(Uk, Yred.T)

    return Yappr


if __name__ == '__main__':
    poddiml = [1, 3, 5, 10, 15, 20, 25]
    nsnapsl = [17, 33, 65]  # , 129]  # number of snapshots

    N = 25
    Nts = 5e2 + 1
    t0, tE = 0.0, 5.0
    Re = 2e3

    solvtol = 1e-4
    poddiml = [2**x for x in range(2, 6)]

    M, A, Ycheck, f, B, tmesh, J = \
        pdu.get_mayf(N=N, Re=Re, t0=t0, tE=tE, Nts=Nts,
                     # krylov='gmres', krpslvprms={'tol': solvtol},
                     krplsprms={})

    M, A, Y, f, B, tmesh, J = pdu.get_mayf(N=N, Re=Re, t0=t0, tE=tE, Nts=Nts)

    # ## check error vs poddims and nsnapshots
    datajsn = 'data/comppodgenpodRe{0}Nts{1}N{2}'.format(Re, Nts, N) +\
        'nsnaps{0}poddims{1}'.format(''.join(map(repr, nsnapsl)),
                                     ''.join(map(repr, poddiml)))
    try:
        print 'try to load data from ', datajsn
        fjs = open(datajsn)
        allersl = json.load(fjs)['allersl']
        fjs = open(datajsn)
        alllegsl = json.load(fjs)['alllegsl']
    except IOError:
        allersl, alllegsl = [], []
        for nsnaps in nsnapsl:
            gpderr, pderr = [], []
            for poddim in poddiml:
                # distance between and the indices of the snapshots
                ssdist = Nts/nsnaps
                ssinds = np.arange(0, Nts, ssdist).astype(int)

                gYappr = solve_podapprx(M, A, Ycheck, f, B, poddim, nsnaps,
                                        tmesh, genpod=True)
                Yappr = solve_podapprx(M, A, Ycheck[:, ssinds], f, B,
                                       poddim, nsnaps, tmesh, genpod=False)
                gerrM = pdu.timespace_diff_norm(tmesh, Y, gYappr, M)
                errM = pdu.timespace_diff_norm(tmesh, Y, Yappr, M)
                gpderr.append(gerrM)
                pderr.append(errM)
            allersl.extend([pderr, gpderr])
            alllegsl.extend(['pod, $k={0}$'.format(nsnaps),
                             'genpod, $k={0}$'.format(nsnaps)])

        jsfile = open(datajsn, mode='w')
        jsfile.write(json.dumps(dict(alllegsl=alllegsl,
                                     allersl=allersl)))
        jsfile.close()
        print 'saved pod genpod convdata to ' + datajsn

    cpu.para_plot(poddiml, allersl, leglist=alllegsl,
                  logscaley=10, tikzfile='ksnpoddimgpdpd.tikz')

