import dolfin
import os
import numpy as np

# import dolfin_navier_scipy.dolfin_to_sparrays as dts
import dolfin_navier_scipy.data_output_utils as dou
import dolfin_navier_scipy.stokes_navier_utils as snu
import dolfin_navier_scipy.problem_setups as dnsps

import distr_control_fenics.cont_obs_utils as cou


dolfin.parameters.linear_algebra_backend = 'uBLAS'


def gopod(problemname='drivencavity',
          N=10, Re=1e2, t0=0.0, tE=1.0, Nts=11, NU=3, NY=3,
          paraoutput=True, multiproc=False,
          krylov=None, krpslvprms={}, krplsprms={}):
    """Main routine for LQGBT

    Parameters
    ----------
    problemname : string, optional
        what problem to be solved, 'cylinderwake' or 'drivencavity'
    N : int, optional
        parameter for the dimension of the space discretization
    Re : real, optional
        Reynolds number, defaults to `1e2`
    t0, tE, Nts : real, real, int, optional
        starting and endpoint of the considered time interval, number of
        time instancses, default to `0.0, 1.0, 11`
    NU, NY : int, optional
        dimensions of components of in and output space (will double because
        there are two components), default to `3, 3`
    krylov : {None, 'gmres'}, optional
        whether or not to use an iterative solver, defaults to `None`
    krpslvprms : dictionary, optional
        to specify parameters of the linear solver for use in Krypy, e.g.,

          * initial guess
          * tolerance
          * number of iterations

        defaults to `None`
    krplsprms : dictionary, optional
        parameters to define the linear system like

          *preconditioner

    """
    femp, stokesmatsc, rhsd_vfrc, rhsd_stbc, data_prfx, ddir, proutdir = \
        dnsps.get_sysmats(problem=problemname, N=N, Re=Re)

    # specify in what spatial direction Bu changes. The remaining is constant
    uspacedep = femp['uspacedep']

    # output
    ddir = 'data/'
    try:
        os.chdir(ddir)
    except OSError:
        raise Warning('need "' + ddir + '" subdir for storing the data')
    os.chdir('..')
    data_prfx = ddir + data_prfx

    # casting some parameters
    NV = len(femp['invinds'])

    # contsetupstr = 'NV{0}NU{1}NY{2}alphau{3}'.format(NV, NU, NY, alphau)
    contsetupstr = 'NV{0}NU{1}NY{2}Re{3}'.format(NV, NU, NY, Re)

    soldict = stokesmatsc  # containing A, J, JT
    soldict.update(femp)  # adding V, Q, invinds, diribcs
    soldict.update(rhsd_vfrc)  # adding fvc, fpr
    soldict.update(fv_stbc=rhsd_stbc['fv'], fp_stbc=rhsd_stbc['fp'],
                   N=N, nu=femp['nu'], data_prfx=data_prfx)
    soldict.update(paraviewoutput=paraoutput)
    soldict.update(krylov=krylov, krplsprms=krplsprms, krpslvprms=krpslvprms)

#
# Prepare for control
#

    # get the control and observation operators
    try:
        b_mat = dou.load_spa(ddir + contsetupstr + '__b_mat')
        u_masmat = dou.load_spa(ddir + contsetupstr + '__u_masmat')
        print 'loaded `b_mat`'
    except IOError:
        print 'computing `b_mat`...'
        b_mat, u_masmat = cou.get_inp_opa(cdcoo=femp['cdcoo'], V=femp['V'],
                                          NU=NU, xcomp=uspacedep)
        dou.save_spa(b_mat, ddir + contsetupstr + '__b_mat')
        dou.save_spa(u_masmat, ddir + contsetupstr + '__u_masmat')
    try:
        mc_mat = dou.load_spa(ddir + contsetupstr + '__mc_mat')
        y_masmat = dou.load_spa(ddir + contsetupstr + '__y_masmat')
        print 'loaded `c_mat`'
    except IOError:
        print 'computing `c_mat`...'
        mc_mat, y_masmat = cou.get_mout_opa(odcoo=femp['odcoo'],
                                            V=femp['V'], NY=NY)
        dou.save_spa(mc_mat, ddir + contsetupstr + '__mc_mat')
        dou.save_spa(y_masmat, ddir + contsetupstr + '__y_masmat')

    # restrict the operators to the inner nodes
    invinds = femp['invinds']
    mc_mat = mc_mat[:, invinds][:, :]
    b_mat = b_mat[invinds, :][:, :]

    # tb_mat = 1./np.sqrt(alphau)

# setup the system for the correction
#
# # compute the uncontrolled steady state Stokes solution
#
    v_ss_stokes, list_norm_nwtnupd = \
        snu.solve_steadystate_nse(vel_pcrd_stps=0, vel_nwtn_stps=0,
                                  clearprvdata=True, **soldict)
    tmesh = np.linspace(t0, tE, Nts)

    soldict.update(trange=tmesh,
                   iniv=v_ss_stokes,
                   lin_vel_point=v_ss_stokes,
                   clearprvdata=True,
                   vel_nwtn_stps=1,
                   return_dictofvelstrs=False,
                   paraviewoutput=True,
                   vfileprfx='results/fullvel',
                   pfileprfx='results/fullp')

    convc_mat_n, rhs_con_n, rhsv_conbc_n = \
        snu.get_v_conv_conts(prev_v=v_ss_stokes, invinds=invinds,
                             V=femp['V'], diribcs=femp['diribcs'],
                             Picard=False)

    convc_mat_z, rhs_con_z, rhsv_conbc_z = \
        snu.get_v_conv_conts(prev_v=0*v_ss_stokes, invinds=invinds,
                             V=femp['V'], diribcs=femp['diribcs'],
                             Picard=False)

    vellist = snu.solve_nse(return_as_list=True, **soldict)
    velar = np.array(vellist)[:, :, 0].T
    rhsv = soldict['fv_stbc'] + soldict['fvc'] + rhsv_conbc_n + rhs_con_n
    rhsp = soldict['fp_stbc'] + soldict['fpr']

    # print 'fvstbc', np.linalg.norm(soldict['fv_stbc'])
    # print 'fvc', np.linalg.norm(soldict['fvc'])
    # print 'rhsvconbc', np.linalg.norm(rhsv_conbc_n)
    # print 'rhscon', np.linalg.norm(rhs_con_n)

    print 'velarshape :', velar.shape

    checkreturns = False
    if checkreturns:
        inivel = velar[:, 0:1]

        ylist = snu.solve_nse(A=soldict['A']+convc_mat_n-convc_mat_z,
                              M=soldict['M'],
                              J=soldict['J'], fvc=rhsv, fpr=rhsp,
                              iniv=inivel,
                              fv_stbc=0*rhsv - rhsv_conbc_z - rhs_con_z,
                              fp_stbc=0*rhsp,
                              lin_vel_point=0*inivel, trange=tmesh,
                              V=femp['V'], Q=femp['Q'],
                              invinds=femp['invinds'],
                              diribcs=femp['diribcs'], N=soldict['N'],
                              nu=soldict['nu'],
                              vel_nwtn_stps=1,
                              return_as_list=True)

        velarcheck = np.array(ylist)[:, :, 0].T

        print np.linalg.norm(velarcheck - velar)
        # print np.linalg.norm(velarcheck[:, 1] - velar[:, 1])

    return (soldict['M'], soldict['A']+convc_mat_n, velar,
            rhsv, b_mat, tmesh, soldict['J'])


if __name__ == '__main__':
    # lqgbt(N=10, Re=500, use_ric_ini=None, plain_bt=False)
    gopod(problemname='drivencavity', N=10,  # use_ric_ini=2e2,
          Re=1.0e2, t0=0.0, tE=2.0, Nts=1e3+1)
