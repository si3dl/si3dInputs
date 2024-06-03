"""
# Copy right Sergio A. Valbuena 2021
# UC Davis - TERC
# February 2023
"""
import os
import numpy as np


def openbc4si3d(idx, simStartDate, t, in_out, bc_type, pathSave, LakeName, **kw):
    """
    Creates a open boundary condition input file for SI3D.
    :param idx: Integer identifying id of open boundary.
    :type idx: int
    :param simStartDate: String identifying the date of start for simulation.
    :type simStartDate: str
    :param t: Time in hours for time series records
    :param in_out: Flag specifying type of open boundary. Inflow = 1, Outflow = 0
    :param bc_type: Flag identifying open boundary type. Water Surface Elevation = 1, Surface Flow = 2, Subsurface Flow = 3
    :param pathSave:
    :param kw:
    :return:
    """

    r = len(t)

    # ----------------------- Creation of file ---------------------------------
    os.chdir(pathSave)
    filename = 'openbc0' + str(idx) + '.txt'
    fid = open(filename, 'w+')
    fid.write('%s' % 'Open boundary condition file for si3d model' + LakeName + ' -\n')
    fid.write('%s' % 'Time is in hours from day ' + simStartDate + ' UTC    - ' + '\n')

    if in_out == 1:
        fid.write('%s\n' % 'Data for Inflow open boundary condition')
        if bc_type == 1:
            fid.write('%s\n' % '(1) Time [hrs],  (2) wse [m], (3) T [oC], (4...) Tracers ')
            v2 = kw['wse']
            v3 = kw['T']

        elif bc_type == 2:
            fid.write('%s\n' % '(1) Time [hrs],  (2) Qin [m3/s], (3) T [oC], (4...) Tracers ')
            v2 = kw['Qin']
            v3 = kw['T']

        elif bc_type == 3:
            print('UNDER DEVELOPMENT. FUNCTION NOT WORKING AT THE MOMENT. PLEASE USE 1 OR 2')
            exit()

        fid.write('%s\n' % '-------------------------------------------------- ')
        fid.write('%s' % '   npts  =   ' + str(r) + '\n')

    elif in_out == 0:
        fid.write('%s\n' % 'Data for Outflow open boundary condition')
        if bc_type == 1:
            fid.write('%s\n' % '(1) Time [hrs],  (2) wse [m],  (3) T [oC],  (4...) Tracers ')
            v2 = kw['wse']
            v3 = -0.01 * np.ones((r, 1))

        elif bc_type == 2:
            fid.write('%s\n' % '(1) Time [hrs],  (2) Qout [m3/s],  (3) T [oC],  (4...) Tracers ')
            v2 = kw['Qout']
            v3 = -0.01 * np.ones((r, 1))

        elif bc_type == 3:
            print('UNDER DEVELOPMENT. FUNCTION NOT WORKING AT THE MOMENT. PLEASE USE 1 OR 2')
            exit()

        fid.write('%s\n' % '-------------------------------------------------- ')
        fid.write('%s' % '   npts  =   ' + str(r) + '\n')
    else:
        print('¡¡¡ERROR!!! The flag for inflow or outflow was not specified correctly. Please use 1 for inflow and 0 for outflow')
        exit()

    if 'concTr' in kw:
        rows, nTracer = np.shape(kw['concTr'])
        if in_out == 1:
            v4 = kw['concTr']
        elif in_out == 0:
            v4 = -0.01 * np.ones((rows, nTracer))
        for i in range(0, r):
            fid.write('%10.4f %10.4f %10.4f' % (t[i], v2[i], v3[i]))
            for j in range(0, nTracer):
                fid.write('%10.4f' % v4[i, j])
            fid.write('\n')
    else:
        for i in range(0, r):
            fid.write('%10.4f %10.4f %10.4f' % (t[i], v2[i], v3[i]))
            fid.write('\n')

    fid.close()
    return
