"""
initcond4si3d.py
This script serves to create the files needed for the SI3D model runs. The code is based on previous matlab versions created by Alicia Cortes and Francisco Rueda.
Functions that are present within this script are:

1. initCond4si3d
    This function writes the initial condition file 'si3d_init.txt' for si3d simulations. The code considers constant and variable thickness layers, and the same for the temperature profiles. The use of the function is shown next for each of the scenarios. (4)
    Constant thickness and constant Temperature initCond4si3d(LakeName,SimStartDate,DeltaZ,TempProf,PathSave,NTracers,H,dz,Tc)
    Constant thickness and variable temperature profile initCond4si3d(LakeName,SimStartDate,DeltaZ,TempProf,PathSave,NTracers,H,dz,z_CTD,T_CTD)
    Variable thickness and constant Tempretaure initCond4si3d(LakeName,SimStartDate,DeltaZ,TempProf,PathSave,NTracers,H,dz,Tc,spacingMethod,dz0s,dz0b,dzxs,dzxb)
    Variable thickness and variable temperature profile initCond4si3d(LakeName,SimStartDate,DeltaZ,TempProf,PathSave,NTracers,H,dz,z_CTD,T_CTD,spacingMethod,dz0s,dz0b,dzxs,dzxb)

    NOTE PLEASE keep in mind that for the use of variable temperature profile, the CTD information must cover the whole lake depth. In case the CTD profile is incomplete (z_CTD[-1] < H) then it must be rearranged prior to using this function.

2. LayerGenerator
    This function writes the layer file for the initial conditions 'si3d_layer.txt' for si3d simulations. The code is only used when the option of variable thickness of layers is toggled. This function writes the number of layers and the depth at each layer, which is needed when ibathyf < 0

For a better understanding on the use of these functions, the reader is directed to the corresponding repositories that make use of the functions in here. "surfBondCond.py", InitConditions.py", and ""bathymetry.py"

Copy right Sergio A. Valbuena 2021
UC Davis - TERC
February 2021
"""
import os
import numpy as np


def initCond4si3d(LakeName, SimStartDate, DeltaZ, TempProf, PathSave, NTracers, **kw):
    """
    Creates an initial condition file for SI3D
    :param LakeName:
    :param SimStartDate:
    :param DeltaZ:
    :param TempProf:
    :param PathSave:
    :param NTracers:
    :param kw:
    :return:
    """

    if DeltaZ:
        if TempProf:
            z = np.arange(0 + kw['dz'] / 2, kw['H'], kw['dz'])
            T = kw['Tc'] * np.ones(len(z))
            z *= -1
            dummy2 = 'Source: From constant values                     - '
        else:
            z = np.arange(0 + kw['dz'] / 2, kw['H'] + kw['dz'], kw['dz'])
            T = np.interp(z, kw['z_CTD'], kw['T_CTD'])
            z *= -1
            dummy2 = 'Source: From CTD_Profile                         - '
        dummy1 = 'Depths (m) not used   Temp (oC)                  - '
        if NTracers != 0:
            dummy1 = 'Depths (m) not used;  Temp (oC);  WQ [mg/m3], SS [kg/m3], Hg [ng/m3], Tracers (M/V) --> - '
    else:
        # Length of initial grid for creating the unevenly spaced grid
        N = 1000
        gridks = np.arange(1, N + 1, 1)
        if TempProf:
            if kw['spacingMethod'] == 'exp':
                gridDZ = kw['dz0s'] * kw['dzxs'] ** gridks
                gridZ = np.cumsum(gridDZ)
                igdepth = np.where(gridZ >= kw['H'])
                km = igdepth[0][0]
                kml = km + 2
                surf = np.array([-100, -100])
                zlevel = np.concatenate((surf, gridZ[0:km + 1]))
                _ = LayerGenerator(zlevel, kml, PathSave)
                zz = np.zeros(len(zlevel) - 1)
                zz[1:] = zlevel[2:]
                zi = -(zz[0:-1] + zz[1:]) / 2
                z = zi
            elif kw['spacingMethod'] == 'sbconc':
                gridkb = np.arange(1, N + 1, 1)
                gridDZs = kw['dz0s'] * kw['dzxs'] ** gridks
                gridZs = np.cumsum(gridDZs)
                gridDZb = kw['dz0b'] * kw['dzxb'] ** gridkb
                gridZb = np.zeros(len(gridDZb) + 1)
                gridZb[0] = kw['H']
                gridZb[1:] = kw['H'] - np.cumsum(gridDZb)
                idepths = gridZs <= kw['Hn']
                idepthb = gridZb >= kw['Hn']
                gridZbb = gridZb[idepthb]
                gridZss = gridZs[idepths]
                gridZbb = sorted(gridZbb)
                gridZ = np.concatenate((gridZss, gridZbb))
                km = len(gridZ)
                kml = km + 2
                surf = np.array([-100, -100])
                zlevel = np.concatenate((surf, gridZ))
                _ = LayerGenerator(zlevel, kml, PathSave)
                zz = zlevel[1:]
                zz[0] = 0
                zi = -(zz[0:-1] + zz[1:]) / 2
                z = zi
            elif kw['spacingMethod'] == 'surfvarBotconsta':
                gridkb = np.arange(1, N + 1, 1)
                gridDZs = kw['dz0s'] * kw['dzxs'] ** gridks
                gridZs = np.cumsum(gridDZs)
                idepths = gridZs <= kw['Hn']
                gridZss = gridZs[idepths]
                Href = gridZss[-1]
                gridZb = np.arange(Href + kw['dzc'], kw['H'] + kw['dzc'], kw['dzc'])
                if gridZb[-1] > kw['H']:
                    gridZbb = gridZb
                    gridZbb[-1] = kw['H']
                else:
                    gridZbb = np.concatenate((gridZb, kw['H']))
                gridZ = np.concatenate((gridZss, gridZbb))
                gridZ = np.round(gridZ, 2)
                km = len(gridZ)
                kml = km + 2
                surf = np.array([-100, -100])
                zlevel = np.concatenate((surf, gridZ))
                _ = LayerGenerator(zlevel, kml, PathSave)
                zz = zlevel[1:]
                zz[0] = 0
                zi = -(zz[0:-1] + zz[1:]) / 2
                z = zi
            T = np.interp(z, kw['z_CTD'], kw['T_CTD'])
            dummy2 = 'Source: From constant values                     - '
        else:
            if kw['spacingMethod'] == 'exp':
                gridDZ = kw['dz0s'] * kw['dzxs'] ** gridks
                gridZ = np.cumsum(gridDZ)
                igdepth = np.where(gridZ >= kw['H'])
                km = igdepth[0][0]
                kml = km + 2
                surf = np.array([-100, -100])
                zlevel = np.concatenate((surf, gridZ[0:km + 1]))
                _ = LayerGenerator(zlevel, kml, PathSave)
                zz = np.zeros(len(zlevel) - 1)
                zz[1:] = zlevel[2:]
                zi = -(zz[0:-1] + zz[1:]) / 2
                z = zi
            elif kw['spacingMethod'] == 'sbconc':
                gridkb = np.arange(1, N + 1, 1)
                gridDZs = kw['dz0s'] * kw['dzxs'] ** gridks
                gridZs = np.cumsum(gridDZs)
                gridDZb = kw['dz0b'] * kw['dzxb'] ** gridkb
                gridZb = np.zeros(len(gridDZb) + 1)
                gridZb[0] = kw['H']
                gridZb[1:] = kw['H'] - np.cumsum(gridDZb)
                idepths = gridZs <= kw['H'] * (1 - 1 / kw['n'])
                idepthb = gridZb >= kw['H'] * (1 - 1 / kw['n'])
                gridZbb = gridZb[idepthb]
                gridZss = gridZs[idepths]
                gridZbb = sorted(gridZbb)
                gridZ = np.concatenate((gridZss, gridZbb))
                km = len(gridZ)
                kml = km + 2
                surf = np.array([-100, -100])
                zlevel = np.concatenate((surf, gridZ))
                _ = LayerGenerator(zlevel, kml, PathSave)
                zz = zlevel[1:]
                zz[0] = 0
                zi = -(zz[0:-1] + zz[1:]) / 2
                z = zi
            elif kw['spacingMethod'] == 'surfvarBotconsta':
                gridkb = np.arange(1, N + 1, 1)
                gridDZs = kw['dz0s'] * kw['dzxs'] ** gridks
                gridZs = np.cumsum(gridDZs)
                gridDZb = kw['dz0b'] * kw['dzxb'] ** gridkb
                gridZb = np.zeros(len(gridDZb) + 1)
                gridZb[0] = kw['H']
                gridZb[1:] = kw['H'] - np.cumsum(gridDZb)
                idepths = gridZs <= kw['Hn']
                gridZss = gridZs[idepths]
                Href = gridZss[-1]
                gridZb = np.arange(Href + kw['dzc'], kw['H'] + kw['dzc'], kw['dzc'])
                if gridZb[-1] > kw['H']:
                    gridZbb = gridZb
                    gridZbb[-1] = kw['H']
                else:
                    gridZbb = np.concatenate((gridZb, kw['H']))
                gridZ = np.concatenate((gridZss, gridZbb))
                gridZ = np.round(gridZ, 2)
                km = len(gridZ)
                kml = km + 2
                surf = np.array([-100, -100])
                zlevel = np.concatenate((surf, gridZ))
                _ = LayerGenerator(zlevel, kml, PathSave)
                zz = zlevel[1:]
                zz[0] = 0
                zi = -(zz[0:-1] + zz[1:]) / 2
                z = zi
            z *= -1
            T = np.interp(z, kw['z_CTD'], kw['T_CTD'])
            z *= -1
            dummy2 = 'Source: From CTD_Profile                         - '
        dummy1 = 'Depths (m)   Temp (oC)                           - '

    dummy3 = '         z          T'
    if NTracers != 0:
        dummy1 = 'Depths (m) not used;  Temp (oC);  WQ [mg/m3], SS [kg/m3], Hg [ng/m3], Tracers (M/V) --> - '
        # To create array of header names for z, T, and constituents
        name_Tr = kw['name_Tr']
        for i in range(0, len(name_Tr)):
            while len(name_Tr[i]) <= 10:
                name_Tr[i] = ' ' + name_Tr[i]
        dummy3 = dummy3 + ''.join(name_Tr)
    # ----------------------- Creation of file ---------------------------------
    os.chdir(PathSave)
    fid = open('si3d_init.txt', 'w+')
    fid.write('%s\n' % 'Initial condition file for si3d model            - ')
    fid.write('%s' % LakeName + '             - ' + '\n')
    fid.write('%s' % 'Simulation starting on ' + SimStartDate + ' UTC    - ' + '\n')
    fid.write('%s\n' % dummy1)
    fid.write('%s\n' % dummy2)
    fid.write('%s\n' % dummy3)
    fid.write('%s\n' % '----------------------------------------------------------------- ')

    if NTracers == 0:
        fid.write('%10.2f %10.4f \n' % (z[0], T[0]))
        for i in range(0, len(T)):
            fid.write('%10.2f %10.4f \n' % (z[i], T[i]))
        fid.write('%10.2f %10.4f \n' % (z[i], T[i]))
    else:
        _, cols1 = kw['z_Tr'].shape
        _, cols2 = kw['conc_Tr'].shape
        if cols1 != NTracers or cols2 != NTracers or cols1 != cols2:
            print('¡¡¡ERROR!!! The number of tracers ', str(NTracers),
                  ' does not match number of columns in the variables for the depths ', str(cols1),
                  ' and concentrations of the tracers ', str(cols2))
            exit()
        tracers = np.empty((len(z), NTracers)) * np.nan
        for i in range(0, NTracers):
            tracers[:, i] = np.interp(-z, kw['z_Tr'][:, i], kw['conc_Tr'][:, i])
        fid.write('%10.2f %10.4f' % (z[0], T[0]))
        for j in range(0, NTracers):
            fid.write('%11.4f' % tracers[0, j])
        fid.write('\n')
        for i in range(0, len(T)):
            fid.write('%10.2f %10.4f' % (z[i], T[i]))
            for j in range(0, NTracers):
                fid.write('%11.4f' % tracers[i, j])
            fid.write('\n')
        fid.write('%10.2f %10.4f' % (z[i], T[i]))
        for j in range(0, NTracers):
            fid.write('%11.4f' % tracers[i, j])
        fid.write('\n')

    fid.close()
    return T, z


def LayerGenerator(zlevel, kml, PathSave):
    """
    This function is only used when the layer thickness is variable
    :param zlevel:
    :param kml:
    :param PathSave:
    :return:
    """
    os.chdir(PathSave)
    fid = open('si3d_layer.txt', 'wt+')
    fid.write('%s\n' % 'Depths to top of layers in Si3D Grid            ')
    fid.write('%s\n' % '** used if ibathyf in si3d_inp.txt is set to < 0       ')
    fid.write('%s\n' % '------------------------------------------------------ ')
    fid.write('%s' % '   km1   =        ' + str(kml) + '\n')
    for i in range(0, kml):
        fid.write('%10.2f %10.4f \n' % (i + 1, zlevel[i]))

    fid.close()
    Layer = 'Layer file created in the folder' + PathSave + ' with name si3d_layer.txt'
    print(Layer)
    return Layer
