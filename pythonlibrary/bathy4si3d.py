"""
bathy4si3d.py
This script serves to create the files needed for the SI3D model runs. The code is based on previous matlab versions created by Alicia Cortes and Francisco Rueda.
Functions that are present within this script are:
1. bathy4si3d
    This function writes the bathymetry file 'h' for si3d simulations. The code considers canonical and real basins. The use of this function is found next for each of the basins considered:
    if the basin is a real lake use the functions as: bathy4si3d(BasinType,SimName,dx,xg,yg,zg)
    where xg, yg, and zg are 2-D matrices that contain the grid dimensions for the horizontal dimension based on a x=0,y=0 origin, and for the vetical dimension uses the depth of the lake with origin z=0 at the lake's surface amd NEGATIVE z values.
    if the basin is rectangular use the functions as: bathy4si3d(BasinType,SimName,dx,L,B,H)
    if the basin is spherical use the functions as: bathy4si3d(BasinType,SimName,dx,D,H)
    if the basin is cylindrical use the functions as: bathy4si3d(BasinType,SimName,dx,D,H)
    Where L,B,H are the dimensions of length, width, and depth for the rectangular basin. D,H are the diameter and depth respectively for the cylindrical and spherical basins.
    The file name will have the description of the grid size dx and the type of basin. This is for referencing but the user must change this name to 'h' to be able to run simulations in Si3D model

For a better understanding on the use of these functions, the reader is directed to the corresponding repositories that make use of the functions in here. "surfBondCond.py", InitConditions.py", and ""bathymetry.py"

Copy right Sergio A. Valbuena 2021
UC Davis - TERC
February 2021
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as Dt


def bathy4si3d(BasinType, SimName, dx, PathSave, *args):
    """
    Creates a bathymetry file for SI3D.
    :param BasinType: Integer corresponding to basin type 1: Lake, type 2: rectangular, or type 3: circular.
    :type BasinType: int [1-3]
    :param SimName:
    :type SimName: str
    :param dx:
    :param PathSave:
    :param args:
    :return: (x, y, z) tuple of numpy 2D meshgrids corresponding to x, y, z coordinates of bathymetry data
    """
    dxsave = ' (dx= ' + str(dx) + '),'
    Entry = SimName + dxsave
    if len(Entry) != 27:
        print(
            'ERROR! The length of the header for the bathymetry files must be 27 characters. The current number is ' + str(
                len(Entry)) + ' characters. Please change the SimName accorddingly to match the length')
        sys.exit()
    else:
        print('GOOD! The length of the header is the right length (27)')
    if BasinType == 1:
        basin = 'Lake'
        mindepth = 0
        X = args[0]
        Y = args[1]
        zg = args[2]
        zg = np.flipud(zg)
        idata = zg > mindepth
        zg[idata] = mindepth
        zz = -99 * np.ones(np.shape(zg))
        idata = ~np.isnan(zg)
        zz[idata] = zg[idata] * (-10)
        Z = zz
    elif BasinType == 2:
        basin = 'rectangular'
        L = args[0]
        B = args[1]
        H = args[2]
        x = np.arange(dx, L + dx, dx)
        y = np.arange(dx, B + dx, dx)
        X, Y = np.meshgrid(x, y)
        z = -99 * np.ones(np.shape(X))
        z[:, :] = H * (10)
        Z = z
        del x, y, z
    elif BasinType == 3:
        basin = 'circular'
        R = args[0] / 2
        H = args[1]
        x = np.arange(0, 2 * R + 2 * dx, dx)
        y = np.arange(0, 2 * R + 2 * dx, dx)
        X, Y = np.meshgrid(x, y)
        z = np.empty((len(x), len(y)))
        C = H / R
        for i in range(len(x)):
            for j in range(len(x)):
                A = R ** 2 - (x[i] - R) ** 2 - (y[j] - R) ** 2
                if A < 0:
                    z[i, j] = np.nan
                else:
                    z[i, j] = C * A ** (1 / 2)
        Z = z[0:-1, 0:-1]
        Z = Z * 10
        idata = np.isnan(Z)
        Z[idata] = -99
    else:
        H = 1

    os.chdir(PathSave)
    ny, nx = np.shape(Z)
    filename = 'h'
    fid = open(filename, "w+")
    fid.write("%s" % Entry + '   imx =  ' + str(nx) + ',jmx =  ' + str(ny) + ',ncols = ' + str(nx))
    fid.write("\n")
    H1 = 'HV       V'
    for i in range(1, nx):
        H1 = H1 + '   V'

    fid.write("%s\n" % H1)
    fid.write('%s' % '     ')
    for i in range(1, nx):
        fid.write("%5.0f" % (i + 1))
    fid.write("%5.0f\n" % (nx + 1))
    for i in range(1, ny + 1):
        fid.write("%5d" % (ny - i + 2))
        for item in Z[i - 1, :]:
            fid.write("%5.0f" % item)
        fid.write('\n')

    fid.close()
    print('The bathymetry file was save in ' + PathSave + ' as ' + filename)
    return X, Y, Z