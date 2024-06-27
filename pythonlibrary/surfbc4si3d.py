"""
surfbc4si3d.py
This script serves to create the files needed for the SI3D model runs. The code is based on previous matlab versions created by Alicia Cortes and Francisco Rueda.
Functions that are present within this script are:

1. Surfbc4si3d
    This function writes the surface boundary condition file 'surfbc.txt' for si3d simulations. The code includes 3 different methods to input forcing conditions that are coherent with si3d capabilities. 1) The heat budget is estimated among 3 different possibilities. The resulting net heat source is saved within the file along with other atmospheric conditions like wind speed, air temperature, atmospheric pressure, eta and others. 2) Using the shortwave date from site and uses cloud cover to obtain the needed parameters to run a heatbudget model built within si3d. 3) Uses the same input parameters are 2) with the difference that the longwave incoming radiation is used rather than the approximation used by using cloud cover estimations.
    The proper use of this function is shown next for each of the 3 possibilities:
    1) surfbc4si3d(LakeName,surfbcType,days,hr,mins,year,dt,PathSave,HeatBudgetMethod,eta,Hswn,Hlwin,Hlwout,Ta,Pa,RH,Cl,cw,u,v,WaTemp,Pa_P,esMethod)
    2) surfbc4si3d(LakeName,surfbcType,days,hr,mins,year,dt,PathSave,eta,Hswn,Ta,Pa,RH,Cl,cw,u,v)
    3) surfbc4si3d(LakeName,surfbcType,days,hr,mins,year,dt,PathSave,eta,Hswn,Ta,Pa,RH,Hlwin,cw,u,v)
    Where the definition of the variables is:
    u,v horizontal wind velocity components.
    Ta stands for air temperature, Pa for atmospheric pressure, RH relative humidity, eta the light penetration coefficient (secchi depth dependent), CL is cloud cover, WaTemp is the surface water temperature.
    Pa_P is the ratio of the atmospheric pressute at site in comparison to sea pressure, Hswn is the net shortwave radiation, Hlwin and HLwout stand for the incoming and outgoing longwave radiation. Finally, cw stands for wind drag coefficient.

For a better understanding on the use of these functions, the reader is directed to the corresponding repositories that make use of the functions in here. "surfBondCond.py", InitConditions.py", and ""bathymetry.py"

Copy right Sergio A. Valbuena 2021
UC Davis - TERC
February 2021
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as Dt


def surfbcW4si3d(caseStudy, Time, dt, PathSave, cw, u, v):
    """
    :param caseStudy:
    :param Time:
    :param dt:
    :param PathSave:
    :param cw:
    :param u:
    :param v:
    :return:
    """
    os.chdir(PathSave)
    r = len(Time)
    days = Time / 24
    daystart = Time[0]

    fid = open('surfbcW.txt', 'wt+')
    fid.write('%s\n' % 'Surface boundary condition file for si3d model')
    fid.write('%s' % caseStudy + ' simulations \n')
    fid.write('%s' % 'Time is given in hours from the start date used within the input.txt \n')
    fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time cw ua va')
    fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + caseStudy + ' Met Data \n')
    fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + '\n')
    fid.write('%s' % '   npts = ' + str(r) + '\n')
    for i in range(0, r):
        a0 = (days[i] - daystart) * 24
        a1 = cw[i]  # **** Wind drag coefficient
        a2 = u[i]  # **** Wind speed in the EW direction
        a3 = v[i]  # **** Wind speed in the NS direction
        format = '%10.4f %10.4f %10.4f %10.4f \n'
        fid.write(format % (a0, a1, a2, a3))
    return


def surfbc4si3d(show, LakeName, surfbcType, days, hr, mins, year, dt, PathSave, **kw):
    """
    Function to create surface boundary condition using variable options for the heat
    sources. This function preprocess the meteorological parameters and creates a si3d_surfbc.txt file
    for SI3D. The file has the inputs for the heatbudget method chosen.
    :param show:
    :param LakeName:
    :param surfbcType:
    :param days:
    :param hr:
    :param mins:
    :param year:
    :param dt:
    :param PathSave:
    :param args:
    :return:
    """
    os.chdir(PathSave)
    r = len(days)
    daystart = days[0]
    # To write the file surfbc for the numerical simulation in si3d
    fid = open('si3d_surfbc.txt', 'wt+')
    fid.write('%s\n' % 'Surface boundary condition file for si3d model')
    fid.write('%s' % LakeName + ' simulations \n')
    fid.write('%s' % 'Time is given in hours from ' + str(hr[0]) + ':' + str(mins[0]) + ' hrs on julian day ' + str(
        days[0]) + ',' + str(year) + '\n')

    if surfbcType == 1:
        HeatBudgetMethod = kw['HM']
        eta = kw['eta']
        Hswn = kw['Hswn']
        Hlwin = kw['Hlwin']
        Hlwout = kw['Hlwout']
        Ta = kw['Ta']
        Pa = kw['Pa']
        RH = kw['RH'] / 100
        Cl = kw['Cl']
        cw = kw['cw']
        u = kw['u']
        v = kw['v']
        WaTemp = kw['waT']
        TimeSim = kw['t_Sim']

        if HeatBudgetMethod == 'Chapra1995':
            rho0 = 1000
            Lw = 2.6e-6
            cChapra = args[13]
            CbPa_P = 0.61 * cChapra
            esMethod = args[14]
            # Vapor Pressure
            if esMethod == 1:
                es = 6.11 * np.exp(17.3 * Ta / (Ta + 237.3))
            elif esMethod == 2:
                es = 6.11 * np.exp(7.5 * Ta / (Ta + 237.3))
            elif esMethod == 3:
                es3 = 10 ** (9.286 - (2322.38 / (Ta + 273.15)))
            ea = es * RH
            # Longwave radiation
            Hlwn = Hlwin - Hlwout
            # Latent Heat Flux (Negative as it exits)
            fwind = 1.02e-9 * (u ** 2 + v ** 2) ** 0.5
            Hl = -rho0 * Lw * fwind * (es - ea)
            # Sensible heat flux (negative due to consideration of difference between water temp and air temp)
            Hs = -rho0 * Lw * fwind * CbPa_P * (WaTemp - Ta)
            # Net Heat Flux
            Hn = Hswn + Hlwn + Hs + Hl
        elif HeatBudgetMethod == 'AirSea':
            print('UNDER DEVELOPMENT, THIS FUNCTION DOES NOT WORK')
            print('The file has not been created')
            exit()
        elif HeatBudgetMethod == 'TERC':
            print('UNDER DEVELOPMENT, THIS FUNCTION DOES NOT WORK')
            print('The file has not been created')
            exit()
        fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time attc Hsw Hn cw ua va')
        fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + LakeName + ' Met Data ' + str(year) + '\n')
        fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + 'HeatBudget = ' + HeatBudgetMethod + '\n')
        fid.write('%s' % '   npts = ' + str(r) + '\n')
        for i in range(0, r):
            a0 = (days[i] - daystart) * 24
            a1 = eta[i]
            a2 = Hswn[i]
            a3 = Hn[i]
            a4 = cw[i]
            a5 = u[i]
            a6 = v[i]
            format = '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n'
            fid.write(format % (a0, a1, a2, a3, a4, a5, a6))
    elif surfbcType == 2:
        fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time attc Hsw Ta Pa hr cc cw ua va')
        fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + LakeName + ' Met Data ' + str(year) + '\n')
        fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + '\n')
        fid.write('%s' % '   npts = ' + str(r) + '\n')
        eta = kw['eta']
        Hswn = kw['Hswn']
        Ta = kw['Ta']
        Pa = kw['Pa']
        RH = kw['RH'] / 100
        Cl = kw['Cl']
        cw = kw['cw']
        u = kw['u']
        v = kw['v']
        TimeSim = kw['t_Sim']
        for i in range(0, r):
            a0 = (days[i] - daystart) * 24
            a1 = eta[i]  # light attenuation coefficient
            a2 = Hswn[i]  # Penetrative component of heat flux (albedo already taken into account)
            a3 = Ta[i]  # Air temperature
            a4 = Pa[i]  # Atmospheric pressure
            a5 = RH[i]  # relative humidty (fraction)
            a6 = Cl[i]  # cloud cover (fraction)
            a7 = cw[i]  # **** Wind drag coefficient
            a8 = u[i]  # **** Wind speed in the EW direction
            a9 = v[i]  # **** Wind speed in the NS direction
            if a4 >= 100000:
                format = '%10.4f %10.4f %10.4f %10.4f %10.3f %10.4f %10.4f %10.4f %10.4f %10.4f \n'
            else:
                format = '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n'
            fid.write(format % (a0, a1, a2, a3, a4, a5, a6, a7, a8, a9))
        if show:
            fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
            fig2, (ax5, ax6, ax7, ax8) = plt.subplots(nrows=4, ncols=1)
            fig1.set_size_inches(6, 8)
            fig2.set_size_inches(6, 8)
            ax1.plot(TimeSim, eta)
            ax2.plot(TimeSim, Hswn)
            ax3.plot(TimeSim, Cl)
            ax4.plot(TimeSim, (u ** 2 + v ** 2) ** 0.5)
            ax5.plot(TimeSim, Ta)
            ax6.plot(TimeSim, Pa)
            ax7.plot(TimeSim, RH)
            ax8.plot(TimeSim, cw)

            ax1.set_ylabel(r'$eta$')
            plt.tight_layout()
            ax2.set_ylabel(r'$Hswn\ [Wm^{-2}]$')
            ax3.set_ylabel(r'$Cloud Cover$')
            ax4.set_ylabel(r'$Wspd\ [ms^{-1}]$')
            ax5.set_ylabel(r'$Ta\ [^{\circ}C]$')
            ax6.set_ylabel(r'$Atm\ P\ [Pa]$')
            ax7.set_ylabel(r'$RH$')
            ax8.set_ylabel(r'$Wind\ Drag$')
            ax8.set_xlabel('day of year')
            plt.tight_layout()
            plt.show()
        else:
            print('No plot')
    elif surfbcType == 3:
        fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time attc Hsw Ta Pa hr Hlw cw ua va')
        fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + LakeName + ' Met Data ' + str(year) + '\n')
        fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + '\n')
        fid.write('%s' % '   npts = ' + str(r) + '\n')
        eta = kw['eta']
        Hswn = kw['Hswn']
        Ta = kw['Ta']
        Pa = kw['Pa']
        RH = kw['RH'] / 100
        Hlwin = kw['Hlwin']
        cw = kw['cw']
        u = kw['u']
        v = kw['v']
        TimeSim = kw['t_Sim']
        for i in range(0, r):
            a0 = (days[i] - daystart) * 24
            a1 = eta[i]  # light attenuation coefficient
            a2 = Hswn[i]  # Penetrative component of heat flux (albedo already taken into account)
            a3 = Ta[i]  # Air temperature
            a4 = Pa[i]  # Atmospheric pressure
            a5 = RH[i]  # relative humidty (fraction)
            a6 = Hlwin[i]  # Longwave radiation in
            a7 = cw[i]  # **** Wind drag coefficient
            a8 = u[i]  # **** Wind speed in the EW direction
            a9 = v[i]  # **** Wind speed in the NS direction
            if a4 >= 100000:
                format = '%10.4f %10.4f %10.4f %10.4f %10.3f %10.4f %10.4f %10.4f %10.4f %10.4f \n'
            else:
                format = '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n'
            fid.write(format % (a0, a1, a2, a3, a4, a5, a6, a7, a8, a9))

        if show:
            fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
            fig2, (ax5, ax6, ax7, ax8) = plt.subplots(nrows=4, ncols=1)
            fig1.set_size_inches(6, 8)
            fig2.set_size_inches(6, 8)
            ax1.plot(TimeSim, eta)
            ax2.plot(TimeSim, Hswn)
            ax3.plot(TimeSim, Hlwin)
            ax4.plot(TimeSim, (u ** 2 + v ** 2) ** 0.5)
            ax5.plot(TimeSim, Ta)
            ax6.plot(TimeSim, Pa)
            ax7.plot(TimeSim, RH)
            ax8.plot(TimeSim, cw)

            ax1.set_ylabel(r'$eta$')
            plt.tight_layout()
            ax2.set_ylabel(r'$Hswn\ [Wm^{-2}]$')
            ax3.set_ylabel(r'$Hlwin\ [Wm^{-2}]$')
            ax4.set_ylabel(r'$Wspd\ [ms^{-1}]$')
            ax5.set_ylabel(r'$Ta\ [^{\circ}C]$')
            ax6.set_ylabel(r'$Atm\ P\ [Pa]$')
            ax7.set_ylabel(r'$RH$')
            ax8.set_ylabel(r'$Wind\ Drag$')
            ax8.set_xlabel('day of year')
            plt.tight_layout()
            plt.show()
        else:
            print('No plot')
    elif surfbcType == 10:
        n_stat = kw['n_stat']
        fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time attc Hsw Ta Pa hr cc cw ua va')
        fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + LakeName + ' Met Data ' + str(year) + '\n')
        fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + '\n')
        fid.write('%s' % ' No. Stats ' + str(n_stat) + '\n')
        fid.write('%s ' % ' Grid Locs')
        eta = kw['eta']
        Hswn = kw['Hswn']
        Ta = kw['Ta']
        Pa = kw['Pa']
        RH = kw['RH'] / 100
        Cl = kw['Cl']
        u = kw['u']
        v = kw['v']
        TimeSim = kw['t_Sim']
        imet = kw['imet']
        jmet = kw['jmet']

        for j in range(0, n_stat):
            format = '%10.2f %10.2f '
            fid.write(format % (imet[j], jmet[j]))
        fid.write('%s' % '\n')
        fid.write('%s' % '   npts = ' + str(r) + '\n')

        for i in range(0, r):
            a0 = (days[i] - daystart) * 24
            a1 = eta[i]  # light attenuation coefficient
            a2 = Pa[i]  # Atmospheric pressure
            format = '%10.4f %10.4f %10.3f '
            fid.write(format % (a0, a1, a2))
            for sta in range(0, n_stat):
                format = '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f '
                a3 = Hswn[i]  # Penetrative component of heat flux (albedo already taken into account)
                a4 = Ta[i, sta]  # Air temperature
                a5 = RH[i, sta]  # relative humidty (fraction)
                a6 = Cl[i, sta]  # **** Cloud cover
                a7 = u[i, sta]  # **** Wind speed in the EW direction
                a8 = v[i, sta]  # **** Wind speed in the NS direction
                fid.write(format % (a3, a4, a5, a6, a7, a8))
            fid.write('%s' % '\n')

    elif surfbcType == 11:
        n_stat = kw['n_stat']
        fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time attc Hsw Ta Pa hr cc cw ua va')
        fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + LakeName + ' Met Data ' + str(year) + '\n')
        fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + '\n')
        fid.write('%s' % ' No. Stats ' + str(n_stat) + '\n')
        fid.write('%s ' % ' Grid Locs')
        eta = kw['eta']
        Hswn = kw['Hswn']
        Ta = kw['Ta']
        Pa = kw['Pa']
        RH = kw['RH'] / 100
        Hlwin = kw['Hlwin']
        u = kw['u']
        v = kw['v']
        TimeSim = kw['t_Sim']
        imet = kw['imet']
        jmet = kw['jmet']

        for j in range(0, n_stat):
            format = '%10.2f %10.2f '
            fid.write(format % (imet[j], jmet[j]))
        fid.write('%s' % '\n')
        fid.write('%s' % '   npts = ' + str(r) + '\n')

        for i in range(0, r):
            a0 = (days[i] - daystart) * 24
            a1 = eta[i]  # light attenuation coefficient
            a2 = Pa[i]  # Atmospheric pressure
            format = '%10.4f %10.4f %10.3f '
            fid.write(format % (a0, a1, a2))
            for sta in range(0, n_stat):
                format = '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f '
                a3 = Hswn[i, sta]  # Penetrative component of heat flux (albedo already taken into account)
                a4 = Ta[i, sta]  # Air temperature
                a5 = RH[i, sta]  # relative humidty (fraction)
                a6 = Hlwin[i, sta]  # **** Cloud cover
                a7 = u[i, sta]  # **** Wind speed in the EW direction
                a8 = v[i, sta]  # **** Wind speed in the NS direction
                fid.write(format % (a3, a4, a5, a6, a7, a8))
            fid.write('%s' % '\n')

    elif surfbcType == 20:
        fid.write('%s\n' % '   Time in   // Data format is (10X,G11.2,...) Time cw ua va')
        fid.write('%s' % '   ' + str(dt) + '-min    // SOURCE = ' + LakeName + ' Met Data \n')
        fid.write('%s' % ' intervals  (Note : file prepared on ' + str(Dt.date.today()) + '\n')
        fid.write('%s' % '   npts = ' + str(r) + '\n')

        cw = kw['cw']
        u = kw['u']
        v = kw['v']
        for i in range(0, r):
            a0 = (days[i] - daystart) * 24
            a1 = cw[i]  # **** Wind drag coefficient
            a2 = u[i]  # **** Wind speed in the EW direction
            a3 = v[i]  # **** Wind speed in the NS direction
            format = '%10.4f %10.4f %10.4f %10.4f \n'
            fid.write(format % (a0, a1, a2, a3))

    fid.close()
    return


def HeatBudget(HeatBudgetMethod, eta, Hswn, Hlwin, Hlwout, Ta, Pa, RH, Cl, cw, u, v, WaTemp, cChapra, esMethod):
    if HeatBudgetMethod == 'Chapra1995':
        rho0 = 997
        Lv = 2.5e6
        CbPa_P = 0.61 * cChapra
        wspd = (u ** 2 + v ** 2) ** 0.5
        # Vapor Pressure
        if esMethod == 1:
            es = 6.11 * np.exp(17.3 * Ta / (Ta + 237.3))
            esw = 6.11 * np.exp(17.3 * WaTemp / (WaTemp + 237.3))
        elif esMethod == 2:
            es = 6.11 * np.exp(7.5 * Ta / (Ta + 237.3))
            esw = 6.11 * np.exp(7.5 * WaTemp / (WaTemp + 237.3))
        elif esMethod == 3:
            es = 10 ** (9.286 - (2322.38 / (Ta + 273.15)))
            esw = 10 ** (9.286 - (2322.38 / (WaTemp + 273.15)))
        ea = es * RH
        # Longwave radiation
        Hlwn = Hlwin - Hlwout
        # Latent Heat Flux (Negative as it exits)
        fwind = 1.02e-9 * wspd
        Hl = -rho0 * Lv * fwind * (esw - ea)
        # Sensible heat flux (negative due to consideration of difference between water temp and air temp)
        Hs = -rho0 * Lv * fwind * CbPa_P * (WaTemp - Ta)
        # Net Heat Flux
        Hn = Hswn + Hlwn + Hs + Hl
    elif HeatBudgetMethod == 'AirSea':
        print('UNDER DEVELOPMENT, THIS FUNCTION DOES NOT WORK')
        print('The file has not been created')
        exit()
    elif HeatBudgetMethod == 'TERC':
        print('UNDER DEVELOPMENT, THIS FUNCTION DOES NOT WORK')
        print('The file has not been created')
        exit()
    return Hswn, Hlwn, Hl, Hs, Hn
