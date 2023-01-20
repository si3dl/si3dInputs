% si3dInputs.m
% This script serves to create the files needed for the SI3D model runs. The code is based on previous matlab versions created by Alicia Cortes and Francisco Rueda.
% Functions that are present within this script are:

% Last Modified Jan 2023 by M. Swann

% 1. bathy4si3d
%     This function writes the bathymetry file 'h' for si3d simulations. 
% example function call - bathy4si3d(LakeName,pathname,savepath,bathyfile,shoreline,dx);
%     Inputs
%           - LakeName = name of lake
%           - pathname = path to where raw bathy files are storede
%           - savepath = path to where si3d_input files will be saved
%           - bathyfile = name of bathy file (csv)
%           - shoreline = name of shoreline file (csv)
%           - dx - horizontal grid size in meters

% 2. initcond4si3d
%   This function creates the initial condition file (si3d_init.txt) and the
%   layer file (s3id_layer.txt) based on input tab deliminated profiles

%   NB - input profile must be specified as tab deliminteated file
%   NB - input profile must specify temps to bottom of lake

% example function call - %initcond4si3d(LakeName,initcondfile,savepath,SimStartDate,spacing_method,TempProf,dz,dzmin,dzmax);
%   Inputs
%           - LakeName = lake Name
%           - fpath - path to initial condition profile
%           - savepath = path to where si3d_input files will be saved
%           - SimStartDate = startdate speciefied as DATETIME (Local Time)
%           - spacing_mod (constant / expo) : constant (constant dz), expo (expoentialy
%               increasing dx from surface from dzmin to dzmax
%           - TempProf (variable / constant): variable (variable temp based on input profile) / constant
%               (taken as average of input profile)
%           - dz = constant dz [m]
%           - dzmin = minimum dz for exponential spacing
%           - dzmax max dz for exponential spacing

% 3. surfbc4si3d
%   This function creates the surfbace boundary condition file (surbc.txt)
%   NB  -  Input met file must be atab delinated text file specifying met variables in each column 
%   in the following order with the first row as headings
%       Column 1 - dt - datetime in the format YYYY-MM-DD HH:mm:ss. Date should
%           be in local time for cloud cover calculations
%       Column 2 - attc [dimensionless] - light attenuation coefficient (1.7/secchi depth) can
%           be variable
%       Column 3 - Hsw [w m^-2] - incoming shorwtwave radiation 
%       Column 4 - Ta [C] - surface air temperature
%       Column 5 - Pa [kPa] - atmospheric pressure at lake surface
%       Column 6 - hr [%] - relative humidity as percent (e.g. 100 no 1 for   100%)
%       Column 7 - cw [dimensionless] - wind drag 9can be variable
%       Column 8 - WS [m/s] - wind speed at 10 m
%       Column 9 - WDir [degrees] - wind direction (0 degrees is north)

% example function call - surfbc4si3d(LakeName,metfile,savepath,surfbcType,SimStartDate,duration,dt,lat,lon,lsm,hemisphere,spinup);
%   Inputs
%           - LakeName = lake Name
%           - metfile = path to met input file
%           - savepath = path to where si3d_input files will be saved
%           - surfbcType = 2 (Run-Time I) / 3 (Run-Time II) currently
%           suppoorted
%           - SimStartDate = startdate speciefied as DATETIME (Local Time)
%           - duration = duration of run in seconds
%           - dt - timestep of met data output
%           - lat - latitude of study site [degrees]
%           - lon - longitude of study site [degrees]
%           - lsm - local stanard meridian
%           - hemisphere - hemispehre of study site ('N','S') need for
%           correct astronomic calculations
%           - spinup time - model spin time in seconds, winds will be set
%           to 0.001 m/s during this period. 3 days is recommended

clear all; close all; clc; tic

%% Example User Inputs

% path/ file names
LakeName = 'Llanquihue';
pathname = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example'; % path to init cond file files
savepath = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/si3d_inputs'; % location to save input files
initcondfile = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example/profile_example.txt';
bathyfile = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example/Llanquihue_Bathy_50m.csv';
shorelinefile = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example/Llanquihue_Shoreline.csv';
metfile = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example/surfbc_example.txt';
% variable inputs

% bathy
dx = 400; % grid size for bathy h file

% init cond
TempProf = 'variable'; % variable or uniform
start_temp = 14.5; % starting initial temp if TempProf is UNIFORM
spacing_method = 'constant'; % %=
dz = 2; % vertical grid size for CONSTANT layering
dzmin = 0.5; % minimum dz for exponential spaced layers
dzmax = 10; % maximum dz for exponetially spaced layers
SimStartDate = datetime(2022,05,01);
    
% surfbc
surfbcType = 3; % Type 2 = runtime I ; Type 3 = runtime II
duration = 5*24*3600; % duration in seconds
dt = 3600; % intended timestep of met data [s]
lat = 41.1; % latitude in degrees
lon = 72.8; % longitude [degrees]
lsm = 50; % local standard meridian
hemisphere = 'S'; % N/S; % if S the decliniation of the sun in fliped in cloudiness calcs
spinup = 3*3600*24; % spinup time in seconds. Winds will be set to 0.001 m/s

%% Example function calls

% bathymetry
bathy4si3d(LakeName,pathname,savepath,bathyfile,shorelinefile,dx);
toc
disp('bathy input file created')

% init condition
initcond4si3d(LakeName,initcondfile,savepath,SimStartDate,spacing_method,TempProf,dz,dzmin,dzmax);
toc
disp('init condition files created')

% surfbc
surfbc4si3d(LakeName,metfile,savepath,surfbcType,SimStartDate,duration,dt,lat,lon,lsm,hemisphere,spinup);
toc
disp('surfbc file created')