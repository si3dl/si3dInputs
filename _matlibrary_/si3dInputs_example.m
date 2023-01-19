% si3dInputs.m
% This script serves to create the files needed for the SI3D model runs. The code is based on previous matlab versions created by Alicia Cortes and Francisco Rueda.
% Functions that are present within this script are:

% Last Modified Jan 2023 by M. Swann

% 1. bathy4si3d_git
%     This function writes the bathymetry file 'h' for si3d simulations. 
%     if the basin is a real lake use the functions as: bathy4si3d(LakeName,pathname,savepath,bathyfile,shoreline,dx);
%     Inputs
%           - LakeName = name of lake
%           - pathname = path to where raw bathy files are storede
%           - savepath = path to where si3d_input files will be saved
%           - bathyfile = name of bathy file (csv)
%           - shoreline = name of shoreline file (csv)
%           - dx - horizontal grid size in meters

% 2. initcond4si3d_git
%   This function creates the initial condition file (si3d_init.txt) nad the
%   layer file (s3id_layer.txt) based on input tab deliminated profiles
%   NB - input profile must between specified as tab deliminteated file

%   Four possible input formats:
%   NB - input profile must specify temps to bottom of lake
%   Inputs
%           - LakeName = lake Name
%           - fpath - path to initial condition profile
%           - savepath = path to where si3d_input files will be saved
%           - SimStartDate = startdate speciefied as DATETIME
%           - spacing_mod (constant / expo) : constant (constant dz), expo (expoentialy
%               increasing dx from surface
%           - TempProf (variable / constant): variable (variable temp based on input profile) / constant
%               (taken as average of input profile)
%           - dz = constant dz [m]
%           -  dzmin = minimum dz for exponential spacing
%           - dzmax max dz for exponential spacing

clear all; close all; clc; tic

%% User Inputs

% path/ file names
LakeName = 'Llanquihue';
pathname = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example'; % path to init cond file files
savepath = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/si3d_inputs'; % location to save input files
fpath = '/Users/micahswann/Documents/GitHub/si3dInputs/_matlibrary_/Example/profile_example.txt';
bathyfile = 'Llanquihue_Bathy_50m.csv';
shoreline = 'Llanquihue_Shoreline.csv';

% variable inputs
dx = 50;
TempProf = 'variable'; % variable or uniform
start_temp = 14.5; % starting initial temp if TempProf is UNIFORM
spacing_method = 'expo'; % %=
dz = 2; % vertical grid size for CONSTANT layering
dzmin = 0.5; % minimum dz for exponential spaced layers
dzmax = 10; % maximum dz for exponetially spaced layers
SimStartDate = datetime(2022,05,01);
    

%% running example

% bathymetry
bathy4si3d(LakeName,pathname,savepath,bathyfile,shoreline,dx);
toc
disp('bathy input file created')

% init condition
initcond4si3d(LakeName,fpath,savepath,SimStartDate,spacing_method,TempProf,dz,dzmin,dzmax);

toc
disp('init condition files created')