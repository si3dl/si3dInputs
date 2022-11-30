% si3dInputs.m
% This script serves to create the files needed for the SI3D model runs. The code is based on previous matlab versions created by Alicia Cortes and Francisco Rueda.
% Functions that are present within this script are:


% 1. bathy4si3d
%     This function writes the bathymetry file 'h' for si3d simulations. 
%     if the basin is a real lake use the functions as: bathy4si3d(lakename,fpath,savepath,bathyfile,shoreline,dx);

% 2. initcond4si3d2
%   This function creates the initial condition file (si3d_init.txt) nad the
%   layer file (s3id_layer.txt) based on input tab deliminated profiles
%   NB - input profile must between specified as tab deliminteated file
%   Four possible inputs
%   NB - input profile must specify temps to bottom of lake

% initcond4si3d2(LakeName,fpath,savepath,SimStartDate,spacing_method,TempProf,dz,dzmin,dzmax);
% LakeName = lake Name
% fpath - path to initial condition profile
% savepath - path to where inputs will be saved
% SimStartDate = startdate speciefied as DATETIME
% spacing_mod: constant (constant dz), expo (expoentialy increasing dx size
% from surface
% TempProf: variable (variable temp based on input profile) / constant
% (taken as average of input profile)
% dz = constant dz [m]
% dzmin = minimum dz for exponential spacing
% dzmax max dz for exponential spacing
clear all; close all; clc


LakeName = 'Llanquihue';
pathname = '/Users/micahswann/Documents/GitHub/si3dInputs/_matllibrary_/Example'; % path to init cond file files
savepath = '/Users/micahswann/Dropbox/03-CHILELAGOS_TAHOE/32- Si3D/2_Llanquihue/3_Inputs/init_cond'; % location to save files
TempProf = 'variable'; % variable or uniform
start_temp = 14.5; % starting initial temp if TempProf is UNIFORM
spacing_method = 'expo'; % %=
dz = 2; % vertical grid size for CONSTANT layering
dzmin = 0.5; % minimum dz for exponential spaced layers
dzmax = 10; % maximum dz for exponetially spaced layers
fpath = '/Users/micahswann/Desktop/profile_example.txt';
SimStartDate = datetime(2021,01,01);



initcond4si3d(LakeName,fpath,savepath,SimStartDate,spacing_method,TempProf,dz,dzmin,dzmax);