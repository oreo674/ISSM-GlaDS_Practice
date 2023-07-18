% Interpolate bed and surface elevation onto numerical mesh

% File structure
mesh_fname = '../KWsmooth.mesh.mat';
bed_fname = 'bed_DEM.dat';
sfc_fname = 'sfc_DEM.dat';
% bed_fname, sfc_fname files have three columns:
% x, y, elevation
% where elevation is bed or surface elevation

% Read from files
bed = importdata(bed_fname);
sfc = importdata(sfc_fname);
md_struct = load(mesh_fname);

% TBD...
%
% Interpolate bed elevation, ice thickness onto mesh
% given by md_struct.md.mesh;
%
% Save these fields in your format of choice (*.mat, *.csv, *.txt, ...)
% Read these in ../run_KWsmooth.m