% Matlab Script to create the Valley mesh for tests E and F of SHMIP
% O. Gagliardini, Septembre 2016

% create valley.geo
clear;


lc_out=40.0;

A=dlmread('contour.dat');

fid1=fopen('valley.geo','w');
fprintf(fid1,'Mesh.Algorithm=5; \n');

As=size(A,1);

np=0;
for ii=1:As
    np=np+1;
    fprintf(fid1,'Point(%g)={%14.7e,%14.7e,0.0,%g}; \n',np,A(ii,1),A(ii,2),lc_out);
end

fprintf(fid1,'Spline(1)={');
for ii=1:As
  fprintf(fid1,'%g,',ii);
end
fprintf(fid1,'%g}; \n',1);

fprintf(fid1,'Line Loop(2)={1}; \n');
fprintf(fid1,'Plane Surface(3) = {2}; \n');
fprintf(fid1,'Physical Line(4) = {1}; \n');
fprintf(fid1,'Physical Surface(5) = {3}; \n');

fclose(fid1)

% create teterousse.msh using gmsh
system "gmsh teterousse.geo -1 -2"
% convert teterousse.gmsh in an Elmer type mesh
system "ElmerGrid 14 2 teterousse.msh -autoclean"

