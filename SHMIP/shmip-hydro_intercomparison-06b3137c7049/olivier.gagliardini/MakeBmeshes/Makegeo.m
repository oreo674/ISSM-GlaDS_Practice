% Matlab Script to create meshes for suite B tests of SHMIP
% O. Gagliardini, September 2016

% create mesh_XX.geo
clear;

for bb=1:5

    Moulin=dlmread(['./B' num2str(bb) '_M.xy']);
    fid1=fopen(['mesh_B' num2str(bb) '.geo'],'w');

    fprintf(fid1,'L = 100.0e3 ; \n');              
    fprintf(fid1,'W = 20.0e3 ; \n');

    fprintf(fid1,'Mesh.Algorithm = 2 ; \n');

    fprintf(fid1,'lc =  500.0 ; \n');

    fprintf(fid1,'Point(1) = {0.0,0.0,0.0,lc}; \n');
    fprintf(fid1,'Point(2) = { L ,0.0,0.0,lc}; \n');
    fprintf(fid1,'Translate {0.0, W , 0.0} {Duplicata{Point{2} ; } }\n');
    fprintf(fid1,'Translate {0.0, W , 0.0} {Duplicata{Point{1} ; } }\n');


    ms=size(Moulin,1);
    np=4;
    for ii=1:ms
        np=np+1;
        fprintf(fid1,'Point(%g)={%14.7e,%14.7e,0.0,lc}; \n',np,Moulin(ii,1),Moulin(ii,2));
    end

    fprintf(fid1,'Line(1) = {1,2} ; \n');
    fprintf(fid1,'Line(2) = {2,3} ; \n');
    fprintf(fid1,'Line(3) = {4,3} ; \n');
    fprintf(fid1,'Line(4) = {4,1} ; \n');

    fprintf(fid1,'Line Loop(5) = {1,2,-3,4}; \n');
    fprintf(fid1,'Plane Surface(10) = {5}; \n');

    for ii=1:ms
        jj = 4+ii ;
        fprintf(fid1,'Point{%g} In Surface{10}; \n',jj);
    end

    for ii=1:ms
        jj = 4+ii ;
        fprintf(fid1,'Physical Point(%g) = {%g}; \n',ii,jj);
    end

    fprintf(fid1,'Physical Line(%g) = {1}; \n',ms+1);
    fprintf(fid1,'Physical Line(%g) = {2}; \n',ms+2);
    fprintf(fid1,'Physical Line(%g) = {3}; \n',ms+3);
    fprintf(fid1,'Physical Line(%g) = {4}; \n',ms+4);
    fprintf(fid1,'Physical Surface(%g) = {10}; \n',ms+5);

    fclose(fid1)

end

% create .msh using gmsh
%system "gmsh mesh_70moulins.geo -1 -2"
% convert .msh in an Elmer type mesh
%system "ElmerGrid 14 2 mesh_70moulins.msh -autoclean"


