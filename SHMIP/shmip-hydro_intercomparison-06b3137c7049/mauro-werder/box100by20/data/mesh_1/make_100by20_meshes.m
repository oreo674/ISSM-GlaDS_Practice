function box100by20 = make_100by20_meshes()
% box100by20 = make_100by20_meshes()
%
% Makes the box meshes of different density.
%
% Dirichlet on front, no flux the rest.

% the path of this mfile (used for paths below)
mfiledir = [fileparts(mfilename('fullpath')), '/'];

% make the box and write it to a file which the meshing library (Triangle) knows
boundary_xy = [0, 0;
               100e3, 0;
               100e3, 20e3;
               0, 20e3];

% boundary marks>0 on edge:
bmark = [1;2;2;1];         % just a mark which is given to the nodes on the boundary
bmark_edge = [2;2;2;1];  % just a mark which is given to the edges on the boundary

maxareas = 250000000*[0.1, 0.01, 5e-3, 2.5e-3, 1e-3, 0.5e-3]; % area of triangle , 1e-4

% cell array holding all the meshes
box100by20 = {};
for ii=1:length(maxareas)
    box100by20{ii} = make_mesh(boundary_xy, bmark, bmark_edge, maxareas(ii));
end

save([mfiledir, '/box100by20.mat'], 'box100by20');

