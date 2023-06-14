addpath('..')

%% mesh sizes
domainarea = 5.6853e+06;
areabox100by20 = 100e3*20e3;
maxareas =   250000000 * domainarea/areabox100by20 * [0.1, 0.01, 5e-3, 2.5e-3, 1e-3, 0.5e-3]; % area of triangle
ns = 6000./sqrt(maxareas)*2;

valley_meshes = {};
for ii = 1:length(maxareas)
    % make the outlines with the right amount of points
    [x,y] = glacier_outline(ns(ii));

    xy_outline = [x',y'];
    % set bmark to 2
    cborder_bmark = xy_outline(:,1)*0 + 2;
    cborder_bmark_edge = cborder_bmark;
    % except at outlet
    cborder_bmark(1) = 1;

    % make the meshes
    valley_meshes{ii} = make_mesh(xy_outline, cborder_bmark, cborder_bmark_edge, maxareas(ii), ...
                                  [],[],[],true);
end

save('valley_meshes.mat', 'valley_meshes');

rmpath('..')
