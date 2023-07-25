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
md = md_struct.md;

% Interpolate bed elevation, ice thickness onto mesh
% given by md_struct.md.mesh;
bed_interpolant = scatteredInterpolant(bed(:, 1), bed(:, 2), bed(:, 3), 'nearest');
sfc_interpolant = scatteredInterpolant(sfc(:, 1), sfc(:, 2), sfc(:, 3), 'nearest');

bed_elev = bed_interpolant(md.mesh.x, md.mesh.y);
sfc_elev = sfc_interpolant(md.mesh.x, md.mesh.y);

% Set ice thickness to at least 30 m everywhere
bed_elev = min(bed_elev, sfc_elev - 30);

% Save elevations to text files
save('bed_mesh.txt', 'bed_elev', '-ascii')
save('sfc_mesh.txt', 'sfc_elev', '-ascii')

% BONUS: look at hydraulic potential on boundary nodes
rhow = 1000;
rhoi = 910;
g = 9.81;
phi_bed = rhow*g*bed_elev;
p_ice = rhoi*g*(sfc_elev - bed_elev);
phi = phi_bed + 0.8*p_ice;
plot(phi(md.mesh.vertexonboundary==1))

% Plot potential on boundary, label nodes
figure
bndry = md.mesh.vertexonboundary==1;
bndry_x = md.mesh.x(bndry);
bndry_y = md.mesh.y(bndry);
bndry_index = string(1:length(bndry_y));
bndry_phi = phi(bndry);
scatter(bndry_x, bndry_y, 30, bndry_phi)
text(bndry_x, bndry_y, bndry_index)
hold on
[phi_min, phi_min_index] = min(bndry_phi);
plot(bndry_x(phi_min_index), bndry_y(phi_min_index), 'rx')

bndry_indices = find(bndry);
outlet_index = bndry_indices(phi_min_index);
disp('Outlet index:')
disp(outlet_index)

% Plot potential just on terminus lobe - requires manual intervention
terminus_nodes = [181, 327, 182, 204, 183, 184, 217, 1, 218, 2, 250,...
    259, 3, 260, 4, 261, 5, 279, 6];
figure
plot(bndry_phi(terminus_nodes))
text(1:length(terminus_nodes), bndry_phi(terminus_nodes), string(terminus_nodes))
grid on

disp('Multiple outlet indices:')
outlet_indices = bndry_indices([5, 204]);
disp(outlet_indices)
