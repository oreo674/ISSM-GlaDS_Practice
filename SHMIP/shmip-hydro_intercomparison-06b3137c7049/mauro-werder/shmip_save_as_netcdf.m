function shmip_save_as_netcdf(os, filename, timeinds, force)
% save_model_run_netcdf(out_struct, filename, timeinds, force)
%
% on unstructured mesh
%
% force -- if true overwrite existing file (default false)

if ~exist('force', 'var') || isempty(force)
    force = false;
end
if force
    if exist(filename, 'file')
        delete(filename)
    end
end

if ~exist('timeinds', 'var') || isempty(timeinds)
    timeinds = 1:size(os.h_sheets, 2);
elseif strcmp(timeinds, 'end')
    timeinds = size(os.h_sheets, 2);
end

para = os.para;
ps = para.scale;
pp = para.physical;
tri = para.dmesh.tri;
nt = length(timeinds);
time = para.time.out_t(timeinds) - para.time.out_t(timeinds(1));

cell_midpoints = 1/3*(tri.nodes(tri.connect(:,1),:) + tri.nodes(tri.connect(:,2),:) + tri.nodes(tri.connect(:,3),:));

h_sheet = os.h_sheets(:,timeinds)*ps.h;
phi = os.phis(:,timeinds)*ps.phi;
sigma = para.input.sigma(para.dmesh.tri.nodes, time);
p = bsxfun(@minus, phi, os.fields.nodes.phi_m*ps.phi);
h_store = bsxfun(@times, p/pp.rho_w/pp.g_grav,sigma);

N_eff = pp_get_N(os.phis, os.fields.nodes.phi_0, tri)*ps.phi;
N = N_eff(:,timeinds);
S_channel = os.S_channels(:,timeinds)*ps.S;

ice_thick = os.fields.nodes.ice_thick*ps.z;
bed_ele = os.fields.nodes.bed_ele*ps.z;

% discharge
[dis_c, vel_c] = pp_discharge_edge_scaled(os.S_channels, os.phis, os.para.mesh.tri, os.para.scaled_phys);
Q = dis_c(:,timeinds)*ps.Q;
[dis, grad_phi_el, elements, connect_el, vel] = pp_discharge_element(os.h_sheets, os.phis, os.para.mesh.tri, os.para.scaled_phys);
q = dis(:,:,timeinds)*ps.q;
q = squeeze(sqrt(q(:,1,:).^2+q(:,2,:).^2));

% NETCDF dimensions
index1 = tri.n_nodes;
index2 = tri.n_elements;
index_ch = tri.n_edges;
dim = 2;
dims = { 'time', nt; % / time dimension
         'dim', dim;% / spatial dimensions
         'n_nodes_ch' 2 ; % nodes per channel
         'index1', index1; % / index into coords1
         'index2', index2; % / index into coords2
         'index_ch', index_ch; % / channel index
         };

mySchema.Name   = '/';
mySchema.Format = 'classic';
for i = 1:size(dims,1)
    n = dims{i,1};
    v = dims{i,2};
    mySchema.Dimensions(i).Name   = n;
    mySchema.Dimensions(i).Length = v;
end
ncwriteschema(filename, mySchema);

%% create dimensional variable time
nccreate(filename, 'time', ...
          'Dimensions',{'time' nt});
ncwrite(filename, 'time', time);
ncwriteatt(filename, '/time', 'units', 's')
ncwriteatt(filename, '/time', 'long_name', 'time')


%% dependent variables
% Mesh:
nccreate(filename, 'coords1', ...
          'Dimensions',{'index1' index1, 'dim' dim});
ncwrite(filename, 'coords1', tri.nodes);
ncwriteatt(filename, '/coords1', 'units', 'm');
ncwriteatt(filename, '/coords1', 'long_name', 'node coordinates');

nccreate(filename, 'coords2', ...
          'Dimensions',{'index2' index2, 'dim' dim});
ncwrite(filename, 'coords2', cell_midpoints);
ncwriteatt(filename, '/coords2', 'units', 'm');
ncwriteatt(filename, '/coords2', 'long_name', 'cell midpoint coordinates');

nccreate(filename, 'coords_ch', ...
          'Dimensions',{'index_ch' index_ch, 'dim' dim});
ncwrite(filename, 'coords_ch', tri.edge_midpoints);
ncwriteatt(filename, '/coords_ch', 'units', 'm');
ncwriteatt(filename, '/coords_ch', 'long_name', 'channel midpoint coordinates');

nccreate(filename, 'connect_ch', ...
          'Dimensions',{'index_ch' index_ch, 'n_nodes_ch' 2});
ncwrite(filename, 'connect_ch', tri.connect_edge-1);
ncwriteatt(filename, '/connect_ch', 'units', '')
ncwriteatt(filename, '/connect_ch', 'long_name', 'channel connectivity')



% constant vars
nccreate(filename, 'B', ...
          'Dimensions',{'index1' index1});
ncwrite(filename, 'B',bed_ele);
ncwriteatt(filename, '/B', 'units', 'm')
ncwriteatt(filename, '/B', 'long_name', 'bed elevation')

nccreate(filename, 'H', ...
          'Dimensions',{'index1' index1});
ncwrite(filename, 'H', ice_thick);
ncwriteatt(filename, '/H', 'units', 'm')
ncwriteatt(filename, '/H', 'long_name', 'ice thickness')


% time dependent vars
nccreate(filename, 'N', ...
          'Dimensions',{'index1' index1, 'time' nt});
ncwrite(filename, 'N', N);
ncwriteatt(filename, '/N', 'units', 'Pa')
ncwriteatt(filename, '/N', 'long_name', 'effective pressure')

nccreate(filename, 'h', ...
          'Dimensions',{'index1' index1, 'time' nt});
ncwrite(filename, 'h', h_sheet);
ncwriteatt(filename, '/h', 'units', 'm')
ncwriteatt(filename, '/h', 'long_name', 'water sheet thickness')

nccreate(filename, 'hstore', ...
          'Dimensions',{'index1' index1, 'time' nt});
ncwrite(filename, 'hstore', h_store);
ncwriteatt(filename, '/hstore', 'units', 'm')
ncwriteatt(filename, '/hstore', 'long_name', 'stored water effective layer thickness')

nccreate(filename, 'q', ...
          'Dimensions',{'index2' index2, 'time' nt});
ncwrite(filename, 'q', q);
ncwriteatt(filename, '/q', 'units', 'm^2/s')
ncwriteatt(filename, '/q', 'long_name', 'water sheet discharge')

nccreate(filename, 'S', ...
          'Dimensions',{'index_ch' index_ch, 'time' nt});
ncwrite(filename, 'S', S_channel);
ncwriteatt(filename, '/S', 'units', 'm^2')
ncwriteatt(filename, '/S', 'long_name', 'channel cross-sectional area')

nccreate(filename, 'Q', ...
          'Dimensions',{'index_ch' index_ch, 'time' nt});
ncwrite(filename, 'Q', Q);
ncwriteatt(filename, '/Q', 'units', 'm^3/s')
ncwriteatt(filename, '/Q', 'long_name', 'channel discharge')




%% global attributes
[~,run_name] = fileparts(filename);
ncwriteatt(filename, '/', 'title', ['werder_',run_name])
ncwriteatt(filename, '/', 'meshtype', 'unstructured')
ncwriteatt(filename, '/', 'channels_on_edges', 'yes')
ncwriteatt(filename, '/', 'institution', 'Mauro A Werder, ETHZ')
ncwriteatt(filename, '/', 'source', ['GlaDS: ',  para.model.current_glads_git_revision, ' SHMIP: ', strtrim(git('rev-parse --verify HEAD'))])
ncwriteatt(filename, '/', 'references', 'http://shmip.bitbucket.io/')
