function para = get_para_valley(meshnr, bed_para, c_t, source_scenario_nr);
%  para = get_para_seasonal_valley(meshnr, bed_para, source=);

if ~exist('source_scenario_nr', 'var') || isempty(source_scenario_nr)
    source_scenario_nr = 6;
end

if ~exist('c_t', 'var') || isempty(c_t)
    c_t = -7.5000e-08;
end

para = get_default_para();
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;

%% Model description
% a descriptive string for the model-run
pm.model_run_descript = ['Hydro intercomparsion: valley'];

%% some directories (not added to matlab path
% set all relative to either pm.dir.model_runs or pm.dir.glads
[~, dir_of_this_mfile] = get_mfile_name();
pm.dir.model_runs = dir_of_this_mfile;
pm.dir.data = [pm.dir.model_runs, 'data', '/'];
pm.save_filename_root = '';
pm.dir.model_save = [dir_of_this_mfile, 'output/'];

% mesh file
pm.file.mesh = [pm.dir.data, 'mesh/valley_meshes.mat'];

pm.path.topo_fns = pm.dir.data;
pm.path.topo_fns2 = [pm.dir.data, '../../../input_functions/topography/'];
add_paths(pm.path);



% TIME
%%%%%%
pt.start = 0;        % start time
pt.end   = 1000*pp.day;  % end time
pt.out_t = pt.start:10*pp.day:pt.end;

%pt.out_t = pt.start:pp.day:10*pp.day;

% MESH
%%%%%%
dmesh = load(pm.file.mesh);
% get desired mesh out (thanks to matlab's syntax this is butt-ugly)
fns = fieldnames(dmesh);
dmesh = getfield(dmesh,fns{1});
dmesh = dmesh{meshnr};


%  PYHSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%
pp.e_v = 1e-3;

pp.cond_s = 5e-3;  %  conductivity of sheet
pp.cond_c = 1e-1;  %  conductivity of channels

pp.c_t_c = c_t;

%  NUMERICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{5};

%% BC
% Dirichlet BC for phi: applied at nodes where bmark is odd
pin.bc_dirichlet_phi = make_anon_fn('@(xy, time, bmark, phi_0, phi_m) double(0*phi_m)');
% Flux BC for phi and h_w: i.e. set phi or h_w such that this flux is
% given. Applied at edges where bmark_edge is even
% zero flux:
pin.bc_flux = make_anon_fn('@(xy, time, bmark_edge) double(zeros(sum(~logical(mod(bmark_edge,2)) & bmark_edge>0),1))');

%% IC
% initial sheet thickness
pin.ic_h = make_anon_fn('@(xy, time) double(0.05 + 0*xy(:,1))');
% initial channel cross sectional area
pin.ic_S =  make_anon_fn('@(xy, time) double(0.0 + 0*xy(:,1))');

%% Geometry and derived (constant) function
% these are all functions of xy and time only
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation(xy,time,bed_para))', bed_para);
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness(xy,time,bed_para))', bed_para);

%% Source functions
source = 2*steady_runoff(source_scenario_nr); % default == twice scenario A6
pin.source_term_s = make_anon_fn('@(xy,time) double(xy(:,1)*0 + source)', source);
tmpvec = zeros(dmesh.tri.n_nodes,1);
pin.source_term_c = make_anon_fn('@(time) double(tmpvec*0)', tmpvec);
%

%% Storage terms: no storage
% englacial storage (set to zero)
pin.sigma = make_anon_fn('@(xy, time) double(xy(:,1).*0 + pp.e_v)', pp);  % englacial void ratio

% SCALING PARAMETERS
%%%%%%%%%%%%%%%%%%%%
% now that we have the dmesh we can get the length scale
ps.x = 0.5*(abs(diff(dmesh.x_extent)) + abs(diff(dmesh.y_extent))); % length (m)
ps.x_offset = dmesh.x_extent(1);
ps.y_offset = dmesh.y_extent(1);

% the rest of the scales are calculated:
ps = set_default_scales(ps, pp, dmesh);
% if any scales are to be changed from what set_default_scales sets
% them to do this here:

% SCALE EVERYTHING
%%%%%%%%%%%%%%%%%%
% scaling of the rest
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);

% END
%%%%%
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
