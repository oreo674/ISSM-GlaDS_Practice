function para = get_para_box100by20(ch, meshnr, melt_scenario, bedslope, press_melt);
% para = get_para_box100by20(ch, meshnr, melt_scenario, bedslope);

addpath('..')
para = get_default_para();
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;
rmpath('..')

if ~exist('bedslope', 'var') || isempty(bedslope)
    bedslope = 0;
end
if ~exist('press_melt', 'var') || isempty(press_melt)
    press_melt = true;
end

%% Model description
% a descriptive string for the model-run
pm.model_run_descript = ['Hydro intercomparsion: box100by20'];

%% some directories (not added to matlab path
% set all relative to either pm.dir.model_runs or pm.dir.glads
[~, dir_of_this_mfile] = get_mfile_name();
pm.dir.model_runs = dir_of_this_mfile;
pm.dir.data = [pm.dir.model_runs, 'data', '/'];
pm.save_filename_root = '';
% directory to save model output
pm.dir.model_save = 'outputmat/'; % to be set

% mesh file
pm.file.mesh = [pm.dir.data, 'mesh_1/box100by20.mat'];

pm.path.topo_fns = pm.dir.data;
add_paths(pm.path);


% TIME
%%%%%%
pt.start = 0;        % start time
pt.end   = 8000*pp.day;  % end time
pt.out_t = pt.start:100*pp.day:pt.end;

% MESH
%%%%%%
dmesh = load(pm.file.mesh);
% get desired mesh out (thanks to matlab's syntax this is butt-ugly)
fns = fieldnames(dmesh);
dmesh = getfield(dmesh,fns{1});
dmesh = dmesh{meshnr};

%  PYHSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%
if ~ch
    pp.l_c = 0;       % width of sheet which contributes to channel melt (m) (edges)
end
pp.e_v = 0*1e-4;

pp.cond_s = 5e-3;  %  conductivity of sheet
pp.cond_c = 1e-1;  %  conductivity of channels

if press_melt
    pp.c_t_c = -7.5e-8;
else
    pp.c_t_c = 0;
end

%  NUMERICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{2};

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
% pin.bed_elevation = make_anon_fn('@(xy, time) double(xy(:,1)*0)');
% $$$ thick_para = 6;
% $$$ thick_off = 5e3;
% $$$ thick_nz = 1;
%pin.ice_thickness = make_anon_fn('@(xy, time) double(thick_para*sqrt(xy(:,1)+thick_off) - thick_para*sqrt(thick_off) + thick_nz )', thick_para, thick_off, thick_nz );

pin.bed_elevation = make_anon_fn('@(xy, time) double(sqrttopo_bed(xy(:,1), xy(:,2), bedslope))', bedslope );
pin.ice_thickness = make_anon_fn('@(xy, time) double(sqrttopo_ice(xy(:,1), xy(:,2), bedslope))', bedslope );

%% Source functions
source = steady_runoff(melt_scenario);
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
