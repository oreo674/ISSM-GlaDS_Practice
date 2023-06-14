function para = get_def_para();
%  para = get_def_para();
%
% Returns the default parameters to run the hydro-inter-comparison runs.

error('this does not work')

para = get_default_para();
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;

%% Model description
% a descriptive string for the model-run
pm.model_run_descript = ['Hydro intercomparsion'];
% save-file root (empty means no saveing)
pm.save_filename_root = '';
% save the model every so often
pm.backup_time_steps = 1000;
% IC from previous model run?
pm.IC_from_file = 0;

% TIME
%%%%%%
pt.start = 0;        % start time
pt.end   = 5000*pp.day;  % end time
pt.out_t = pt.start:50*pp.day:pt.end;

% MODEL
%%%%%%%

% how much output will be given
pm.verbosity = 2;
pm.plot_verbosity = 0;

% git revision of model runs dir
pm.git_revision_model_runs = strtrim(git('rev-parse --verify HEAD'));

%% some directories (not added to matlab path
% set all relative to either pm.dir.model_runs or pm.dir.glads
[~, dir_of_this_mfile] = get_mfile_name();
pm.dir.model_runs = dir_of_this_mfile;
pm.dir.data = [pm.dir.model_runs, 'data', '/'];
% directory to save model output
pm.dir.model_save = ''; % to be set


%% paths (which will be added to the matlab path)
% all paths in pm.path will be added to the matlab path (unless ==[])

% path of sourceterm functions
pm.path.sourceterm_fns = [];
% path of topography functions
pm.path.topo_fns = [];
% ICs
% directory of IC functions
pm.path.IC_fns = [];
% directory of BC functions
pm.path.BC_fns = [];


%% some file names
% the path and name of problem specific parameter mfile (to be set there)
pm.file.para_mfile = add_mfile_name_to_cellarr(pm.file.para_mfile);  %this is a cell array with all the parameter m-files


%  PYHSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%
u_bed_choice = {'const', 'fn_xy_t', 'fn_phi'};
pp.flags.u_bed = u_bed_choice{1};

% load defaults:
addpath([dir_of_this_mfile, '../parameters/'])
physical_parameters;
other_parameters;

pp.rho_w = rho_w; %  density of water
pp.rho_i = rho_i; %  density of ice
pp.g_grav = g_grav;  %  gravitational acceleration
pp.L_fusion = L_fusion; % latent heat of fusion
pp.c_w = c_w; % heat capacity of water
pp.c_t_c = -c_t; % pressure melting coefficient
pp.n_glen = n_glen; % Glen's exponent

pp.creep_const_s = A;
pp.creep_const_s_soft = A;
pp.creep_const_c = A;

% these can either be constant in space or defined at nodes or edges
pp.l_bed = l_r;     % bed undulation wave length (=l_r in write-up) (nodes)
pp.l_c = l_c;       % width of sheet which contributes to channel melt (m) (edges)
pp.h_bed = h_r;     % bed undulation height (=h_r in write-up) (nodes)

% water flow
% $$$ pp.alpha_s = alpha;  % 1 is linear, 5/4 default
% $$$ pp.beta_s = beta;   % 2 is linear, 3/2 default (turbulent flow)
% $$$ pp.alpha_c = alpha;
% $$$ pp.beta_c = beta;
pp.cond_s = k_s;  %  conductivity of sheet
pp.cond_c = k_c;  %  conductivity of channels

pp.e_v = e_v;

%  NUMERICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{2};

% SCALEd PARAMETERS
%%%%%%%%%%%%%%%%%%%%
% $$$
% $$$
% $$$ %% INPUT FUNCTIONS
% $$$ %%%%%%%%%%%%%%%%%
% $$$
% $$$ %% BC
% $$$ % Dirichlet BC for phi: applied at nodes where bmark is odd
% $$$ pin.bc_dirichlet_phi = make_anon_fn('@(xy, time, bmark, phi_0, phi_m) double(0*phi_m)', pp);
% $$$ % Flux BC for phi and h_w: i.e. set phi or h_w such that this flux is
% $$$ % given. Applied at edges where bmark_edge is even
% $$$ % zero flux:
% $$$ pin.bc_flux = make_anon_fn('@(xy, time, bmark_edge) double(zeros(sum(~logical(mod(bmark_edge,2)) & bmark_edge>0),1))');
% $$$
% $$$ %% IC
% $$$ % initial sheet thickness
% $$$ pin.ic_h = make_anon_fn('@(xy, time) double(0.05 + 0*xy(:,1))');
% $$$ % initial channel cross sectional area
% $$$ pin.ic_S =  make_anon_fn('@(xy, time) double(0.0 + 0*xy(:,1))');
% $$$
% $$$ %% Geometry and derived (constant) function
% $$$ % these are all functions of xy and time only
% $$$
% $$$ %% Source functions
% $$$
% $$$
% $$$ %% Storage terms: no storage
% $$$ % englacial storage (set to zero)
% $$$ pin.sigma = make_anon_fn('@(xy, time) double(xy(:,1).*0 + pp.e_v)', pp);  % englacial void ratio
% $$$
% $$$ % SCALING PARAMETERS
% $$$ %%%%%%%%%%%%%%%%%%%%
% $$$ % now that we have the dmesh we can get the length scale
% $$$ ps.x = 0.5*(abs(diff(dmesh.x_extent)) + abs(diff(dmesh.y_extent))); % length (m)
% $$$ ps.x_offset = dmesh.x_extent(1);
% $$$ ps.y_offset = dmesh.y_extent(1);
% $$$
% $$$ % the rest of the scales are calculated:
% $$$ ps = set_default_scales(ps, pp, dmesh);
% $$$ % if any scales are to be changed from what set_default_scales sets
% $$$ % them to do this here:
% $$$
% $$$ % SCALE EVERYTHING
% $$$ %%%%%%%%%%%%%%%%%%
% $$$ % scaling of the rest
% $$$ [psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);

% END
%%%%%
para.model = pm;
para.numeric = pn;
para.physical = pp;
para.input = pin;
para.mesh = mesh;
para.dmesh = dmesh;
para.scale = ps;
para.time = pt;