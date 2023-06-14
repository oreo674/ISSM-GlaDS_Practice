function para = get_para_moulins(meshnr, moulin_scenario)

% para = get_para_box100by20(ch, meshnr, melt_scenario, bedslope);

addpath('..')
para = get_para_box100by20(true, meshnr, 1, 0, 1); % this includes small basal melt of B1
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;
rmpath('..')

%% Model description
% a descriptive string for the model-run
pm.model_run_descript = ['Hydro intercomparsion: box100by20 with moulins'];


pp.e_v = 1e-4;

%% Source functions
% pin.source_term_s is set above to a small basal melt
moulins = csvread(['../../../input_functions/source/B',num2str(moulin_scenario), '_M.csv']);
tmpvec = zeros(dmesh.tri.n_nodes,1);
inds = map_point_onto_mesh(moulins(:,2:3), dmesh.tri.nodes);
tmpvec(inds) = moulins(:,4);
% $$$ pmd.moulins.inds = inds;
% $$$ pmd.moulins.x_sec_areas = inds*0;
pin.source_term_c = make_anon_fn('@(time) double(tmpvec)', tmpvec);

%% Storage terms: no storage
% englacial storage (set to zero)
pin.sigma = make_anon_fn('@(xy, time) double(xy(:,1).*0 + pp.e_v)', pp);  % englacial void ratio

%  NUMERICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{5};


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
