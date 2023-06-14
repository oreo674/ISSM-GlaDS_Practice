function para = get_para_valley_seasonal(meshnr, dT)

para = get_para_valley(meshnr, 300/6e3);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;

% IC from other run
base_scenario = 1;
infile = ['output/basal1_mesh',num2str(meshnr),'.mat'];
load(infile, 'para', 'S_channels', 'h_sheets');
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;

%% Model description
% a descriptive string for the model-run
pm.model_run_descript = ['Hydro intercomparsion: valley seasonal'];

% TIME
%%%%%%
pt.start = 0;        % start time
pt.end   = 2*365*pp.day;  % end time
pt.out_t = 0:pp.day:pt.end;

pp.source_dT = dT;

%% Source functions
ele_fn = make_anon_fn('@(xy) double(pin.ice_thickness(xy,0) + pin.bed_elevation(xy,0) )', pin);
pin.source_term_s = make_anon_fn('@(xy, time) double(seasonal_runoff(xy, time, dT, ele_fn))', ele_fn, dT);
tmpvec = zeros(dmesh.tri.n_nodes,1);
pin.source_term_c = make_anon_fn('@(time) double(tmpvec*0)', tmpvec);

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
