function para = get_para_seasonal(meshnr, dT, e_v);
%  para = get_para_seasonal(meshnr, dT, e_v);

addpath('..')
para = get_para_box100by20(1,meshnr,1);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;
rmpath('..')

% IC from other run
base_scenario = 1;
infile = ['../steady/outputmat/sqrt_ch',num2str(base_scenario), '_mesh',num2str(meshnr),'.mat'];
load(infile, 'para', 'S_channels', 'h_sheets');
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;


% directory to save model output
pm.dir.model_save = 'output/'; % to be set

% TIME
%%%%%%
pt.start = 0;        % start time
pt.end   = 4*365*pp.day;  % end time
%pt.out_t = 0*pp.day+1:pp.day:pt.end;
pt.out_t = [0:pp.day*10:3*365*pp.day-1, 3*365*pp.day:pp.day:pt.end];

pp.source_dT = dT;
if ~exist('e_v', 'var') || isempty(e_v)
    pp.e_v = 0*1e-4;
else
    pp.e_v = e_v;
end

%% Source functions
ele_fn = pin.ice_thickness;
pin.source_term_s = make_anon_fn('@(xy, time) double(seasonal_runoff(xy, time, dT, ele_fn))', ele_fn, dT);
tmpvec = zeros(dmesh.tri.n_nodes,1);
pin.source_term_c = make_anon_fn('@(time) double(tmpvec*0)', tmpvec);

%% Storage term
pin.sigma = make_anon_fn('@(xy, time) double(xy(:,1).*0 + pp.e_v)', pp);  % englacial void ratio

%% ICs
pm.IC_from_file = 1;
pm.file.IC_file = infile;
pm.IC_file_timestep = 'end'; % from what time step to take the IC
pm.hotstart = 0;  % hotstart is when one contiuous a model run as opposed to start a new one with IC from a old one

%  NUMERICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
if pp.e_v>0
    pn.ts.stepper =  steppers{5};
else
    pn.ts.stepper =  steppers{2};
end

% note, below are in scaled time:
pn.ts.ode113.opts.MaxStep = 0.5;
pn.ts.ode15s.opts.MaxStep = 0.5;

% SCALING PARAMETERS
%%%%%%%%%%%%%%%%%%%%
ps = set_default_scales(ps, pp, dmesh);

% SCALE EVERYTHING
%%%%%%%%%%%%%%%%%%
% scaling of the rest
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);

% END
%%%%%
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
