function para = get_para_diurnal(meshnr, relative_amp, IC_scenario, e_v);
% para = get_para_diurnal(meshnr, relative_amp, IC_scenario, e_v);

addpath('..')
para = get_para_box100by20(1,meshnr,1);
clear para;
rmpath('..')

% IC from other run
base_scenario = 5;
infile = ['../steady/outputmat/sqrt_moulins', num2str(base_scenario), '_mesh', num2str(meshnr),'.mat'];
load(infile, 'para', 'S_channels', 'h_sheets');
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);
clear para;

% directory to save model output
pm.dir.model_save = 'output/'; % to be set

% TIME
%%%%%%
pt.start = 0;        % start time
pt.end   = 50*pp.day;  % end time
pt.out_t = 0*pp.day:pp.hour:pt.end;
%pt.out_t = [0, 48*pp.day:pp.hour:pt.end];

pp.relative_amp = relative_amp;

if ~exist('e_v', 'var') || isempty(e_v)
    if strcmp(pn.ts.stepper,'ode15s')
        pp.e_v = 1e-7; % this seems to have 0 impact on the result
    else
        pp.e_v = 0;
    end
else
    pp.e_v = e_v;
end

%% Source functions
switch IC_scenario
  case 'A'
    steady_input = pin.source_term_s(dmesh.tri.nodes, 0);
    pin.source_term_s = make_anon_fn('@(xy, time) double(diurnal_runoff(time, relative_amp, steady_input, true))', relative_amp, steady_input);
    tmpvec = zeros(dmesh.tri.n_nodes,1);
    pin.source_term_c = make_anon_fn('@(time) double(tmpvec*0)', tmpvec);
  case 'B'
    % keep pin.source_term_s (basal melt)
    steady_input = pin.source_term_c(0);
    pin.source_term_c = make_anon_fn('@(time) double(diurnal_runoff(time, relative_amp, steady_input, false))', relative_amp, steady_input);
  otherwise
    error('Need to use IC from scenario B or E')
end

%% Storage terms: no storage
% englacial storage (set to zero)
pin.sigma = make_anon_fn('@(xy, time) double(xy(:,1).*0 + pp.e_v)', pp);  % englacial void ratio

% ICs (note only IC for h_sheet and S_channel are used, not phi)
% path where IC are stored

% $$$ % initial sheet thickness
% $$$ h_ic = h_sheets(:,end)*ps.h;
% $$$ pin.ic_h = make_anon_fn('@(xy, time) double(h_ic)', h_ic);
% $$$ % initial channel cross sectional area
% $$$ S_ic = S_channels(:,end)*ps.S;
% $$$ pin.ic_S =  make_anon_fn('@(xy, time) double(S_ic)', S_ic);

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
pn.ts.ode113.opts.MaxStep = 0.005; % in scaled units
pn.ts.ode113.opts.RelTol = 1e-8;
pn.ts.ode113.opts.AbsTol = 1e-8;
pn.ts.ode15s.opts.MaxStep = 0.005;
pn.ts.ode15s.opts.RelTol = 1e-8;
pn.ts.ode15s.opts.AbsTol = 1e-8;

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
