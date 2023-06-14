addpath(fullfile(pwd,'..'))

if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end
if ~exist('resave', 'var') || isempty(resave)
    % load mat-file and resafe .nc
    resave = false;
end

e_v = 0; % 1e-4;
if e_v~=0
    ev_str = '_ev';
else
    ev_str = '';
end
resave = 1
relative_amp_s = [1/4, 3/4, 1, 2];
for i=1:length(relative_amp_s)
    if resave
        steady_out = load(['output/sqrt_moulins_diurnal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.mat']);
        end_ = size(steady_out.phis,2);
        inds = end_-24:end_-1;
        shmip_save_as_netcdf(steady_out, ['outputnc/sqrt_moulins_diurnal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.nc'], inds, true);
    else
        relative_amp = relative_amp_s(i);
        para = get_para_diurnal(meshnr, relative_amp, 'B', e_v);

        steady_out = run_model(para);
        %% unwrap the model output structure:
        [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
        [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
        %% save
        save_model_run(steady_out, ['sqrt_moulins_diurnal',num2str(i), '_mesh', num2str(meshnr), ev_str ,'.mat']);
        end_ = size(steady_out.phis,2);
        inds = end_-24:end_-1;
        shmip_save_as_netcdf(steady_out, ['outputnc/sqrt_moulins_diurnal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.nc'], inds, true);
    end
end
