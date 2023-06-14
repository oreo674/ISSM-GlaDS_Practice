addpath(fullfile(pwd,'..'))

if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end
if ~exist('resave', 'var') || isempty(resave)
    % load mat-file and resafe .nc
    resave = false;
end

e_v =  1e-4;
if e_v~=0
    ev_str = '_ev';
else
    ev_str = '';
end
resave = 1
dTs = [-4 -2 0 2 4]; %(0:4)-3
for i=1:length(dTs)
    if resave
        ['output/sqrt_seasonal',num2str(i), '_mesh', num2str(meshnr), ev_str ,'.mat']
        ['outputnc/sqrt_seasonal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.nc']
        steady_out = load(['output/sqrt_seasonal',num2str(i), '_mesh', num2str(meshnr), ev_str ,'.mat']);
        end_ = size(steady_out.phis,2);
        inds = end_-365:end_-1;
        shmip_save_as_netcdf(steady_out, ['outputnc/sqrt_seasonal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.nc'], inds, true);
    else
        para = get_para_seasonal(meshnr, dTs(i), e_v);
        steady_out = run_model(para);
        %% unwrap the model output structure:
        [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
        [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
        %% save
        save_model_run(steady_out, ['sqrt_seasonal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.mat']);
        end_ = size(steady_out.phis,2);
        inds = end_-364:end_;
        shmip_save_as_netcdf(steady_out, ['outputnc/sqrt_seasonal',num2str(i), '_mesh', num2str(meshnr), ev_str,'.nc'], inds, true);
    end
end
