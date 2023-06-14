if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end
if ~exist('resave', 'var') || isempty(resave)
    % load mat-file and resafe .nc
    resave = false;
end

c_t = 0 %-7.5000e-08;
if c_t~=0
    ct_str = '';
else
    ct_str = '_no-ct';
end
resave = 1
bed_paras = [300/6e3, 0, -0.1, -0.5, -0.7];
for i=1:length(bed_paras)
    if resave
        ['output/valley_steady',num2str(i), '_mesh', num2str(meshnr), ct_str, '.mat']
        ['outputnc/valley_steady',num2str(i), '_mesh', num2str(meshnr), ct_str, '.nc']
        steady_out = load(['output/valley_steady',num2str(i), '_mesh', num2str(meshnr), ct_str, '.mat']);
        shmip_save_as_netcdf(steady_out, ['outputnc/valley_steady',num2str(i), '_mesh', num2str(meshnr), ct_str, '.nc'], 'end', true);
    else
        ind = i
        para = get_para_valley(meshnr, bed_paras(i), c_t);

        steady_out = run_model(para);
        %% unwrap the model output structure:
        [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
        [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
        %% save
        save_model_run(steady_out, ['valley_steady',num2str(i), '_mesh', num2str(meshnr), ct_str, '.mat']);
        shmip_save_as_netcdf(steady_out, ['outputnc/valley_steady',num2str(i), '_mesh', num2str(meshnr), ct_str, '.nc'], 'end', true);
    end
end
