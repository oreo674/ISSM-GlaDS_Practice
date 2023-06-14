if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end

i=1;
bed_paras = [300/6e3];
for i=1
    para = get_para_valley(meshnr, bed_paras(i), 1);

    steady_out = run_model(para);
    %% unwrap the model output structure:
    [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
    [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
    %% save
    save_model_run(steady_out, ['basal',num2str(i), '_mesh', num2str(meshnr) ,'.mat']);
    save_model_run_netcdf(steady_out, ['outputnc/basal',num2str(i),'.nc'], [], true);
    i = i+1;
end
