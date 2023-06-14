addpath(fullfile(pwd,'..'))

if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end
if ~exist('resave', 'var') || isempty(resave)
    % load mat-file and resafe .nc
    resave = false;
end

channels = true;

for m_scenario = 1:5
    if resave
        steady_out = load(['outputmat/sqrt_moulins',num2str(m_scenario),'_mesh',num2str(meshnr),'.mat']);
        shmip_save_as_netcdf(steady_out, ['outputnc/sqrt_moulins',num2str(m_scenario),'_mesh',num2str(meshnr),'.nc'], 'end', true);
    else
        ind = m_scenario
        para = get_para_moulins(meshnr, m_scenario);
        steady_out = run_model(para);
        %% unwrap the model output structure:
        [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
        [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
        N_eff = pp_get_N(phis, fields.nodes.phi_0, para.mesh)*para.scale.phi/1e6;

        save_model_run(steady_out, ['sqrt_moulins',num2str(m_scenario),'_mesh',num2str(meshnr),'.mat']);
        shmip_save_as_netcdf(steady_out, ['outputnc/sqrt_moulins',num2str(m_scenario),'_mesh',num2str(meshnr),'.nc'], 'end', true);
    end
end
