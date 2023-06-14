addpath(fullfile(pwd,'..'))

if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end
if ~exist('resave', 'var') || isempty(resave)
    % load mat-file and resafe .nc
    resave = false;
end


channels = false;

for scenario = 1:6
    if resave
        steady_out = load(['outputmat/sqrt_noch',num2str(scenario),'_mesh',num2str(meshnr),'.mat']);
        save_model_run_netcdf(steady_out, ['outputnc/sqrt_noch',num2str(scenario),'.nc'], length(steady_out.para.time.out_t), true);
    else % run
        para = get_para_box100by20(channels, meshnr, scenario);
        steady_out = run_model(para);
        %% unwrap the model output structure:
        [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
        [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
        N_eff = pp_get_N(phis, fields.nodes.phi_0, para.mesh)*para.scale.phi/1e6;
        if channels
            error('no channels!')
        else
            save_model_run(steady_out, ['sqrt_noch',num2str(scenario),'_mesh',num2str(meshnr),'.mat']);
            save_model_run_netcdf(steady_out, ['outputnc/sqrt_noch',num2str(scenario),'.nc'], length(steady_out.para.time.out_t), true);
        end
    end
end
