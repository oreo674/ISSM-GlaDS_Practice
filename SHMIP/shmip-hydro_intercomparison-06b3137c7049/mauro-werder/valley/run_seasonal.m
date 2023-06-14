if ~exist('meshnr', 'var') || isempty(meshnr)
    meshnr = 3;
end
if ~exist('resave', 'var') || isempty(resave)
    % load mat-file and resafe .nc
    resave = false;
end
if ~exist('make_ic', 'var') || isempty(make_ic)
    make_ic = false;
end

if ~resave && make_ic
    % make IC
    disp('Starting on IC')
    bed_para = [300/6e3];
    para = get_para_valley(meshnr, bed_para, 1);

    steady_out = run_model(para);
    %% unwrap the model output structure:
    [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
    [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
    %% save
    disp('finished IC')
    save_model_run(steady_out, ['basal1_mesh', num2str(meshnr) ,'.mat']);
% $$$ disp('finished 2')
% $$$ save_model_run_netcdf(steady_out, ['outputnc/basal',num2str(i),'.nc'], [], true);
    disp('saved IC')
end

% do runs
resave=1
dTs = [-6 -3 0 3 6] %(0:4)-3
for i=1:length(dTs)
    if resave
        steady_out = load(['output/valley_seasonal',num2str(i), '_mesh', num2str(meshnr) ,'.mat']);
        end_ = size(steady_out.phis,2);
        inds = end_-365:end_-1;
        shmip_save_as_netcdf(steady_out, ['outputnc/valley_seasonal',num2str(i), '_mesh', num2str(meshnr),'.nc'], inds, true);
    else
        ind = i
        para = get_para_valley_seasonal(meshnr, dTs(i));

        steady_out = run_model(para);
        %% unwrap the model output structure:
        [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(steady_out);
        [pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd] = unwrap_all_para(para);
        %% save
        save_model_run(steady_out, ['valley_seasonal',num2str(i), '_mesh', num2str(meshnr) ,'.mat']);
        end_ = size(steady_out.phis,2);
        inds = end_-365:end_-1;
        shmip_save_as_netcdf(steady_out, ['outputnc/valley_seasonal',num2str(i), '_mesh', num2str(meshnr),'.nc'], inds, true);
    end
end
