error('outdated?')

channels = true;
meshnr = 3;

for scenario = 1:6
    steady_out = load(['scenario1/ch',num2str(scenario),'.mat']);
    save_model_run_netcdf(steady_out, ['scenario1/ch',num2str(scenario),'.nc'], [], true);
end

for scenario = 1:6
    steady_out = load(['scenario1/noch',num2str(scenario) ,'.mat']);
    save_model_run_netcdf(steady_out, ['scenario1/noch',num2str(scenario) ,'.nc'], [], true);
end
