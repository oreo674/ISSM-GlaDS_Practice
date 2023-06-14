% Runs all the scenarios.

addpath(fullfile(pwd,'../input_functions/source/matlab/'))
addpath(fullfile(pwd,'../input_functions/topography/matlab/'))
addpath(fullfile(pwd,'box100by20'))
addpath(fullfile(pwd,'box100by20/steady/'))
addpath(fullfile(pwd,'box100by20/seasonal/'))
addpath(fullfile(pwd,'box100by20/diurnal/'))
addpath(fullfile(pwd,'valley/'))

% These are global variables.
meshnr = 4
resave = true % just load the .mat files and re-save the .nc files

% Scenarios

scens = {
    % sqrt
    {'run_A', fullfile(pwd,'box100by20/steady/')}
    {'run_B', fullfile(pwd,'box100by20/steady/')}
    {'run_C', fullfile(pwd,'box100by20/diurnal/')}
    {'run_D', fullfile(pwd,'box100by20/seasonal/')}
    % valley'
    {'run_E', fullfile(pwd,'valley/')}
    {'run_F', fullfile(pwd,'valley/')}
        };

torun = 1:length(scens);
for i=torun
    scen = scens{i}{1};
    cf = scens{i}{2};
    cd(cf);
% $$$     pwd
    scen
    eval(scen)
end
cd ..