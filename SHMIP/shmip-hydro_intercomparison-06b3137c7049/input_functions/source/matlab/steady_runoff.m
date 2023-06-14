function m=steady_runoff(scenario_nr)
% m=steady_runoff(scenario_nr)
%
% Returns the steady inputs for model run suites A and B:
%
% | A1       | 7.93e-11 m/s |
% | A2       | 1.59e-09 m/s |
% | A3       | 5.79e-09 m/s |
% | A4       | 2.5e-08 m/s |
% | A5       | 4.5e-08 m/s |
% | A6       | 5.79e-07 m/s |
%
% Input:
% - number for the scenario
%
% Output:
% - source (m/s)

%tmp = [7.93e-11, 1.59e-09, 5.79e-09,  1.74e-08,  5.79e-08,  5.79e-07];
tmp = [7.93e-11, 1.59e-09, 5.79e-09,  2.5e-08,  4.5e-08,  5.79e-07];
m = tmp(scenario_nr);
