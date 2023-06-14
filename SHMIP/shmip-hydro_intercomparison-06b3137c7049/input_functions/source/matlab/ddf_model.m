function runoff = ddf_model(ele, temp0)
% A degree day factor (DDF) melt model.
%
% ele -- elevation (m)
% temp0 -- temperature (C)
%
% out:
% runoff -- water runoff per unit area (m/s)
%
% Contains hard-coded:
% - degree day factor: 0.01 m/d/K
% - lapse_rate = -0.0075 K/m
% - elevation_zero = 0m


day = 24*60*60;
lapse_rate = -0.0075;
DDF = 0.01/day; % degree day factor
elevation_zero = 0;

runoff = ((ele-elevation_zero)*lapse_rate+temp0)*DDF;
runoff(runoff<0) = 0;
