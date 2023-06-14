function runoff = ddf_model_diurnal(ele, temp0, time)
% runoff = ddf_model_dirunal(ele, temp0, time)
%
% A degree day factor (DDF) melt model which modulates the mean daily melt with a
% diurnally varying signal.
%
% ele -- elevation (m)
% temp0 -- "mean daily" temperature (C)
% time -- time (sec)
%
% out:
% runoff -- water runoff per unit area (m/s)
%
% Contains hard-coded:
% - degree day factor: 0.01 m/d/K
% - lapse_rate = -0.0075 K/m
% - elevation_zero = 0m
% - ratio_amp = 1
%     -- ratio diurnal amplitude/mean runoff


day = 24*60*60;
lapse_rate = -0.0075;
DDF = 0.01/day; % degree day factor
elevation_zero = 0;
ratio_amp = 1;

runoff = ((ele-elevation_zero)*lapse_rate+temp0)*DDF;
runoff(runoff<0) = 0;

% add diurnals
amps = runoff*ratio_amp;
runoff = runoff - amps.*cos(2*pi/day * time);

% if runoff is negative, set to zero (only possible if ratio_amp>1)
runoff(runoff<0) = 0;
