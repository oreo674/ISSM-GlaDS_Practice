function temp = seasonal_temp(time, deltaT)
% temp = seasonal_temp(time, deltaT)
%
% Returns temperature at zero elevation mimicking a seasonal variation using a sine curve.
%
% time -- time since 1 Jan (sec)
% deltaT -- temperature offset: shift all temperatures by this much
%
% out:
% temp -- temperature at 0 elevation (C)
%
% Hard coded:
% - max_temp: maximal temperature (before deltaT offset is applied): 11C
% - mean_temp: average temperature

year = 365*24*60*60;

% this mimics Kangerlussuaq
max_temp = 11;
mean_temp = -5;

amp = max_temp-mean_temp;

temp = -amp*cos(2*pi/year*time) + mean_temp + deltaT;
