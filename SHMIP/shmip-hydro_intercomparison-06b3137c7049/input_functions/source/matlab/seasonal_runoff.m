function runoff = seasonal_runoff(xy, time, deltaT, ele_fn)
% runoff = seasonal_runoff(xy, time, deltaT, ele_fn)
%
% Calculates the runoff for the seasonal temperature forcing.

basal = steady_runoff(1);

% zero runoff if time is negative: for spin-up
if time<0
    runoff = xy(:,1) * 0;
end
temp = seasonal_temp(time, deltaT);
ele = ele_fn(xy);
runoff = ddf_model(ele, temp);

runoff = runoff + basal;
