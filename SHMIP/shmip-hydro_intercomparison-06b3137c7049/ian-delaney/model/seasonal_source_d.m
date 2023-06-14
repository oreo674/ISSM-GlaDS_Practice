function out = seasonal_source_d(t, p)
% Makes a seasonally varying source.

temp = -16 * cos(2*pi*t/p.year)-5 + p.dt;

% lapse it to all elevations
out = (p.lr .* p.H +temp)*p.ddf;

out = out.*p.width;

out(out<0) = 0;

% add basal melt
out = out + p.basal;
