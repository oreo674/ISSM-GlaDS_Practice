function out = moulin_source(t, p)
% Makes a diurnally varying source.  para must contain:
%  amp -- relative amplitude

%t/p.day
if t<0 % for spin-up
    out = p.M_steady;
else
    out = p.M_steady .* (1 - p.amp.*sin(2.*pi.*t./p.day));
    out(out<0) = 0;
end

% add basal melt
out = out + 7.93e-11 * p.width;