function m=diurnal_runoff(time, relative_amp, steady_state_runoff, basal_melt)
% m=diurnal_runoff(time, relative_amp, steady_state_runoff, basal_melt)
%
% Puts a diurnal signal onto a steady one.  If time is negative, no diurnal signal will be
% added (for spin-up).  If runoff is negative, it is set to zero.  Some small basal melt
% is always present if basal_melt==true.
%
% Input:
%  - time (s)
%  - relative_amp: amp = relative_amp*steady_state_runoff
%  - steady_state_runoff: runoff without diurnal variations
%  - basal_melt: include basal melt if == true

if basal_melt
    basal = steady_runoff(1);
else
    basal = 0;
end

day = 24*60*60;

if time<0
    m = steady_state_runoff;
else
    m = steady_state_runoff*(1 - relative_amp*sin(2*pi*time/day)) ;
end

m(m<0) = 0;

% add basal
m = m + basal;