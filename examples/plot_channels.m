% addpath('/home/tim/SFU-code/ISSM-GlaDS-test/box-steady') % Just for my prototyping
set_paths;
load BoxSteady_template

% Important to take abs(Q) since Q can be >0 and <0 depending on the
% orientation of the edge
cdata = abs(md.results.TransientSolution(end).ChannelDischarge);

figure
hold on
% Set colormap (try changing this to see what looks good)
cmap = flipud(bone);

plotchannels(md, cdata, 'colormap', cmap, 'linewidth', 1)

% Good habits for making figures
axis image
xlim([0, 100e3])
ylim([0, 25e3])
xlabel('Distance upglacier (m)')
ylabel('Distance across-glacier (m)')
title('Channel discharge (m^3 s^{-1})')


