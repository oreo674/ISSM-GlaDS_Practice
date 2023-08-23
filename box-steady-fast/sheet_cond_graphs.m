%% Graphing good spatial maps

results = load('std_parameter.mat')
md = results.md;
P = [md.results.TransientSolution.EffectivePressure]/1e6;
md.mesh.x = md.mesh.x/1e3;
md.mesh.y = md.mesh.y/1e3;
nodes = [md.mesh.x, md.mesh.y];


close all
connect = zeros(md.mesh.numberofelements, 3);
for ii=1:md.mesh.numberofelements
    [row, col] = find(md.mesh.vertexconnectivity(:, 1:100)==ii);
    connect(ii, :) = row;
end

figure
ax = axes;
hold on
% Plot map of sheet thickness
patch('Faces', connect, 'Vertices', nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', P(:, end), 'FaceColor', 'flat')
axis image
cmocean('balance')
xlim([0, 100])
ylim([0, 25])

cb = colorbar(ax, 'EastOutside');
cb.Label.String = 'N (MPa)';
clim([-8, 8])
Q_cmap = cmocean('turbid');
Qmin = 1; Qmax = 20;
plotchannels(md, abs(md.results.TransientSolution(end).ChannelDischarge),...
        'colormap', Q_cmap, 'min', Qmin, 'max', Qmax)
ax2 = axes;
hold on
ax2.Visible=false;
colormap(ax2, Q_cmap);
cb2 = colorbar(ax2, 'NorthOutside');
cb2.Label.String = 'Q (m^3 s^{-1})';
cb2.Ticks = [0, 5, 10, 15, 20];
clim(ax2, [Qmin, Qmax])

linkaxes([ax, ax2])

xlabel(ax, 'Distance from terminus (km)')
ylabel(ax, 'Distance across glacier (km)')
yticks(ax, [0, 12.5, 25])
title(ax, 'Spatial Map')

%print('D5LastMelt', '-dpng', '-r600')



%%


close all
figure('Units', 'inches', 'Position', [2, 2, 10, 5])
    plotmodel(md,'data', md.results.TransientSolution(end).EffectivePressure,'title','Final N [Pa]')






