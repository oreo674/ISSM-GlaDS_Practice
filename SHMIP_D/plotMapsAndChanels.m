%% Working example of multiple map panels with two colorbars
% See last example here:
% https://www.mathworks.com/help/matlab/creating_plots/graph-with-multiple-x-axes-and-y-axes.html
​
results = load('D5check.mat');
md = results.md;
md.mesh.x = md.mesh.x/1e3;
md.mesh.y = md.mesh.y/1e3;
nodes = [md.mesh.x, md.mesh.y];
N = [md.results.TransientSolution.EffectivePressure]/1e6;
​
% Messy but need to compute the nodes adjacent to each element
connect = zeros(md.mesh.numberofelements, 3);
for ii=1:md.mesh.numberofelements
    [row, col] = find(md.mesh.vertexconnectivity(:, 1:100)==ii);
    connect(ii, :) = row;
end
​
figure
​
% Need to assign a name to the tiled layout for the colorbar later
T = tiledlayout(2, 1);
​
% Plot effective pressure at first timestep
ax1 = nexttile(1);
hold on
patch('Faces', connect, 'Vertices', nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', N(:, 10), 'FaceColor', 'flat')
axis image
cmocean('deep')
xlim(ax1, [0, 100])
ylim(ax1, [0, 25])
yticks(ax1, [0, 12.5, 25])
clim(ax1, [0, 0.1])
​
% Add channels
Q_cmap = cmocean('turbid');
Qmin = 1; Qmax = 250;
plotchannels(md, abs(md.results.TransientSolution(end).ChannelDischarge),...
        'colormap', Q_cmap, 'min', Qmin, 'max', Qmax)
​
​
% Plot effective pressure at next timestep
ax2 = nexttile(2);
hold on
patch('Faces', connect, 'Vertices', nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', N(:, 20), 'FaceColor', 'flat')
axis image
cmocean('deep')
xlim(ax2, [0, 100])
ylim(ax2, [0, 25])
yticks(ax2, [0, 12.5, 25])
clim(ax2, [0, 0.1])
clim(ax2, [0, 0.1])
​
% Add channels
Q_cmap = cmocean('turbid');
Qmin = 1; Qmax = 250;
plotchannels(md, abs(md.results.TransientSolution(end-30).ChannelDischarge),...
        'colormap', Q_cmap, 'min', Qmin, 'max', Qmax)
​
% Effective pressure colorbar
cb = colorbar(ax2);
cb.Layout.Tile = 'east';
cb.Label.String = 'N (MPa)';
​
% Add secondary invisible axes for the second colorbar
cax = axes(T);
cax.Visible=false;
colormap(cax, Q_cmap);
cb2 = colorbar(cax);
cb2.Layout.Tile = 'North';
cb2.Label.String = 'Q (m^3 s^{-1})';
cb2.Ticks = [1, 50, 100, 150, 200, 250];
clim(cax, [Qmin, Qmax])