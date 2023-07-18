model_name = 'KWsyntheticv2';
output_fname = [model_name, '.mat'];
md_struct = load(output_fname);
md = md_struct.md;
md.mesh.x = md.mesh.x/1e3;
md.mesh.y = md.mesh.y/1e3;
nodes = [md.mesh.x, md.mesh.y];

% PLOT GEOMETRY
figure;
xc = 0:0.1:30;
ff_bed = scatteredInterpolant(md.mesh.x, md.mesh.y, md.geometry.bed);
ff_thick = scatteredInterpolant(md.mesh.x, md.mesh.y, md.geometry.thickness);
bed_center = ff_bed(xc, 0*xc);
thick_center = ff_thick(xc, 0*xc);
plot(xc, bed_center);
hold on
plot(xc, thick_center + bed_center);
grid on
xlabel('x (km)')
ylabel('z (m)')
yyaxis right
plot(xc, thick_center)
ylabel('Thickness (m)')
legend({'Bed', 'Surface', 'Thickness'})


% Simple scatter plot
N = [md.results.TransientSolution.EffectivePressure];
h_sheet = [md.results.TransientSolution.HydrologySheetThickness];
phi = [md.results.TransientSolution.HydraulicPotential];
time = [md.results.TransientSolution.time];
phi_base = md.materials.rho_water*md.constants.g*md.geometry.bed;
% pw = phi - phi_base;
floatation = pw./(N + pw);
figure;
scatter(md.mesh.x, pw(:, end)/1e6)
xlabel('x (km)')
ylabel('p_w (MPa)')


% Messy but need to compute the nodes adjacent to each element
connect = zeros(md.mesh.numberofelements, 3);
for ii=1:md.mesh.numberofelements
    [row, col] = find(md.mesh.vertexconnectivity(:, 1:100)==ii);
    connect(ii, :) = row;
end

figure('Units', 'inches', 'Position', [2, 3, 8, 4])

% Need to assign a name to the tiled layout for the colorbar later
T = tiledlayout(2, 1, 'Padding', 'Compact', 'TileSpacing', 'compact');

% WATER PRESSURE
ax1 = nexttile(1);
hold on
patch('Faces', connect, 'Vertices', nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', floatation(:, end), 'FaceColor', 'flat')
axis image
% cmocean('balance', 'Pivot', 0)
cmocean('dense')
clim(ax1, [0, 1])
xlim(ax1, [0, 30])
ylim(ax1, [-2.5, 2.5])
yticks(ax1, [0, 12.5, 25])
% clim(ax1, [0, 0.1])

% Add channels
Q_cmap = cmocean('turbid');
Qmin = 0.1; Qmax = 5;
plotchannels(md, abs(md.results.TransientSolution(end-10).ChannelDischarge),...
        'colormap', Q_cmap, 'min', Qmin, 'max', Qmax)

% SHEET THICKNESS
ax2 = nexttile(2);
hold on
patch('Faces', connect, 'Vertices', nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', h_sheet(:, end), 'FaceColor', 'flat')
axis image
cmocean('matter')
xlim(ax2, [0, 30])
ylim(ax2, [-2.5, 2.5])
clim(ax2, [0, 0.1])

% Add channels
plotchannels(md, abs(md.results.TransientSolution(end).ChannelDischarge),...
        'colormap', Q_cmap, 'min', Qmin, 'max', Qmax)

% Effective pressure colorbar
cb = colorbar(ax1);
% cb.Layout.Tile = 'east';
cb.Label.String = 'p_w';

cbh = colorbar(ax2);
cbh.Label.String = 'h_s';
% cbh.Layout.Tile = 

% Add secondary invisible axes for the second colorbar
cax = axes(T);
cax.Visible=false;
colormap(cax, Q_cmap);
cb2 = colorbar(cax);
cb2.Layout.Tile = 'North';
cb2.Label.String = 'Q (m^3 s^{-1})';
cb2.Ticks = [Qmin, 0.5, 1, 2, 3, 4, Qmax];
clim(cax, [Qmin, Qmax])

figure;
plot(time, mean(pw, 1))
grid on
xlabel('Year')
ylabel('Mean p_w (MPa)')