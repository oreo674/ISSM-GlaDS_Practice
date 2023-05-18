set_paths;
load BoxSteady_template


cdata = abs(md.results.TransientSolution(end).ChannelDischarge);
cmap = flipud(bone);

figure
hold on
% Plot on the current axes (gca)
%  from model outputs md
%  plot data cdata (channel discharge in this example)
%  colormap = cmap (set to bone)
%  Set colorbar limits to [0, 100]
%  'vmin', 0 shows all channels. 'vmin', 1 would mask any channels
%  with discharge less than 1, for example
edge_plot_ISSM(gca, md, cdata, cmap, [0, 100], 'vmin', 0)
xlabel('x (m)')
ylabel('y (m)')
title('Q (m^3 s^{-1})')

%% An example of how to do this manually
% figure
% hold on
% for ii=1:dmesh.tri.n_edges
%     n1 = dmesh.tri.connect_edge(ii, 1);
%     n2 = dmesh.tri.connect_edge(ii, 2);
%     xs = [dmesh.tri.nodes(n1, 1), dmesh.tri.nodes(n2, 1)];
%     ys = [dmesh.tri.nodes(n1, 2), dmesh.tri.nodes(n2, 2)];
%     if cdata(ii)>0.1
%         color = [1, 0, 0];
%         lw = 2;
%     else
%         color = [0, 0, 0];
%         lw = 0.5;
%     end
%     plot(xs, ys, 'Color', color, 'LineWidth', lw)
% end
% axis image
