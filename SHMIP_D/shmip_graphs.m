shmip_3 = load('test.mat');
close all;
N_3 = [shmip_3.md.results.TransientSolution.EffectivePressure]/1e6;
% h = [shmip_3.md.results.TransientSolution.HydrologySheetThickness];
% q = [shmip.md.results.TransientSolution.];
Q_3 = abs([shmip_3.md.results.TransientSolution.ChannelDischarge]);
tt3 = [shmip_3.md.results.TransientSolution.time];


shmip4 = load('D2R_1e-4.mat');
N_4 = [shmip4.md.results.TransientSolution.EffectivePressure];
% Q_4 = abs([shmip_4.md.results.TransientSolution.ChannelDischarge]);
tt4 = [shmip4.md.results.TransientSolution.time];

%%% SENSITIVITY TESTING
% figure;
% xlabel('Time');
% ylabel('Mean Effective Pressure (MPa)');
% plot(tt3, mean(N_3,1)/1e6);
% 
% hold on;
plot(tt4, mean(N_4,1)/1e6);
xlabel('Time');
ylabel('Effective Pressure MPa');
title('Mean N vs t for D2');
% legend({'ResTol = 1e-4','ResTol = 1e-3'},'Location','southeast');


%% scatter plotting N vs x for various times (5 x 5 = 25 days)
% shmip3x = shmip_3.md.mesh.x/1e3;
% figure;
% subplot(3,3,1)
% scatter(shmip3x, N_3(:, 1780));
% xlabel('Distance from terminus (km)');
% ylabel('Effectie Pressure (MPa)');
% 
% 
% subplot(3,3,2)
% scatter(shmip3x, N_3(:, 1785));
% xlabel('Distance from terminus (km)');
% ylabel('Effectie Pressure (MPa)');
% 
% subplot(3,3,3)
% scatter(shmip3x, N_3(:, 1790));
% xlabel('Distance from terminus (km)');
% ylabel('Effectie Pressure (MPa)');
% 
% subplot(3,3,4)
% scatter(shmip3x, N_3(:, 1795));
% xlabel('Distance from terminus (km)');
% ylabel('Effectie Pressure (MPa)');
% 
% subplot(3,3,5)
% scatter(shmip3x, N_3(:, 1800));
% xlabel('Distance from terminus (km)');
% ylabel('Effectie Pressure (MPa)');
% 
% subplot(3,3,6)
% scatter(shmip3x, N_3(:, 1805));
% xlabel('Distance from terminus (km)');
% ylabel('Effectie Pressure (MPa)');
% sgtitle('N vs x for 8900-9025 Days')
% print('EffPress_D1','-dpng','-r600')


%% spatial map w channels at 1790
% figure;
% plotmodel(shmip_3.md, 'data', shmip_3.md.results.TransientSolution(1790).EffectivePressure/1e6, ...
%     'colormap', 'hsv', 'xlabel', 'Distance from terminus (m)', 'ylabel', 'Width of glacier (m)', ...
%     'title', 'Spatial N for 8950 Days', 'colorbartitle', 'N (MPa)');
% hold on;
% plotchannels(shmip_3.md, abs(shmip_3.md.results.TransientSolution(1790).ChannelDischarge), ...
%     'min', 0.1);
% print('D1Spat8950Days','-dpng','-r800')
% 
% %% Plotting sheet thickness
% figure;
% plotmodel(shmip_3.md, 'data', shmip_3.md.results.TransientSolution(1790).HydrologySheetThickness, ...
%     'colormap', 'Ala', 'xlabel', 'Distance from terminus (m)', 'ylabel', 'Width of glacier (m)', ...
%     'title', 'Spatial N for 8950 Days', 'colorbartitle', 'h (m)');
% hold on;
% plotchannels(shmip_3.md, abs(shmip_3.md.results.TransientSolution(1790).ChannelDischarge), ...
%     'min', 0.1);



% figure;
% plotmodel(shmip_3.md, 'data', shmip_3.md.results.TransientSolution(1795).EffectivePressure, ...
%     'colormap', 'hsv', 'xlabel', 'Distance from terminus', 'ylabel', 'Width of glacer', ...
%     'title', 'Spatial N for 8975 Days');
% hold on;
% plotchannels(shmip_3.md, abs(shmip_3.md.results.TransientSolution(1795).ChannelDischarge), ...
%     'min', 0.1);



%%% Plotting time series
% subplot(2,1,1);
% plot(tt3, mean(N_3,1)/1e6);
% xlabel('Time');
% ylabel('Mean Effective Pressure MPa')
% title('Mean N vs. t for Test')
% % 
% subtt = tt3(1753:1826);
% subN = N_3(:, 1753:1826)
% 
% subplot(2,1,2)
% plot(tt4, mean(N_4,1)/1e6);
% xlabel('Time');
% ylabel('Mean Effective Pressure MPa');
% title('Mean N vs t for e-3 Conv');

% scatter(shmip_3.md.mesh.x/1e3, N_3(:, 1800)/1e6)
% xlabel('Distance from terminus (km)');
% ylabel('N (MPa)');
% title('N vs x for Day 9000')

% figure;
% plot(tt, max(Q_3, [], 1));
% hold on;
% plot(tt, max(Q_4, [], 1));
% xlabel('Max Q');
% ylabel('Time')
% title('Max Q vs to for D1');
% legend({'ResTol = 1e-3', 'ResTol = 1e-4'}, 'Location', 'southeast');



%%% Numerically reaching SS
% dQ = Q_4(:,2:end) - Q_4(:,1:end-1);
% dt = tt4(2:end) - tt4(1:end-1);
% dQdt = dQ/dt;
% dQdt_max = max(dQdt(:,end));
% % row = find(dQdt==dQdt_max)
% dQdt_fin = abs(dQdt(end))
% if dQdt_fin < 1e-4
%     disp('Glacier is in SS')
% end









% use subplott
