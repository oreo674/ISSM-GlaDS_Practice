close all
shmip1 = load('D2_hopethisworks.mat');
N1 = [shmip1.md.results.TransientSolution.EffectivePressure];
Q1 = [shmip1.md.results.TransientSolution.ChannelDischarge];
tt = [shmip1.md.results.TransientSolution.time];

 shmip2 = load('D2_hopethisworks.mat');
 N2 = [shmip2.md.results.TransientSolution.EffectivePressure];
% Q2 = [shmip2.md.results.TransientSolution.ChannelDischarge];

% shmip3 = load('');
% N3 = [shmip3.md.results.TransientSolution.EffectivePressure];
% Q3 = [shmip3.md.results.TransientSolution.figure;

% shmip4 = load('');
% N4 = [shmip4.md.results.TransientSolution.EffectivePressure];
% Q4 = [shmip4.md.results.TransientSolution.ChannelDischarge];
% 
% shmip5 = load('');
% N5 = [shmip5.md.results.TransientSolution.EffectivePressure];
% Q5 = [shmip5.md.results.TransientSolution.ChannelDischarge];
% 
% 
% time = tt((end-72):end)*365;




% yr_real = load('');
% n = [yr_real.md.results.TransientSolution.EffectivePressure];
% plotmodel(yr_real.md, 'data', yr_real.md.results.TransientSolution(end).EffectivePressure)




%%% Plotting pressure vs x
% subplot(5,1,1);
% scatter(shmip1.md.mesh.x/1e3, N1(:,end)/1e6);
% xlabel('Distance from terminus in km');
% ylabel('N in MPa');
% title('N vs x for D1');
% 
% subplot(5,1,2);
% scatter(shmip1.md.mesh.x/1e3, N2(:,end)/1e6);
% title('N vs x for D2');
% 
% subplot(5,1,3);
% scatter(shmip1.md.mesh.x/1e4, N3(:,end)/1e6);
% title('N vs x for D3');
% 
% subplot(5,1,4);
% scatter(shmip1.md.mesh.x/1e3, N4(:,end)/1e6);
% title('N vs x for D4');
% 
% subplot(5,1,5);
% scatter(shmip1.md.mesh.x/1e3, N5(:,end)/1e6);
% title('N vs x for D5');




%%% Plotting N against t

% N1mean = mean(N1, 1)/1e6;
% N2mean = mean(N2, 1)/1e6;
% N3mean = mean(N3, 1)/1e6;
% N4mean = mean(N4, 1)/1e6;
% N5mean = mean(N5, 1)/1e6;
% 
% 
% subplot(5,1,1);
% plot(time, N1mean((end-72):end));
% title('Time Series for N in D1');
% xlabel('Time in Days');
% ylabel('Effective Pressure (MPa)')
% 
% subplot(5,1,2);
% plot(time, N2mean((end-72):end));
% title('Time Series for N in D2');
% 
% subplot(5,1,3);
% plot(time, N3mean((end-72):end));
% title('Time Series for N in D3');
% 
% subplot(5,1,4);
% plot(time, N4mean((end-72):end));
% title('Time Series for N in D4');
% 
% subplot(5,1,5);
% plot(time, N5mean((end-72):end));
% title('Time Series for N in D5');
x = [shmip1.md.mesh.x]/1e3;
y = [shmip1.md.mesh.y]/1e3;
% 
% scatter(x, y, 5, N1(:,end)/1e6)


% scatter(shmip1.md.mesh.x/1e3,shmip1.md.geometry.surface)
% xlabel('Distance from terminus (km)');
% ylabel('Surface Elevation (m)')
% title('Glacier Surface')
% %print('GlacierScatter','-dpng','-r700')
% 
% plot(shmip1.md.mesh.x/1e3,shmip1.md.geometry.surface)
% xlabel('Distance from terminus (km)');
% ylabel('Surface Elevation (m)')
% title('Glacier Surface')
% print('GlacierPlot','-dpng','-r700')



plotmodel(shmip1.md, 'data', shmip1.md.results.TransientSolution(end).EffectivePressure/1e6, 'xlabel', 'Distance from terminus (m)', ...
    'ylabel', 'Distance across glacier', 'colorbartitle', 'N (MPa)')
print('IncorrectD2','-dpng','-r600')
% hold on;
% plotmodel(shmip1.md, 'data', shmip1.md.results.TransientSolution(end).EffectivePressure)

%subplot(3,2,2)
%plotmodel(shmip2.md, 'data', shmip1.md.results.TransientSolution(end).EffectivePressure)


% plot(tt, max(Q_1, [], 1));
% hold on;
% plot(tt, max(Q_4, [], 1));
% xlabel('Max Q');
% ylabel('Time')
% title('Max Q vs to for D1');
% legend({'ResTol = 1e-3', 'ResTol = 1e-4'}, 'Location', 'southeast');ChannelDischarge];
% 






