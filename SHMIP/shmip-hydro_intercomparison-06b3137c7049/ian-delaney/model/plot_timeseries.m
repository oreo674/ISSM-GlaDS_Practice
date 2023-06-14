[ax,h1,h2]=plotyy(t/para.year, phi(mean(199),:),t/para.year,tau(mean(119),:));
%plot(, 'r');
%ylim([0 40])
% set(ax(2),'YLim',[0 50])
% set(ax(1),'YLim',[0 2e6])
xlabel('Years')
legend('\phi', 'tau');