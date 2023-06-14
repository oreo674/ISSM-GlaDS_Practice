%% plot settings

% sets some stuff
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [papx papy]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [-papx*0.1 -papy*0.1 1.15*papx 1.15*papy]);
    
    tmp = findall(gcf, 'FontUnits','points');
    for oo=1:length(tmp)
        if get(tmp(oo),'FontSize')<fontsize
            set(tmp(oo),'FontSize',fontsize)
        end
    end
    set(gcf, 'Color', [1 1 1]);