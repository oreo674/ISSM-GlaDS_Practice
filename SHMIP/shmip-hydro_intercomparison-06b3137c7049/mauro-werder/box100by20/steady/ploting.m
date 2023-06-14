channels = true;
meshnr = 3;

% $$$ ch = {};
% $$$ for scenario = 1:6
% $$$     ch{end+1} = load(['scenario1/ch',num2str(scenario),'.mat']);
% $$$ 
% $$$ end
% $$$ nc = {};
% $$$ for scenario = 1:6
% $$$     nc{end+1} = load(['scenario1/noch',num2str(scenario) ,'.mat']);
% $$$ end
% $$$ % 

aNnc_o = animate_something_and_channels();

figure
papx = 30;
papy = 20;
fontsize = 12;
plot_settings
for scenario = 1:6
    ax = subplot(3,2,scenario)
    aNnc_o.ax = ax
    aNnc_o.start = 21;
    [para, fields, phis, h_sheets, S_channels, u_beds] = unwrap_model_output_unc(ch{scenario});
    aNnc
    title(['Input ', num2str(para.input.source_term_s(0,0)*para.physical.day*1e3), '(mm/day)'])
end
plot_settings
print('-dpng', '-r300', fign)



