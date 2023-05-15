% Run plot_solution after runme.m

figure('Units', 'inches', 'Position', [2, 2, 10, 5])
plotmodel(md,'data',md.results.TransientSolution(1).HydraulicPotential,'title','Initial hydraulic potential [Pa]',...
    'data',md.results.TransientSolution(end).HydraulicPotential,'title','Final hydraulic potential [Pa]',...
    'data',md.results.TransientSolution(1).HydrologySheetThickness,'title','Initial sheet thickness [m]',...
    'data',md.results.TransientSolution(end).HydrologySheetThickness,'title','Final sheet thickness [m]', ...
    'data',md.results.TransientSolution(1).EffectivePressure,'title','Initial N [Pa]',...
    'data',md.results.TransientSolution(end).EffectivePressure,'title','Final N [Pa]')

h_sheet = [md.results.TransientSolution.HydrologySheetThickness];
phi = [md.results.TransientSolution.HydraulicPotential];
Q = abs([md.results.TransientSolution.ChannelDischarge]);
S = [md.results.TransientSolution.ChannelArea];
tt = [md.results.TransientSolution.time];

figure
plot(tt, mean(h_sheet, 1))
title('h sheet')

figure
plot(tt, mean(phi, 1))
title('phi')

figure
plot(tt, max(Q, [], 1))
title('Q channel')

figure
plot(tt, sum(S, 1))