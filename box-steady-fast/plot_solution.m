% Run plot_solution after runme.m

figure('Units', 'inches', 'Position', [2, 2, 10, 5])
plotmodel(md,'data',md.results.TransientSolution(1).HydraulicPotential,'title','Initial hydraulic potential [Pa]',...
    'data',md.results.TransientSolution(end).HydraulicPotential,'title','Final hydraulic potential [Pa]',...
    'data',md.results.TransientSolution(1).HydrologySheetThickness,'title','Initial sheet thickness [m]',...
    'data',md.results.TransientSolution(end).HydrologySheetThickness,'title','Final sheet thickness [m]')