% Sheet conductivity test
md_base = load('std_parameter.mat');
md_low = load('sheet_cond5e-4.mat')
md_high = load('sheet_cond5e-2.mat')
k_low = [md_low.md.results.TransientSolution.EffectivePressure];
k_base = [md_base.md.results.TransientSolution.EffectivePressure];
k_high = [md_high.md.results.TransientSolution.EffectivePressure];

%Graphing

scatter(md_low.md.mesh.x/1000, k_low(:,end));
hold on;
scatter(md_base.md.mesh.x/1000, k_base(:,end));
hold on;
scatter(md_high.md.mesh.x/1000, k_high(:,end));
xlabel('Distance from terminus (km)');
ylabel('Effective Pressure N');
legend({'k=5e-4','k=5e-3','k=5e-2'},'Location','southeast');
grid on;
%print('sheet_cond_comparison','-dpng','-r600');