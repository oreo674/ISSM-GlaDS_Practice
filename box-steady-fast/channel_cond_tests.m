% Importing models and N for each k_c
md_low = load('5e-3chancond_1e-3sens.mat');
md_high = load('0,5chancond_1e-2sens.mat');
md_base = load('std_parameter');
kc_low = [md_low.md.results.TransientSolution.EffectivePressure];
kc_high=[md_high.md.results.TransientSolution.EffectivePressure];
kc_base = [md_base.md.results.TransientSolution.EffectivePressure];

%Building scatter plot
scatter(md_low.md.mesh.x/1000, kc_low(:,end));
hold on;
scatter(md_base.md.mesh.x/1000, kc_base(:,end));
hold on;
scatter(md_high.md.mesh.x/1000, kc_high(:,end));
xlabel('km');
ylabel('N');
legend({'k_c=0.005','k_c=0.05','k_c=0.5'},'Location','southeast');
grid on;
%print('channel_cond_comparison','-dpng','-r600');
