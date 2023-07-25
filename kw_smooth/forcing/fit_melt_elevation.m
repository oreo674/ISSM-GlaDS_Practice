% Fit melt-elevation curve from KR's melt data

% Read data and remove nan's
elev_raw = load('Elevation_KW.txt');
melt_raw = load('TotalAnnualMelt_KW_nodebris.txt');
mask = ~isnan(melt_raw.*elev_raw);
elev = elev_raw(mask);
melt = melt_raw(mask);

% Average by elevation
zbins = 800:50:3500;
zavg = 0.5*(zbins(2:end) + zbins(1:end-1));
mean_melt = zeros(size(zavg));
nbins = length(zbins) - 1;
for ii=1:nbins
    zmin = zbins(ii); zmax = zbins(ii+1);
    mean_melt(ii) = median(melt(elev>=zmin & elev<=zmax));
end

% Fit smoothing spline
options = fitoptions('Method', 'smoothingSpline', 'SmoothingParam', 1e-6);
[curve, gof, outputs] = fit(zavg', mean_melt', 'SmoothingSpline', options);
save('melt_elevation_fit.mat', 'curve')
curve_fit = curve(zavg);

% Plot it up
figure('Units', 'inches', 'Position', [2, 2, 4, 4])
hold on
scatter(melt, elev, 'Marker', '.')
plot(mean_melt, zavg, 'Color', 'red', 'LineWidth', 2)
plot(curve_fit, zavg, 'LineWidth', 2, 'Color', 'k')

grid on
xlabel('Melt (m w.e.)')
ylabel('Elevation (m asl.)')
legend({'Raw', 'Average', 'Spline'})
xlim([-0.25, 13])
print('melt_elevation_fit', '-dpng', '-r600')