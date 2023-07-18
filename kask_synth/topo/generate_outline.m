function generate_outline(fname)
% PARAMETERS

% Set x_max, must match synthetic_topo.m
x_max = 100e3;
x_upper = 30e3;

% dt controls step size along boundary
dt = 0.01;
t = 0:dt:1;

% Create upper, lower boundaries and merge
x_upper = x_upper/2*(1 - cos(pi*t));
y_upper = synthetic_topo(x_upper);
y_upper(1) = 0; % Force outlet at (0, 0)
x_lower = fliplr(x_upper(2:end-1));
y_lower = fliplr(-y_upper(2:end-1));
x_full = [x_upper, x_lower];
y_full = [y_upper, y_lower];
n_points = length(x_full);

% Write .exp file
fid = fopen(fname, 'w');
fprintf(fid, '## Name:DomainOutline\n');
fprintf(fid, '## Icon:0\n');
fprintf(fid, '# Points Count  Value\n');
fprintf(fid, '%d %.6f\n', n_points, 1.0);
fprintf(fid, '# X pos Y pos\n');
for ii=1:n_points
    fprintf(fid, '%.2f %.2f\n', x_full(ii), y_full(ii));
end
fclose(fid);

% PLOT OUTLINE
figure;
dx = x_upper(2:end) - x_upper(1:end-1);
plot(0.5*(x_upper(2:end) + x_upper(1:end-1)), dx);
title('\Delta x')
grid()
xlabel('x (m)')
ylabel('\Delta x(m)')

f = figure;
plot(x_full, y_full, 'Marker', 'x')
title('Outline')
xlabel('x (m)')
ylabel('y (m)')
grid()
print(f, 'outline', '-dpng', '-r600')