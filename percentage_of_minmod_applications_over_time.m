n = 40;
m = 40;
l = 256;
folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/';
file_name = strcat('flux_limiter_dt_0.005_sigma_1_rho_0.3_alpha_1.45_Dphi_0.0075_', ...
        num2str(n), '_', num2str(m), '_', num2str(l), '_1000.txt');
x = dlmread(strcat(folder, file_name));

figure;
hold on;
max_applications = 1;%n * m * l;
h_x = plot(x(:, 1), x(:, 5) / max_applications, '-', 'LineWidth', 2.5);
h_y = plot(x(:, 1), x(:, 6) / max_applications, '-', 'LineWidth', 2.5);
h_phi = plot(x(:, 1), x(:, 7) / max_applications, '-', 'LineWidth', 2.5);
grid on;
box on;
% ylim([0 1]);
% ylabel('# of flux limiter applications');
ylabel('max absolute velocity');
xlabel('t');
legend([h_x, h_y, h_phi], {'$x$', '$y$', '$\varphi$'}, 'interpreter', 'latex');

set(gca,...
    'Units', 'normalized',...
    'FontUnits', 'points',...
    'FontWeight', 'normal',...
    'FontSize', 24,...
    'FontName', 'Helvetica',...
    'linew', 1);