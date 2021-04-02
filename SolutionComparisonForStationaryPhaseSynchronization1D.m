function [] = SolutionComparisonForStationaryPhaseSynchronization1D

    clc;
    clearvars;
    
    % Show exact solution
    
    figure; hold on;
    parameters = {[1.0, 0.1, 0.94554218642329801]; [1.0, 0.2, 0.8768234220629445]; [1.0, 0.3, 0.77312962537029162]; [1.0, 0.4, 0.5897079877827418]};
    for i = 1 : numel(parameters)
        coupling_strength = parameters{i}(1, 1);
        diffusion = parameters{i}(1, 2);
        order_parameter = parameters{i}(1, 3);
        [x, f] = ExactSolution(coupling_strength, diffusion, order_parameter);
        h_exact = plot(x, f, 'LineWidth', 1, 'Color', [0 0 0]);
    end % i
    set(gca, 'ColorOrderIndex', 1);

    % Show numerical solution
    n = 1;
    m = 1;
    l = 256;
    skip_time = 200 * 1;
    dx = 1 / n;
    dy = 1 / m;
    dphi = 2 * pi / l;

    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/TestCaseStationaryPhaseSynchronization1D/';
    % folder = '/Volumes/Kruk/spc2/spc2FiniteVolumeMethods/';
    numerical_legend_handles = zeros(numel(parameters), 1);
    numerical_legend_names = cell(numel(parameters), 1);
    markers = ['o', 'o', 'o', 'o'];
    for i = 1 : numel(parameters)
        coupling_strength = parameters{i}(1, 1);
        diffusion = parameters{i}(1, 2);
        base_name = strcat('dt_1e-05_sigma_', num2str(coupling_strength), '_rho_1_alpha_0_Dphi_', num2str(diffusion), '_', ...
            num2str(n), '_', num2str(m), '_', num2str(l), '_1000');
        file_name = strcat(folder, base_name, '.bin');
        file_id = fopen(file_name);
        size_of_float = 8;
        status = fseek(file_id, (1 + n * m * l) * skip_time * size_of_float, 'bof'); % skip t_0 steps
        assert(~status);

        continuum_limit_pdf = fread(file_id, [1 + n * m * l, 1], 'double');
        time = continuum_limit_pdf(1, 1);
        density_all = reshape(continuum_limit_pdf(2 : n * m * l + 1, 1), l, n, m);
        density_phi = squeeze(sum(sum(density_all, 2), 3) * dx * dy);
        % density_phi = [density_phi; density_phi
        discretization_grid = 0 : dphi : (l - 1) * dphi;
        mean_direction = mod(angle(mean(exp(1i * discretization_grid) .* density_phi.')), 2 * pi);
        mean_direction_idx = floor(mean_direction / (2 * pi) * l);

        % figure;
        numerical_legend_handles(i, 1) = plot(mod(discretization_grid - mean_direction + pi, 2 * pi), density_phi, markers(i), 'MarkerSize', 6);
        % numerical_legend_handles(i, 1) = plot(discretization_grid, circshift(density_phi, -mean_direction_idx + l / 2));
        numerical_legend_names{i, 1} = strcat('D_\phi=', num2str(diffusion), ', FVM');
    end % i
    axis([0 2*pi 0 1.3]);
    grid on; box on; hold on;
    h_legend = legend([numerical_legend_handles; h_exact], [numerical_legend_names; {'Exact solutions'}], 'FontSize', 20, 'Interpreter', 'Tex', 'Location', 'NorthEast');
    set(gca,...
        'Units', 'normalized',...
        'FontUnits', 'points',...
        'FontWeight', 'normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'linew', 1);
%     set(h_legend, 'Color', 'None');
    
end

function [x, f] = ExactSolution(coupling_strength, diffusion, order_parameter)

    l = 1000; dphi = 2 * pi / l;
    gamma = coupling_strength * order_parameter / diffusion;
    x_0 = pi;
    x = 0 : dphi : (l - 1) * dphi;
    f = exp(gamma * cos(x - x_0)) / (2.0 * pi * besseli(0, gamma));

end