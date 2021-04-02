function [] = SolutionComparisonForNonstationaryPhaseSynchronization1D

    clc;
    clearvars;
    
    % Show exact solution
    figure; hold on;
    parameters = {[1.0, 0.1, 1.0, 0.802199127693474, -0.63583668406608];
        [1.0, 0.15, 1.0, 0.66281545417635168, -0.55291509351945056];
        [1.0, 0.2, 1.0, 0.49900250126841533, -0.48775623710338245];
        [1.0, 0.25, 1.0, 0.26487835523831732, -0.43737613627975269]};
    for i = 1 : numel(parameters)
        sigma = parameters{i}(1, 1);
        D = parameters{i}(1, 2);
        alpha = parameters{i}(1, 3);
        op_magnitude = parameters{i}(1, 4);
        v = parameters{i}(1, 5);
        [grid_points, exact_density] = ExactSolution(sigma, D, alpha, op_magnitude, v);
        h_exact = plot(grid_points, exact_density, 'LineWidth', 1, 'Color', [0 0 0]);
    end % i
    set(gca, 'ColorOrderIndex', 1);

    % Show numerical solution
    n = 1;
    m = 1;
    l = 256;
    skip_time = 1000 * 1;
    dx = 1 / n;
    dy = 1 / m;
    dphi = 2 * pi / l;

    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/';
    % folder = '/Volumes/Kruk/spc2/spc2FiniteVolumeMethods/';
    numerical_legend_handles = zeros(numel(parameters), 1);
    numerical_legend_names = cell(numel(parameters), 1);
    markers = ['o', 'o', 'o', 'o'];
    for i = 1 : numel(parameters)
        sigma = parameters{i}(1, 1);
        D = parameters{i}(1, 2);
        alpha = parameters{i}(1, 3);
        base_name = strcat('dt_1e-05_sigma_', num2str(sigma), '_rho_1_alpha_', num2str(alpha), '_Dphi_', num2str(D), '_', ...
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
        numerical_legend_handles(i, 1) = plot(mod(discretization_grid - mean_direction + pi, 2 * pi), density_phi, 'o', 'MarkerSize', 6);
        % numerical_legend_handles(i, 1) = plot(discretization_grid, circshift(density_phi, -mean_direction_idx + l / 2));
        numerical_legend_names{i, 1} = strcat('D_\phi=', num2str(D), ', FVM');
    end % i
    axis([0 2*pi 0 0.8]);
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

% f is 1 x nw
function [grid_points, f] = ExactSolution(sigma, D, alpha, op_magnitude, v)

    digits(1500);
    nw = 100;
    w_min = 0;
    w_max = 2 * pi;
    dw = (w_max - w_min) / nw;
    grid_points = w_min + (0 : 1 : nw - 0) * dw;
    grid_point_weight = grid_points .* 0 + dw;
%     [grid_point, grid_point_weight] = GenerateGaussLegendrePoints(w_min, w_max, nw);
%     grid_point = w_min + (0 : 1 : nw) * dw;
%     grid_point_weight = grid_point .* 0 + dw / 2;
    
%     sigma = sym(1.0); D = sym(0.1); alpha = sym(1.0); op_magnitude = sym(0.802199127693474); v = sym(-0.63583668406608);
    x_0 = pi;
    
    E = @(x) exp((-v * (x - v * 0 - x_0) + sigma * op_magnitude * cos(x - x_0 + alpha)) / D);
    E_inverse = @(x) exp((v * (x - v * 0 - x_0) - sigma * op_magnitude * cos(x - x_0 + alpha)) / D);
    [IntegralUpToCurrentPoint, I_2pi] = IntegrateRange(E_inverse, 0 + v * 0, 2 * pi + v * 0, nw);
    c_e = exp(2.0 * sym(pi) * v / D);
    
    jw = (0 : nw - 1);
    f_norm = vpa(sum(grid_point_weight(1, jw + 1) .* E(grid_points(1, jw + 1)) .* (1.0 + (c_e - 1.0) * IntegralUpToCurrentPoint(1, jw + 1) / I_2pi), 2));
    f = vpa(E(grid_points(1, jw + 1)) .* (1.0 + (c_e - 1.0) * IntegralUpToCurrentPoint(1, jw + 1) / I_2pi) / f_norm);
    grid_points = grid_points(1, jw + 1);

end

function [IntegralUpToCurrentPoint, I_2pi] = IntegrateRange(f, lower_limit, upper_limit, n_points)

    dx = (upper_limit - lower_limit) / n_points;
%     grid_point = lower_limit + (0 : 1 : n_points - 0) * dx; % rectangle method
%     grid_point_weight = grid_point .* 0 + dx;
%     [grid_point, grid_point_weight] = GenerateGaussLegendrePoints(lower_limit, upper_limit, n_points);
    grid_point = lower_limit + (0 : 1 : n_points - 1) * dx; % trapezoidal rule
    grid_point_weight = grid_point .* 0 + dx / 2;
    
    IntegralUpToCurrentPoint = sym(zeros(1, n_points));
%     IntegralUpToCurrentPoint(1, 2 : end) = vpa(grid_point_weight(1, (0 : n_points - 1) + 1) .* f(grid_point(1, (0 : n_points - 1) + 1)));
    IntegralUpToCurrentPoint(1, 2 : end) = vpa(grid_point_weight(1, (0 : n_points - 2) + 1) .* f(grid_point(1, (0 : n_points - 2) + 1)) ...
        + grid_point_weight(1, (1 : n_points - 1) + 1) .* f(grid_point(1, (1 : n_points - 1) + 1)));
    IntegralUpToCurrentPoint = cumsum(IntegralUpToCurrentPoint, 2);
    
    I_2pi = IntegralUpToCurrentPoint(1, end) + vpa((f(grid_point(1, end)) + f(upper_limit)) * 0.5 * (upper_limit - grid_point(1, end)));

end