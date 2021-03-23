function [] = ErrorConvergenceForNonstationaryPhaseSynchronization1D

    clc;
    clearvars;
    
    n = 1; dx = 1 / n;
    m = 1; dy = 1 / m;
    dt_recip = 1; dt = 1 / dt_recip;
    t_0 = 0; 
    t_1 = t_0;
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/TestCaseNonstationaryPhaseSynchronization1D/';
    size_of_float = 4;
    
    %% Show exact solution
    
    for t_0 = 0 : 0.1 : t_0 + 0
        [grid_points, exact_density] = ExactSolution(t_0);
        figure(1); hold on;
        plot(grid_points, exact_density, 'LineWidth', 2.5, 'Color', [0 0 0]);
        title(sprintf('t=%f', t_0));
        drawnow;
    end % t_0
    
    %% Quantify errors
    
    % num_cells l_1 l_2 l_\infty
    errors = dlmread(strcat(folder, 'error_convergence.txt'));
%     errors = errors(1 : end - 1, :);
    dphis = 2 * pi ./ errors(:, 1);
    
    figure(2);
%     h_l1 = plot(dphis, errors(:, 2), '-*', 'LineWidth', 2.5); hold on;
%     h_l2 = plot(dphis, errors(:, 3), '-*', 'LineWidth', 2.5);
%     h_linf = plot(dphis, errors(:, 4), '-*', 'LineWidth', 2.5);
    h_l1 = loglog(dphis, errors(:, 2), '-*', 'LineWidth', 2.5); hold on;
    h_l2 = loglog(dphis, errors(:, 3), '-*', 'LineWidth', 2.5);
    h_linf = loglog(dphis, errors(:, 4), '-*', 'LineWidth', 2.5);
    x = linspace(dphis(end, 1), dphis(1, 1));
%     loglog(x, x, '--k', 'LineWidth', 2.5);
    loglog(x, 10^(-0.7) * x .^ 2, '--k', 'LineWidth', 2.5);
%     loglog(x, 10^(2) * x .^ 4, '--k', 'LineWidth', 2.5);
    xlim([0.005 0.22]);
	ylim([10^(-6) 10^(-1)]);
    legend([h_l1 h_l2 h_linf], {'L^1 error', 'L^2 error', 'L^\infty error'}, 'interpreter', 'tex', 'location', 'southeast');
    grid on;
    box on;
    
    xticks([10^(-2), 10^(-1)]);
    yticks([10^(-6), 10^(-4), 10^(-2)]);
    set(gca,...
        'Units', 'normalized',...
        'FontUnits', 'points',...
        'FontWeight', 'normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'linew', 1);
    
    error_substr = '<';
    for d = 1 : size(errors, 1) - 1
        error_substr = strcat(error_substr, '%f,');
    end % i
    error_substr = strcat(error_substr, '>');
    sprintf(strcat('L1: ', error_substr, ' |= %f'), diff(log(errors(:, 2))) ./ diff(log(dphis)),...
        (log(errors(end, 2)) - log(errors(1, 2))) / (log(dphis(end, 1)) - log(dphis(1, 1))))
    sprintf(strcat('L2: ', error_substr, ' |= %f'), diff(log(errors(:, 3))) ./ diff(log(dphis)),...
        (log(errors(end, 3)) - log(errors(1, 3))) / (log(dphis(end, 1)) - log(dphis(1, 1))))
    sprintf(strcat('L\\infty: ', error_substr, ' |= %f'), diff(log(errors(:, 4))) ./ diff(log(dphis)),...
        (log(errors(end, 4)) - log(errors(1, 4))) / (log(dphis(end, 1)) - log(dphis(1, 1))))
    
end

%% f is 1 x nw
function [grid_points, f] = ExactSolution(t)

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
    
    sigma = sym(1.0); D = sym(0.18); alpha = sym(0.9); op_magnitude = sym(0.66381359265845585); v = sym(-0.50563442409325432);
    
    E = @(x) exp((-v * (x - v * 0) + sigma * op_magnitude * cos(x - v * t + alpha)) / D);
    E_inverse = @(x) exp((v * (x - v * 0) - sigma * op_magnitude * cos(x - v * t + alpha)) / D);
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