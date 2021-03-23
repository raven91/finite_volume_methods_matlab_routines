function [] = ErrorConvergenceForStationaryPhaseSynchronization1D

    clc;
    clearvars;
    
    n = 1; dx = 1 / n;
    m = 1; dy = 1 / m;
    dt_recip = 1; dt = 1 / dt_recip;
    t_0 = 0; 
    t_1 = t_0;
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/TestCaseStationaryPhaseSynchronization1D/';
    size_of_float = 4;
    
    %% Show exact solution
    
    l = 1000; dphi = 2 * pi / l;
    exact_density = zeros(l, 1);
    for k = 0 : l - 1
        exact_density(k + 1, 1) = ExactSolution(k * dphi, t_0);
    end % k
    figure(1); hold on;
    plot(0 : dphi : (l - 1) * dphi, exact_density, 'LineWidth', 2.5, 'Color', [0 0 0]);
    title(sprintf('t=%f', t_0));
    drawnow;
    
    %% Quantify errors
    
    % num_cells l_1 l_2 l_\infty
    errors = dlmread(strcat(folder, 'error_convergence.txt'));
%     errors = errors(1 : end - 1, :);
    dphis = 2 * pi ./ errors(:, 1)
    
    figure(2);
%     h_l1 = plot(dphis, errors(:, 2), '-*', 'LineWidth', 2.5); hold on;
%     h_l2 = plot(dphis, errors(:, 3), '-*', 'LineWidth', 2.5);
%     h_linf = plot(dphis, errors(:, 4), '-*', 'LineWidth', 2.5);
    h_l1 = loglog(dphis, errors(:, 2), '-*', 'LineWidth', 2.5); hold on;
    h_l2 = loglog(dphis, errors(:, 3), '-*', 'LineWidth', 2.5);
    h_linf = loglog(dphis, errors(:, 4), '-*', 'LineWidth', 2.5);
    x = linspace(dphis(end, 1), dphis(1, 1));
%     loglog(x, x, '--k', 'LineWidth', 2.5);
    loglog(x, 10^(-0.1) * x .^ 2, '--k', 'LineWidth', 2.5);
%     loglog(x, 10^(2) * x .^ 4, '--k', 'LineWidth', 2.5);
	xlim([0.005 0.22]);
    ylim([10^(-7) 10^(-1)]);
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

%%
function [f] = ExactSolution(x, t)

    coupling_strength = 1.0;
    order_parameter = 0.94554218642329801;
    diffusion = 0.1;
    x_0 = 0.0;
    gamma = coupling_strength * order_parameter / diffusion;
    f = exp(gamma * cos(x - x_0)) / (2.0 * pi * besseli(0, gamma));

end