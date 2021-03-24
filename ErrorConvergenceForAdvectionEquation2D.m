function [] = ErrorConvergenceForAdvectionEquation2D

    clc;
    clearvars;
    
    l = 1; dphi = 2 * pi / l;
    dt_recip = 1; dt = 1 / dt_recip;
    t_0 = 0; 
    t_1 = t_0;
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/TestCaseAdvectionEquation2D/';
    size_of_double = 8; 
    
    %% Show exact solution
    
    for t_0 = 0 : 0.1 : 0
        n = 100; dx = 1 / n;
        m = 100; dy = 1 / m;
        density = zeros(n, m);
        for i = 0 : n - 1
            for j = 0 : m - 1
                density(i + 1, j + 1) = ExactSolution(i * dx, j * dy, t_0);
            end % j
        end % i
        figure(1);
        surf(0 : dx : (n - 1) * dx, 0 : dy : (m - 1) * dy, density.', 'edgecolor', 'none');
        axis square;
        title(sprintf('t=%f', t_0));
        drawnow;
    end % t_0
    
    %% Quantify errors
    
    % num_cells_x num_cells_y l_1 l_2 l_\infty
    errors = dlmread(strcat(folder, 'error_convergence.txt'));
%     errors = errors(1 : end - 1, :);
    dxs = 1 ./ errors(:, 1);
    dys = 1 ./ errors(:, 2);
    drs = sqrt(dxs .^ 2 + dys .^ 2)
    
    figure(2);
%     h_l1 = plot(dphis, errors(:, 2), '-*', 'LineWidth', 2.5); hold on;
%     h_l2 = plot(dphis, errors(:, 3), '-*', 'LineWidth', 2.5);
%     h_linf = plot(dphis, errors(:, 4), '-*', 'LineWidth', 2.5);
    h_l1 = loglog(drs, errors(:, 3), '-*', 'LineWidth', 2.5); hold on;
    h_l2 = loglog(drs, errors(:, 4), '-*', 'LineWidth', 2.5);
    h_linf = loglog(drs, errors(:, 5), '-*', 'LineWidth', 2.5);
    x = linspace(10^(-3), 10^(-1));
%     loglog(x, x, '--k', 'LineWidth', 2.5);
    loglog(x, 10^(+3.4) * x .^ 2, '--k', 'LineWidth', 2.5);
%     loglog(x, 10^(2) * x .^ 4, '--k', 'LineWidth', 2.5);
    xlim([10^(-3), 10^(-1)]);
% 	ylim([10^(-10) 10^(-4)]);
    legend([h_l1 h_l2 h_linf], {'L^1 error', 'L^2 error', 'L^\infty error'}, 'interpreter', 'tex', 'location', 'southeast');
    grid on;
    box on;
    
    xticks([10^(-3), 10^(-2), 10^(-1)]);
    yticks([10^(-3), 10^(-1), 10^(+1)]);
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
    sprintf(strcat('L1: ', error_substr, ' |= %f'), diff(log(errors(:, 3))) ./ diff(log(drs)),...
        (log(errors(end, 3)) - log(errors(1, 3))) / (log(drs(end, 1)) - log(drs(1, 1))))
    sprintf(strcat('L2: ', error_substr, ' |= %f'), diff(log(errors(:, 4))) ./ diff(log(drs)),...
        (log(errors(end, 4)) - log(errors(1, 4))) / (log(drs(end, 1)) - log(drs(1, 1))))
    sprintf(strcat('L\\infty: ', error_substr, ' |= %f'), diff(log(errors(:, 5))) ./ diff(log(drs)),...
        (log(errors(end, 5)) - log(errors(1, 5))) / (log(drs(end, 1)) - log(drs(1, 1))))
    
end

function [f] = ExactSolution(x, y, t)

    k1 = 1.0; k2 = 1.0;
    mu1 = 0.0; mu2 = 0.0;
    f = exp(k1 * cos(2.0 * pi * (x - mu1 - t)) + k2 * cos(2.0 * pi * (y - mu2 - t))) / (besseli(0, k1) * besseli(0, k2));

end