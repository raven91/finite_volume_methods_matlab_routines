function [] = ErrorConvergenceForLocalizedChimera3DVsPhi

    clc;
    clearvars;
    
    n = 40; dx = 1 / n;
    m = 40; dy = 1 / m;
    dt_recip = 1; dt = 1 / dt_recip;
    t_0 = 1 * dt;
    t_1 = t_0;
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/TestCaseLocalizedChimera3D/';
    size_of_float = 8;
    
    % num_cells l_1 l_2 l_\infty
    errors = dlmread(strcat(folder, 'error_convergence_vs_phi.txt'));
%     errors = errors(1 : end - 1, :);
    dxs = 1 ./ errors(:, 1);
    dys = 1 ./ errors(:, 2);
    dphis = 2 * pi ./ errors(:, 3);
    
    figure;
%     h_l1 = plot(dphis, errors(:, 2), '-*', 'LineWidth', 2.5); hold on;
%     h_l2 = plot(dphis, errors(:, 3), '-*', 'LineWidth', 2.5);
%     h_linf = plot(dphis, errors(:, 4), '-*', 'LineWidth', 2.5);
    h_l1 = loglog(dphis, errors(:, 4), '-*', 'LineWidth', 2.5); hold on;
    h_l2 = loglog(dphis, errors(:, 5), '-*', 'LineWidth', 2.5);
    h_linf = loglog(dphis, errors(:, 6), '-*', 'LineWidth', 2.5);
    x = linspace(dphis(end - 1, 1), dphis(1, 1));
    loglog(x, 10^(1.5) * x, '--k', 'LineWidth', 2.5);
    loglog(x, 10^(-0.4) * x .^ 2, '--k', 'LineWidth', 2.5);
%     loglog(x, 10^(2) * x .^ 4, '--k', 'LineWidth', 2.5);
	xlim([0.01 0.22]);
    ylim([10^(-5) 10^(1)]);
    legend([h_l1 h_l2 h_linf], {'L^1 error', 'L^2 error', 'L^\infty error'}, 'interpreter', 'tex', 'location', 'southeast');
    grid on;
    box on;
    
    xticks([10^(-2), 10^(-1)]);
    yticks([10^(-5), 10^(-3), 10^(-1), 10^(1)]);
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    
    error_substr = '<';
    for d = 1 : size(errors, 1) - 1
        error_substr = strcat(error_substr, '%f,');
    end % i
    error_substr = strcat(error_substr, '>');
    sprintf(strcat('L1: ', error_substr, ' |= %f'), diff(log(errors(:, 4))) ./ diff(log(dphis)),...
        (log(errors(end, 4)) - log(errors(1, 4))) / (log(dphis(end, 1)) - log(dphis(1, 1))))
    sprintf(strcat('L2: ', error_substr, ' |= %f'), diff(log(errors(:, 5))) ./ diff(log(dphis)),...
        (log(errors(end, 5)) - log(errors(1, 5))) / (log(dphis(end, 1)) - log(dphis(1, 1))))
    sprintf(strcat('L\\infty: ', error_substr, ' |= %f'), diff(log(errors(:, 6))) ./ diff(log(dphis)),...
        (log(errors(end, 6)) - log(errors(1, 6))) / (log(dphis(end, 1)) - log(dphis(1, 1))))
    
end