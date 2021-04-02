function [] = SelfconsistentEquationsComparison

    %% order parameter vs diffusion for various \alpha

    clc;
    clearvars;
    p_eps = 0.00001; % to account for precision loss in output files
    figure; hold on; % R vs D_phi
    
    % results of the hydrodynamic theory
    phase_lags = [0.5; 1; 1.5];
    lower_vals_for_diffusion = [0.25; 0.1; 0.013]; % redefined because of the range of validity of hydrodynamic equations
    upper_vals_for_diffusion = [0.5; 0.3; 0.07];
    hydrodynamic_legend_handles = zeros(size(phase_lags));
    hydrodynamic_legend_names = cell(size(phase_lags));
    for i = 1 : size(phase_lags, 1)
        phase_lag = phase_lags(i, 1);
%         critical_diffusion = 0.5 * cos(phase_lag);
%         diffusions = (0.5 * critical_diffusion : 0.0001 : 1.5 * critical_diffusion);
        diffusions = (lower_vals_for_diffusion(i, 1) : 0.0001 : upper_vals_for_diffusion(i, 1));
        order_parameters = sqrt((16 * diffusions .^ 2 + sin(phase_lag) ^ 2) ./ (4 * diffusions) .* (cos(phase_lag) - 2 * diffusions)); % -> from the hydrodynamic analysis
        hydrodynamic_legend_handles(i, 1) = plot(diffusions, order_parameters, '-k', 'LineWidth', 1);
        hydrodynamic_legend_names(i, 1) = {strcat('\alpha = ', num2str(phase_lag, '%.2g'), ', HT')};
    end % i
    
    % results of solving the system of self-consistent equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/ComparisonToFiniteVolumeMethod/';
    order_parameters = ReadOrderParameterFile(folder);
    order_parameters = sortrows(order_parameters, 2);
    alphas = [0.5; 1; 1.5];
    lower_vals_for_diffusion = [0.0025; 0.0025; 0.0025];
    upper_vals_for_diffusion = [0.5; 0.33; 0.09];
    legend_handles = zeros(size(alphas));
    legend_names = cell(size(alphas));
    for i_alpha = 1 : numel(alphas)
        phase_lag = alphas(i_alpha, 1);
        diff_low = lower_vals_for_diffusion(i_alpha, 1) - p_eps;
        diff_high = upper_vals_for_diffusion(i_alpha, 1) + p_eps;
        alpha_idx = (order_parameters(:, 3) == phase_lag);
        alpha_idx = alpha_idx & (order_parameters(:, 2) >= diff_low) & (order_parameters(:, 2) <= diff_high);
        legend_handles(i_alpha, 1) = plot(order_parameters(alpha_idx, 2), order_parameters(alpha_idx, 4), '.', 'LineWidth', 1, 'MarkerSize', 12);
        legend_names(i_alpha, 1) = {strcat('\alpha = ', num2str(phase_lag, '%.2g'), ', SCE')};
    end % i_sigma
    
    % retuls of integrating continuum limit partial differential equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/SelfconsistentEquations/';
    dts = [0.0001, 0.0001, 0.0001];
    subfolders = ["phaselag_0p5/"; "phaselag_1p0/"; "phaselag_1p5/"];
    phase_lags = [0.5; 1; 1.5];
    diffusions = {[0.0025 : 0.0025 : 0.5];
        [0.0025 : 0.0025 : 0.33];
        [0.0025 : 0.0025 : 0.09]};
    numerical_legend_handles = zeros(size(phase_lags));
    numerical_legend_names = cell(size(phase_lags));
    for i = 1 : size(phase_lags, 1)
        phase_lag = phase_lags(i, 1);
        order_parameters = zeros(1, numel(diffusions{i}));
        for j = 1 : numel(diffusions{i})
            diffusion = diffusions{i}(j);
            file_name = strcat(folder, subfolders(i), 'dt_', num2str(dts(i)), '_sigma_1_rho_0_alpha_', num2str(phase_lag), '_Dphi_', num2str(diffusion), '_1_1_256_1000.txt');
            numerical_scheme_output = dlmread(file_name);
            order_parameters(1, j) = numerical_scheme_output(end, 2);
        end % j
        numerical_legend_handles(i, 1) = plot(diffusions{i}, order_parameters, 'o', 'LineWidth', 1, 'MarkerSize', 8);
%         set(numerical_legend_handles(i, 1), 'MarkerFaceColor', get(numerical_legend_handles(i, 1), 'Color'));
        numerical_legend_names(i, 1) = {strcat('\alpha = ', num2str(phase_lag, '%.2g'), ', FVM')};
        hold on;
    end % i
    
    grid on; box on;
    xlim([0 0.5]);
    ylim([0 1]);
    h_legend = legend([legend_handles; numerical_legend_handles], ...
        [legend_names; numerical_legend_names],...
        'Interpreter', 'Tex', 'Orientation', 'Vertical', 'Location', 'NorthEast', 'NumColumns', 2);
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    set(h_legend, 'FontName', 'Helvetica', 'FontSize', 22, 'Color', 'W');
    set(gcf, 'Color', 'W');
    
    %% velocity vs diffusion for various \alpha
    
    clc;
    clearvars;
    p_eps = 0.00001; % to account for precision loss in output files
    figure; hold on; % v vs D_phi
    
    % results of solving the system of self-consistent equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/ComparisonToFiniteVolumeMethod/';
    velocities = ReadVelocityFile(folder);
    velocities = sortrows(velocities, 2);
    alphas = [0.5; 1; 1.5];
    lower_vals_for_diffusion = [0.0025; 0.0025; 0.0025];
    upper_vals_for_diffusion = [0.5; 0.33; 0.09];
    legend_handles = zeros(size(alphas));
    legend_names = cell(size(alphas));
    for i_alpha = 1 : numel(alphas)
        phase_lag = alphas(i_alpha, 1);
        diff_low = lower_vals_for_diffusion(i_alpha, 1) - p_eps;
        diff_high = upper_vals_for_diffusion(i_alpha, 1) + p_eps;
        alpha_idx = (velocities(:, 3) == phase_lag);
        alpha_idx = alpha_idx & (velocities(:, 2) >= diff_low) & (velocities(:, 2) <= diff_high);
        legend_handles(i_alpha, 1) = plot(velocities(alpha_idx, 2), velocities(alpha_idx, 4), '.', 'LineWidth', 1, 'MarkerSize', 12);
        legend_names(i_alpha, 1) = {strcat('\alpha = ', num2str(phase_lag, '%.2g'), ', SCE')};
    end % i_sigma
    
    % retuls of integrating continuum limit partial differential equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/SelfconsistentEquations/';
    dts = [0.0001, 0.0001, 0.0001];
    subfolders = ["phaselag_0p5/"; "phaselag_1p0/"; "phaselag_1p5/"];
    phase_lags = [0.5; 1; 1.5];
    diffusions = {[0.0025 : 0.0025 : 0.5];
        [0.0025 : 0.0025 : 0.33];
        [0.0025 : 0.0025 : 0.09]};
    numerical_legend_handles = zeros(size(phase_lags));
    numerical_legend_names = cell(size(phase_lags));
    for i = 1 : size(phase_lags, 1)
        phase_lag = phase_lags(i, 1);
        velocities = zeros(1, numel(diffusions{i}));
        for j = 1 : numel(diffusions{i})
            diffusion = diffusions{i}(j);
            file_name = strcat(folder, subfolders(i), 'dt_', num2str(dts(i)), '_sigma_1_rho_0_alpha_', num2str(phase_lag), '_Dphi_', num2str(diffusion), '_1_1_256_1000.txt');
            numerical_scheme_output = dlmread(file_name);
            order_parameter_phases = mod(numerical_scheme_output(:, 3), 2 * pi);
            order_parameter_phases_velocity = ClassAEffectiveParticleDistance(order_parameter_phases(1 : end - 1), order_parameter_phases(2 : end), 2 * pi) / 1;
            velocities(1, j) = mean(order_parameter_phases_velocity(end - 100 : end, 1));
        end % j
        numerical_legend_handles(i, 1) = plot(diffusions{i}, velocities, 'o', 'LineWidth', 1, 'MarkerSize', 8);
%         set(numerical_legend_handles(i, 1), 'MarkerFaceColor', get(numerical_legend_handles(i, 1), 'Color'));
        numerical_legend_names(i, 1) = {strcat('\alpha = ', num2str(phase_lag, '%.2g'), ', FVM')};
        hold on;
    end % i
    
    grid on; box on;
    xlim([0 0.5]);
    ylim([-1 0]);
    h_legend = legend([legend_handles; numerical_legend_handles], ...
        [legend_names; numerical_legend_names],...
        'Interpreter', 'Tex', 'Orientation', 'Vertical', 'Location', 'SouthEast', 'NumColumns', 2);
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    set(h_legend, 'FontName', 'Helvetica', 'FontSize', 22, 'Color', 'W');
    set(gcf, 'Color', 'W');

    
    %% order parameter vs phase lag for different diffusions
    
    clc;
    clearvars;
    p_eps = 0.00001; % to account for precision loss in output files
    figure; hold on; % R vs alpha
    
    % results of the hydrodynamic theory
    diffusions = [0.1; 0.25; 0.4];
    lower_vals_for_phase_lags = [0.9; 0.5; 0.1];
    upper_vals_for_phase_lags = [1.54; 1.19; 0.8];
    hydrodynamic_legend_handles = zeros(size(diffusions));
    hydrodynamic_legend_names = cell(size(diffusions));
    for i = 1 : size(diffusions, 1)
        diffusion = diffusions(i, 1);
%         critical_phase_lag = acos(2 * diffusion);
%         phase_lags = (critical_phase_lag - 0.1 : 0.0001 : critical_phase_lag + 0.01);
        phase_lags = (lower_vals_for_phase_lags(i, 1) : 0.0001 : upper_vals_for_phase_lags(i, 1));
        order_parameters = sqrt((16 * diffusion ^ 2 + sin(phase_lags) .^ 2) / (4 * diffusion) .* (cos(phase_lags) - 2 * diffusion)); % -> from the hydrodynamic analysis
        hydrodynamic_legend_handles(i, 1) = plot(phase_lags, order_parameters, '-k', 'LineWidth', 1);
        hydrodynamic_legend_names(i, 1) = {strcat('D_\phi = ', num2str(diffusion, '%.2g'), ', HT')};
    end % i
    
    % results of solving the system of self-consistent equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/ComparisonToFiniteVolumeMethod/';
    order_parameters = ReadOrderParameterFile(folder);
    order_parameters = sortrows(order_parameters, 2);
    Ds = [0.1; 0.25; 0.4];
    lower_vals_for_phase_lags = [0; 0; 0];
    upper_vals_for_phase_lags = [1.54; 1.19; 0.8];
    legend_handles = zeros(size(Ds));
    legend_names = cell(size(Ds));
%     set(gca, 'ColorOrder', flipud(cbrewer('qual', 'Paired', numel(alphas)+0, 'pchip')), 'NextPlot', 'ReplaceChildren');
    for i_D = 1 : numel(Ds)
        D = Ds(i_D, 1);
        alpha_low = lower_vals_for_phase_lags(i_D, 1) - p_eps;
        alpha_high = upper_vals_for_phase_lags(i_D, 1) + p_eps;
        D_idx = (order_parameters(:, 2) == D);
        D_idx = D_idx & (order_parameters(:, 3) >= alpha_low) & (order_parameters(:, 3) <= alpha_high);
        legend_handles(i_D, 1) = plot(order_parameters(D_idx, 3), order_parameters(D_idx, 4), '.', 'LineWidth', 1, 'MarkerSize', 12);
        legend_names(i_D, 1) = {strcat('D_\phi = ', num2str(D, '%.2g'), ', SCE')};
        hold on;
    end % i_D
    
    % retuls of integrating continuum limit partial differential equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/SelfconsistentEquations/';
    dts = [0.0001, 0.0001, 0.0001];
    subfolders = ["diffusion_0p1/"; "diffusion_0p25/"; "diffusion_0p4/"];
    diffusions = [0.1; 0.25; 0.4];
    phase_lags = {[0 : 0.01 : 1.54]; 
        [0 : 0.01 : 1.19];
        [0 : 0.01 : 0.8]};
    numerical_legend_handles = zeros(size(diffusions));
    numerical_legend_names = cell(size(diffusions));
    for i = 1 : size(diffusions, 1)
        diffusion = diffusions(i, 1);
        order_parameters = zeros(1, numel(phase_lags{i}));
        for j = 1 : numel(phase_lags{i})
            phase_lag = phase_lags{i}(j);
            file_name = strcat(folder, subfolders(i), 'dt_', num2str(dts(i)), '_sigma_1_rho_0_alpha_', num2str(phase_lag), '_Dphi_', num2str(diffusion), '_1_1_256_1000.txt');
            numerical_scheme_output = dlmread(file_name);
            order_parameters(1, j) = numerical_scheme_output(end, 2);
        end % j
        numerical_legend_handles(i, 1) = plot(phase_lags{i}, order_parameters, 'o', 'LineWidth', 1, 'MarkerSize', 8);
%         set(numerical_legend_handles(i, 1), 'MarkerFaceColor', get(numerical_legend_handles(i, 1), 'Color'));
        numerical_legend_names(i, 1) = {strcat('D_\phi = ', num2str(diffusion, '%.2g'), ', FVM')};
        hold on;
    end % i
    
    grid on; box on;
    xlim([0 pi/2]);
    ylim([0 1]);
    h_legend = legend([legend_handles; numerical_legend_handles],...
        [legend_names; numerical_legend_names],...
        'Interpreter', 'Tex', 'Orientation', 'Vertical', 'Location', 'NorthEast', 'NumColumns', 2);
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    set(h_legend, 'FontName', 'Helvetica', 'FontSize', 22, 'Color', 'W');
    set(gcf, 'Color', 'W');
    
    %% velocity vs phase lag for different diffusions
    
    clc;
    clearvars;
    p_eps = 0.00001; % to account for precision loss in output files
    figure; hold on; % v vs alpha
    
    % results of solving the system of self-consistent equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/ComparisonToFiniteVolumeMethod/';
    velocities = ReadVelocityFile(folder);
    velocities = sortrows(velocities, 2);
    Ds = [0.1; 0.25; 0.4];
    lower_vals_for_phase_lags = [0; 0; 0];
    upper_vals_for_phase_lags = [1.54; 1.19; 0.8];
    legend_handles = zeros(size(Ds));
    legend_names = cell(size(Ds));
%     set(gca, 'ColorOrder', flipud(cbrewer('qual', 'Paired', numel(alphas)+0, 'pchip')), 'NextPlot', 'ReplaceChildren');
    for i_D = 1 : numel(Ds)
        D = Ds(i_D, 1);
        alpha_low = lower_vals_for_phase_lags(i_D, 1) - p_eps;
        alpha_high = upper_vals_for_phase_lags(i_D, 1) + p_eps;
        D_idx = (velocities(:, 2) == D);
        D_idx = D_idx & (velocities(:, 3) >= alpha_low) & (velocities(:, 3) <= alpha_high);
        legend_handles(i_D, 1) = plot(velocities(D_idx, 3), velocities(D_idx, 4), '.', 'LineWidth', 1, 'MarkerSize', 12);
        legend_names(i_D, 1) = {strcat('D_\phi = ', num2str(D, '%.2g'), ', SCE')};
        hold on;
    end % i_D
    
    % retuls of integrating continuum limit partial differential equations
    set(gca, 'ColorOrderIndex', 1);
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/SelfconsistentEquations/';
    dts = [0.0001, 0.0001, 0.0001];
    subfolders = ["diffusion_0p1/"; "diffusion_0p25/"; "diffusion_0p4/"];
    diffusions = [0.1; 0.25; 0.4];
    phase_lags = {[0 : 0.01 : 1.54]; 
        [0 : 0.01 : 1.19];
        [0 : 0.01 : 0.8]};
    numerical_legend_handles = zeros(size(diffusions));
    numerical_legend_names = cell(size(diffusions));
    for i = 1 : size(diffusions, 1)
        diffusion = diffusions(i, 1);
        velocities = zeros(1, numel(phase_lags{i}));
        for j = 1 : numel(phase_lags{i})
            phase_lag = phase_lags{i}(j);
            file_name = strcat(folder, subfolders(i), 'dt_', num2str(dts(i)), '_sigma_1_rho_0_alpha_', num2str(phase_lag), '_Dphi_', num2str(diffusion), '_1_1_256_1000.txt');
            numerical_scheme_output = dlmread(file_name);
            order_parameter_phases = mod(numerical_scheme_output(:, 3), 2 * pi);
            order_parameter_phases_velocity = ClassAEffectiveParticleDistance(order_parameter_phases(1 : end - 1), order_parameter_phases(2 : end), 2 * pi) / 1;
            velocities(1, j) = mean(order_parameter_phases_velocity(end - 100 : end, 1));
        end % j
        numerical_legend_handles(i, 1) = plot(phase_lags{i}, velocities, 'o', 'LineWidth', 1, 'MarkerSize', 8);
%         set(numerical_legend_handles(i, 1), 'MarkerFaceColor', get(numerical_legend_handles(i, 1), 'Color'));
        numerical_legend_names(i, 1) = {strcat('D_\phi = ', num2str(diffusion, '%.2g'), ', FVM')};
        hold on;
    end % i
    
    grid on; box on;
    xlim([0 pi/2]);
    ylim([-1 0]);
    h_legend = legend([legend_handles; numerical_legend_handles],...
        [legend_names; numerical_legend_names],...
        'Interpreter', 'Tex', 'Orientation', 'Vertical', 'Location', 'SouthEast', 'NumColumns', 2);
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    set(h_legend, 'FontName', 'Helvetica', 'FontSize', 22, 'Color', 'W');
    set(gcf, 'Color', 'W');

    %% solutions of self-consistent equations
    
    clc;
    clearvars;
    folder = '/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/HomogeneousSolutionsNonzeroLag/max_digits10/10000/';
    order_parameters = ReadOrderParameterFileCombined(folder);
    velocities = ReadVelocityFileCombined(folder);
    
    % plot order parameter magnitude
    figure; hold on;
    sigma = 1;
    [mesh_D_phi, mesh_alpha, mesh_order_parameter] = ParameterVectorMesh(order_parameters, sigma);
    mesh_order_parameter = mesh_order_parameter(1 : end - 2, :);
    s = pcolor(mesh_D_phi, mesh_alpha, mesh_order_parameter);
    s.EdgeColor = 'flat';
    axis tight; box on;
    colormap(flipud(bone));
    colorbar;
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    sigma = 1;
    alpha = 0 : 0.01 : 1.57;
    D_phi = 0.5 * sigma * cos(alpha);
    plot(D_phi, alpha, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
    set(gca, 'ColorOrderIndex', 1);
    phase_lags = [0.5; 1; 1.5];
    diffusions = {[0.0025 : 0.0025 : 0.5];
        [0.0025 : 0.0025 : 0.33];
        [0.0025 : 0.0025 : 0.09]};
    for i = 1 : size(phase_lags, 1)
        phase_lag = phase_lags(i, 1);
        plot(diffusions{i}, ones(size(diffusions{i})) * phase_lag, '--', 'LineWidth', 2);
    end % i
    set(gca, 'ColorOrderIndex', 1);
    diffusions = [0.1; 0.25; 0.4];
    phase_lags = {[0 : 0.01 : 1.54]; 
        [0 : 0.01 : 1.19];
        [0 : 0.01 : 0.8]};
    for i = 1 : size(diffusions, 1)
        diffusion = diffusions(i, 1);
        plot(ones(size(phase_lags{i})) * diffusion, phase_lags{i}, '--', 'LineWidth', 2);
    end % i
    
    % plot group velocity
    figure; hold on;
    sigma = 1;
    [mesh_D_phi, mesh_alpha, mesh_velocity] = ParameterVectorMesh(velocities, sigma);
    mesh_velocity = mesh_velocity(1 : end - 2, :);
    s = pcolor(mesh_D_phi, mesh_alpha, mesh_velocity);
    s.EdgeColor = 'flat';
    axis tight;
    box on;
    colormap((bone));
    colorbar;
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    sigma = 1;
    alpha = 0 : 0.01 : 1.57;
    D_phi = 0.5 * sigma * cos(alpha);
    plot(D_phi, alpha, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
    set(gca, 'ColorOrderIndex', 1);
    phase_lags = [0.5; 1; 1.5];
    diffusions = {[0.0025 : 0.0025 : 0.5];
        [0.0025 : 0.0025 : 0.33];
        [0.0025 : 0.0025 : 0.09]};
    for i = 1 : size(phase_lags, 1)
        phase_lag = phase_lags(i, 1);
        plot(diffusions{i}, ones(size(diffusions{i})) * phase_lag, '--', 'LineWidth', 2);
    end % i
    set(gca, 'ColorOrderIndex', 1);
    diffusions = [0.1; 0.25; 0.4];
    phase_lags = {[0 : 0.01 : 1.54]; 
        [0 : 0.01 : 1.19];
        [0 : 0.01 : 0.8]};
    for i = 1 : size(diffusions, 1)
        diffusion = diffusions(i, 1);
        plot(ones(size(phase_lags{i})) * diffusion, phase_lags{i}, '--', 'LineWidth', 2);
    end % i
    
end

function [order_parameters] = ReadOrderParameterFile(folder)

    file_name = 'order_parameter_magnitude.txt';
    order_parameters = dlmread(strcat(folder, file_name));

end

function [velocities] = ReadVelocityFile(folder)

    file_name = 'velocity.txt';
    velocities = dlmread(strcat(folder, file_name));

end

function [order_parameters] = ReadOrderParameterFileCombined(folder)

    % vec = [sigma, D_phi, alpha, ...]

    file_name_base = 'order_parameter_magnitude';
    file_extension = '.txt';
    order_parameters = [];
    for p = 0 : 19
        file_name = strcat(file_name_base, '_nq', num2str(p), file_extension);
        order_parameters = [order_parameters; dlmread(strcat(folder, file_name))];
    end % p

end

function [velocities] = ReadVelocityFileCombined(folder)

    % vec = [sigma, D_phi, alpha, ...]

    file_name_base = 'velocity';
    file_extension = '.txt';
    velocities = [];
    for p = 0 : 19
        file_name = strcat(file_name_base, '_nq', num2str(p), file_extension);
        velocities = [velocities; dlmread(strcat(folder, file_name))];
    end % p

end

function [mesh_D_phi, mesh_alpha, mesh_vec] = ParameterVectorMesh(vec, sigma)

    % vec = [sigma, D_phi, alpha, ...]
    
    idx = (vec(:, 1) == sigma);
    subvec = vec(idx, :);
    
    D_phi_min = min(subvec(:, 2), [], 1);
    D_phi_max = max(subvec(:, 2), [], 1);
    n_D_phi = numel(unique(subvec(:, 2)));
    d_D_phi = diff(sort(unique(subvec(:, 2)), 1), 1, 1); d_D_phi = d_D_phi(1, 1);
    D_phis = (D_phi_min : d_D_phi : D_phi_max);
    
    alpha_min = min(subvec(:, 3), [], 1);
    alpha_max = max(subvec(:, 3), [], 1);
    n_alpha = numel(unique(subvec(:, 3)));
    d_alpha = diff(sort(unique(subvec(:, 3)), 1), 1, 1); d_alpha = d_alpha(1, 1);
    alphas = (alpha_min : d_alpha : alpha_max).';
    if (numel(alphas) < n_alpha)
        alphas = [alphas; alpha_max];
    end
    
    % D_phi - x axis
    % alpha - y axis
    mesh_D_phi = repmat(D_phis, [numel(alphas), 1]) - d_D_phi;
    mesh_alpha = repmat(alphas, [1, numel(D_phis)]) - d_alpha;
    mesh_vec = zeros(n_alpha, n_D_phi);
    for i = 1 : size(subvec, 1)
        D_phi = subvec(i, 2);
        idx_D_phi = round((D_phi - D_phi_min) / d_D_phi);
        alpha = subvec(i, 3);
        idx_alpha = round((alpha - alpha_min) / d_alpha);
        mesh_vec(idx_alpha + 1, idx_D_phi + 1) = subvec(i, end);
    end % i

end

function [dx] = ClassAEffectiveParticleDistance(x_i, x_j, period)
    
    reciprocal_period = 1.0 / period;
    dx = x_j - x_i;
    dx = dx - fix(dx * 2.0 * reciprocal_period) * period;

end