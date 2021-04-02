function [] = PhaseTransitionForLocalizedChimeraVsDiffusion

	clc;
    clearvars;
    
    n = 40; m = 40; l = 256;
    dx = 1 / n; dy = 1 / m; dphi = 2 * pi / l;
    circular_mesh_x = zeros(1, n, 1);
    circular_mesh_x(1, :, 1) = 0 : dx : (n - 1) * dx;
    circular_mesh_x = repmat(circular_mesh_x, [l, 1, m]);
    circular_mesh_y = zeros(1, 1, m);
    circular_mesh_y(1, 1, :) = 0 : dy : (m - 1) * dy;
    circular_mesh_y = repmat(circular_mesh_y, [l, n, 1]);
    circular_mesh_phi = zeros(l, 1, 1);
    circular_mesh_phi(:, 1, 1) = 0 : dphi : (l - 1) * dphi;
    circular_mesh_phi = repmat(circular_mesh_phi, [1, n, m]);
    
    size_of_real = 8;
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/Continuation/ContinuationDiffusionFromNonhomogeneous/upper_path/';
    folder_info = dir(folder);
    is_simulation_file = contains({folder_info(:).name}, '.bin')';
    diffusion_coefficients = zeros(size(folder_info));
    spatial_order_parameters = zeros(size(folder_info));
    for i = 1 : size(folder_info, 1)
        if is_simulation_file(i)
            base_name = folder_info(i).name;
            key = '_Dphi_';
            key_index = strfind(base_name, key);
            diffusion_coefficients(i, 1) = sscanf(base_name(key_index+length(key) : end), '%g', 1);
            
            file_name = strcat(folder, base_name);
            file_info = dir(file_name);
            file_size = file_info.bytes;
            t_max = floor(file_size / ((1 + n * m * l) * size_of_real));
            file_id = fopen(file_name, 'r');
            try
                t_aver = 10;
                status = fseek(file_id, (1 + n * m * l) * (t_max - t_aver) * size_of_real, 'bof');
                assert(~status);
            catch
                t_aver = 1;
                status = fseek(file_id, (1 + n * m * l) * (t_max - t_aver) * size_of_real, 'bof');
                assert(~status);
            end
            
            spatial_order_parameters(i, 1) = 0;
            for tau = 1 : t_aver
                sprintf('t=%f', fread(file_id, 1, 'double'))
                continuum_limit_pdf = fread(file_id, [n * m * l, 1], 'double');
                spatial_order_parameters(i, 1) = spatial_order_parameters(i, 1) + AngularOrderParameter(continuum_limit_pdf, n, m, l, dx, dy, dphi, circular_mesh_phi) / t_aver;
            end % tau
        end
    end % i
    diffusion_coefficients = diffusion_coefficients(is_simulation_file);
    spatial_order_parameters = spatial_order_parameters(is_simulation_file);
    [diffusion_coefficients, sort_idx] = sort(diffusion_coefficients);
    spatial_order_parameters = spatial_order_parameters(sort_idx);

    figure;
    hold on;
    set(gca, 'ColorOrderIndex', 1);
    plot(diffusion_coefficients, spatial_order_parameters, '*-', 'LineWidth', 1, 'MarkerSize', 8);
%     xlim([0.007, 0.0205]);
    box on;
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    
    % Plot unstable branch of a homogeneous solution
%     diffusion = 0.00125 : 0.000625 : 0.0125;
%     plot(diffusion, diffusion.*0, '--k', 'LineWidth', 1);
%     ylim([-0.05, 0.5]);
%     LoadSolutionOfSelfconistentEquations();
%     ylim([0.68, 0.86]);
    % Plot an initial point of continuation method
%     order_parameter_value = 0.1114;%0.6896;%0.1114;
%     plot(0.0075, order_parameter_value, '*k', 'MarkerSize', 8, 'LineWidth', 1);
%     plot(0.0075, order_parameter_value, 'ok', 'MarkerSize', 8, 'LineWidth', 1);

end

function [res] = AngularOrderParameter(continuum_limit_pdf, n, m, l, dx, dy, dphi, circular_mesh_phi)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    circular_expectation_phi = sum(sum(sum(exp(1i * circular_mesh_phi) .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    res = abs(circular_expectation_phi);

end

function [res] = SpatialOrderParameter(continuum_limit_pdf, n, m, l, dx, dy, dphi, circular_mesh_x, circular_mesh_y)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    circular_expectation_xy = sum(sum(sum(exp(1i * 2 * pi * (circular_mesh_x + circular_mesh_y)) .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    res = abs(circular_expectation_xy);

end

function [res] = DensitySeparationParameter(continuum_limit_pdf, n, m, l, dphi)

   density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
   density_xy = squeeze(sum(density_all, 1)) * dphi;
   max_density = max(density_xy(:));
   min_density = min(density_xy(:));
%    res = (single(1) - single(min_density)) / (single(max_density) - single(min_density));
%    if isnan(res)
%        res = 1;
%    end
   res = abs(max_density - min_density);

end

function [res] = CircularDispersion(continuum_limit_pdf, n, m, l, dx, dy, dphi, circular_mesh_x, circular_mesh_y)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    circular_expectation_xy = sum(sum(sum(exp(1i * 2 * pi * (circular_mesh_x + circular_mesh_y)) .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    mean_resultant_length = abs(circular_expectation_xy);
    circular_variance = 1 - mean_resultant_length; % \in[0,1]
    circular_standard_deviation = sqrt(-2 * log(mean_resultant_length)); % \in[0,\infty)
    res = circular_variance;

end

function [res] = ToroidalOrderParameter(continuum_limit_pdf, n, m, l, dx, dy, dphi, circular_mesh_x, circular_mesh_y)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    % find mean values for x & y
    circular_expectation_x = sum(sum(sum(exp(1i * 2 * pi * circular_mesh_x) .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    circular_expectation_y = sum(sum(sum(exp(1i * 2 * pi * circular_mesh_y) .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    % find the mean vector of bivariate distribution of x & y using a torus
    R = 2; r = 1;
    x_1 = (R + r * cos(2 * pi * circular_mesh_x)) .* cos(2 * pi * circular_mesh_y);
    x_2 = (R + r * cos(2 * pi * circular_mesh_x)) .* sin(2 * pi * circular_mesh_y);
    x_3 = r * sin(2 * pi * circular_mesh_x);
    mu_1 = sum(sum(sum(x_1 .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    mu_2 = sum(sum(sum(x_2 .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    mu_3 = sum(sum(sum(x_3 .* density_all, 1) * dphi, 2) * dx, 3) * dy;
    % find where mean values of x & y lie on the same torus
    x_1 = (R + r * cos(angle(circular_expectation_x))) .* cos(angle(circular_expectation_y));
    x_2 = (R + r * cos(angle(circular_expectation_x))) .* sin(angle(circular_expectation_y));
    x_3 = r * sin(angle(circular_expectation_x));
    % compare how far a bivariate mean is from the torus
    res = norm([mu_1, mu_2, mu_3]) / norm([x_1, x_2, x_3]);

end

function [res] = HighDensityVolume(continuum_limit_pdf, n, m, l, dx, dy, dphi)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    density_xy = sum(density_all, 1) * dphi;
    homogeneous_density = 1;
    res = sum(density_xy(density_xy >= homogeneous_density)) * dx * dy;
    

end

function [res] = RadiusOfLocalizedCluster(continuum_limit_pdf, n, m, l, dx, dy, dphi)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    density_xy = sum(density_all, 1) * dphi;
    homogeneous_density = 1;
    high_density_area = sum(density_xy(:) >= homogeneous_density) * dx * dy;
    if abs(max(density_xy(:)) - min(density_xy(:))) < 1e-2
        res = 1;
    else
        res = sqrt(high_density_area / pi);
    end

end

function [res] = SecondaryAngularLayer(continuum_limit_pdf, n, m, l)

    density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    density_angularly_averaged = mean(density_all, 1);
    density_angularly_averaged = repmat(density_angularly_averaged, [l, 1, 1]);
    max_wrt_angle = max(abs(density_all - density_angularly_averaged), [], 1);
    res = min(max_wrt_angle(:)) / max(max_wrt_angle(:));

end

function [res] = CircularSkewness(continuum_limit_pdf, n, m, l, dx, dy, dphi, circular_mesh_phi)

	density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
    circular_expectation_phi = sum(sum(sum(density_all .* exp(1i * circular_mesh_phi), 1) * dphi, 2) * dx, 3) * dy;
    mean_resultant_length = abs(circular_expectation_phi);
    circular_variance = 1 - mean_resultant_length; % \in[0,1]
    mean_direction = angle(circular_expectation_phi);
    numerator = sum(sum(sum(density_all .* sin(2 * (circular_mesh_phi - mean_direction)), 1) * dphi, 2) * dx, 3) * dy;
    res = numerator / circular_variance ^ 1.5;

end

function [] = LoadSolutionOfSelfconistentEquations()

    folder = '/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/RegionForBifurcationAnalysis/';
    order_parameters = ReadOrderParameterFileCombined(folder);
    nondimensionalized_sigma = 4;
    diffusion_min = 0.005;
    diffusion_max = 0.0125 * nondimensionalized_sigma;
    phase_lag = 1.45;
    order_parameters = order_parameters(order_parameters(:, 3) == phase_lag, :);
    order_parameters = order_parameters(order_parameters(:, 2) >= diffusion_min & order_parameters(:, 2) <= diffusion_max, :);
    hold on;
    plot(order_parameters(:, 2) / nondimensionalized_sigma, order_parameters(:, 4), '--k', 'LineWidth', 1);

end

function [order_parameters] = ReadOrderParameterFileCombined(folder)

    % vec = [sigma, D_phi, alpha, ...]

    file_name_base = 'order_parameter_magnitude';
    file_extension = '.txt';
    order_parameters = [];
    for p = 0 : 105
        file_name = strcat(file_name_base, '_', num2str(p), file_extension);
        order_parameters = [order_parameters; dlmread(strcat(folder, file_name))];
    end % p

end