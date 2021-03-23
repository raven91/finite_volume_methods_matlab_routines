function [] = fvmContinuumLimitPdeMarginalDistributions

    clc;
    clearvars;
    
    n = 40; m = 40; l = 256;
    dx = 1 / n; dy = 1 / m; dphi = 2 * pi / l;
    
    folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/';
%     folder = '/Volumes/Kruk/spc2/spc2FiniteVolumeMethods/SpatiallyNonhomogeneousDynamics/';
    base_name = strcat('dt_0.005_sigma_1_rho_0.3_alpha_1.51_Dphi_0.0075_', ...
        num2str(n), '_', num2str(m), '_', num2str(l), '_1000');
    file_name = strcat(folder, base_name, '.bin');
    size_of_real = 8;
    file_info = dir(file_name);
    file_size = file_info.bytes;
    t_max = 1;%floor(file_size / ((1 + n * m * l) * size_of_real))
    file_id = fopen(file_name);
    status = fseek(file_id, (1 + n * m * l) * (t_max - 1) * size_of_real, 'bof'); % skip t_0 steps
    assert(~status);
    
    movie_folder = '/Users/nikita/Documents/Projects/spc2/VideoStorage/';
    fname_movie = strcat(movie_folder, base_name, '.mp4');
    mov = VideoWriter(fname_movie, 'MPEG-4');
    mov.FrameRate = 50;
    mov.Quality = 20;
    open(mov);
    
    for t = 0 : 0
        integration_time = fread(file_id, 1, 'double')
        continuum_limit_pdf = fread(file_id, [n * m * l, 1], 'double');
        density_all = reshape(continuum_limit_pdf(1 : n * m * l, 1), l, n, m);
        
        PlotAngularDensity(density_all, l, dx, dy, dphi);
%         pause(1);
%         drawnow;
%         frame = getframe(gcf);
%         writeVideo(mov, frame);
%         [M, I] = max(density_phi, [], 1);
        
        PlotOneDimSpatialDensity(density_all, n, m, dx, dy, dphi);
%         drawnow;
        
        PlotTwoDimSpatialDensity(density_all, n, m, dx, dy, dphi);
%         time_title = sprintf('t=%.2f\n', integration_time);
%         title(time_title);
%         drawnow;
%         frame = getframe(gcf);
%         writeVideo(mov, frame);
        
%         PlotTwoDimMomentumField(density_all, n, m, l, dx, dy, dphi);

%         PlotSpatialDensityProfileInTheDirectionOfMotion(density_all, n, m, l, dx, dy, dphi);
    end % t
    fclose(file_id);
    
end

function [density_phi] = PlotAngularDensity(density_all, l, dx, dy, dphi)

    
    density_phi = squeeze(sum(sum(density_all, 2), 3) * dx * dy);
%     density_phi = squeeze(density_all(:, 30, 30));
%     [~, index] = max(density_phi);
%     density_phi = circshift(density_phi, l / 2 - index);
%     density_phi = 1 + sin([0 : dphi : (l - 1) * dphi] - t_0); % analytic solution

        figure;
%         density_phi = density_phi / (sum(density_phi) * dphi);
%         semilogy(0 : dphi : (l - 1) * dphi, density_phi, 'LineWidth', 2.5, 'Color', [174 0 21] ./ 255);
%         set(gca, 'YTickLabel', []);
%         plot(0 : dphi : (l - 1) * dphi, density_phi, 'LineWidth', 2.5, 'Color', [174 0 21] ./ 255);
        [~, mi] = min(density_phi);
        plot((0 : dphi : (l - 1) * dphi), circshift(density_phi, -mi), 'k-', 'LineWidth', 4);
%         plot3(0 : dphi : (l - 1) * dphi, (0 : dphi : (l - 1) * dphi) * 0 + t, density_phi,...
%             'LineWidth', 2.5);
        grid on;
        box on;
%         xlim([0 2*pi]);
%         ylim([0 4]);
%         ytickformat('%.2f')
%         ylabel('f(\varphi)', 'Interpreter', 'Tex');
%     set(gca, 'YScale', 'log');
        set(gca,...
            'Units', 'Normalized',...
            'FontUnits', 'Points',...
            'FontWeight', 'Normal',...
            'FontSize', 30,...
            'FontName', 'Helvetica',...
            'LineW', 1);
%         set(gcf, 'Color', 'w');

end

function [density_x, density_y] = PlotOneDimSpatialDensity(density_all, n, m, dx, dy, dphi)

    figure;
    density_x = squeeze(sum(sum(density_all, 1), 3) * dphi * dy);
%         [~, index] = max(density_x);
%         density_x = circshift(density_x, n / 2 - index);
        h_x = plot(0 : dx : (n - 1) * dx, density_x / sum(density_x) / dx, 'LineWidth', 4);
%         plot3(0 : dx : (n - 1) * dx, (0 : dx : (n - 1) * dx) * 0 + t, density_x / sum(density_x) / dx,...
%             'LineWidth', 2.5);
%         ylabel('f(x)', 'Interpreter', 'Tex');
        
    hold on;
    
    density_y = squeeze(sum(sum(density_all, 1), 2) * dphi * dx);
%         [~, index] = max(density_y);
%         density_y = circshift(density_x, n / 2 - index);
%         figure; hold on;
        h_y = plot(0 : dy : (m - 1) * dy, density_y / sum(density_y) / dy, 'LineWidth', 4);
%         plot3(0 : dy : (m - 1) * dy, (0 : dy : (n - 1) * dy) * 0 + t, density_y / sum(density_y) / dy,...
%             'LineWidth', 2.5);
%         ylabel('f(y)', 'Interpreter', 'Tex');
        
    box on; grid on;
    legend([h_x, h_y], {'f(x)', 'f(y)'}, 'Interpreter', 'Tex');
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);

end

function [density_xy] = PlotTwoDimSpatialDensity(density_all, n, m, dx, dy, dphi)

    density_xy = squeeze(sum(density_all, 1) * dphi).';
%     [~,max_idx] = max(density_xy(:));
%     [max_row, max_col] = ind2sub(size(density_xy), max_idx);
%     density_xy = CenterTwoDimDensity(density_xy, max_row, max_col);
    
    [mesh_x, mesh_y] = meshgrid(0 : dx : (n - 0) * dx, 0 : dy : (m - 0) * dy);
    density_xy = PadWithPeriodicBoundaries(density_xy, n, m);
    
    figure;
%     surf(mesh_x, mesh_y, density_xy);% * 0 + mean(mean(density_xy)));
    h_pc = pcolor(mesh_x, mesh_y, density_xy);
    set(h_pc, 'EdgeColor', 'None');
%     shading interp;
    colorbar;
    colormap('jet');
    colormap(gca, TransformedColormap);
    
    axis square tight;
%         xlabel('x', 'Interpreter', 'Tex');
%         ylabel('y', 'Interpreter', 'Tex');
%         set(gca, 'YTick', [], 'XTick', []);
    set(gca,...
        'Units', 'normalized',...
        'FontUnits', 'points',...
        'FontWeight', 'normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'linew', 1);
%         'XTickLabel', {'0', '1'},...
%         'YTickLabel', {'0', '1'},...
%         'XTick', [1, n],...
%         'YTick', [1, m]);
    set(gcf, 'Color', 'W');

end

function [density_xy] = CenterTwoDimDensity(density_xy, row, col)

    density_xy = circshift(density_xy, -row + size(density_xy, 1) / 2, 1);
    density_xy = circshift(density_xy, -col + size(density_xy, 2) / 2, 2);

end

function [new_density_xy] = PadWithPeriodicBoundaries(density_xy, n, m)

	new_density_xy = zeros(m + 1, n + 1);
    new_density_xy(1 : m, 1 : n) = density_xy;
    new_density_xy(m + 1, 1 : n) = density_xy(1, 1 : n);
    new_density_xy(1 : m, n + 1) = density_xy(1 : m, 1);
    new_density_xy(m + 1, n + 1) = density_xy(1, 1);

end

function [density_along_parametrized_curve] = PlotSpatialDensityProfileInTheDirectionOfMotion(density_all, n, m, l, dx, dy, dphi)

    density_xy = squeeze(sum(density_all, 1) * dphi).';
    [~,max_idx] = max(density_xy(:));
    [max_row, max_col] = ind2sub(size(density_xy), max_idx);
    density_xy = CenterTwoDimDensity(density_xy, max_row, max_col);
    
    % approximate spatial derivatives with centered differences
    projected_derivative_wrt_x = (circshift(density_all, -1, 2) - circshift(density_all, 1, 2)) / (2 * dx);
    projected_derivative_wrt_x = squeeze(sum(projected_derivative_wrt_x, 1) * dphi).';
    projected_derivative_wrt_x = CenterTwoDimDensity(projected_derivative_wrt_x, max_row, max_col);
    projected_derivative_wrt_y = (circshift(density_all, -1, 3) - circshift(density_all, 1, 3)) / (2 * dy);
    projected_derivative_wrt_y = squeeze(sum(projected_derivative_wrt_y, 1) * dphi).';
    projected_derivative_wrt_y = CenterTwoDimDensity(projected_derivative_wrt_y, max_row, max_col);
    
    % planes with constant \varphi in UxUxT
    orientation_planes = repmat((0 : dphi : (l - 1) * dphi).', [1, n, m]);
    % \int_T cos(\varphi)f(r,\varphi,t) d\varphi
    momentum_x = squeeze(sum(cos(orientation_planes) .* density_all, 1) * dphi).';
    momentum_x = CenterTwoDimDensity(momentum_x, max_row, max_col);
    % \int_T sin(\varphi)f(r,\varphi,t) d\varphi
    momentum_y = squeeze(sum(sin(orientation_planes) .* density_all, 1) * dphi).';
    momentum_y = CenterTwoDimDensity(momentum_y, max_row, max_col);
    momentum_direction = atan2(momentum_y, momentum_x);
    
    [~,max_idx] = max(density_xy(:));
    [max_row, max_col] = ind2sub(size(density_xy), max_idx);
    x_max = (max_col - 1) * dx;
    y_max = (max_row - 1) * dy;
    phi_max = momentum_direction(max_row, max_col)
%     phi_max = mod(phi_max - pi / 2, 2 * pi);
    
    % the curve is parametrized as the one that goes though (x_max,y_max),
    % i.e., the center of the localized cluster,
    % and is parallel to (cos(phi_max),sin(phi_max)), i.e., the direction
    % of the localized cluster
    parametrization = -0.5 : 0.01 : 0.5;
    x = @(s) mod(x_max + cos(phi_max) * s, 1);
    y = @(s) mod(y_max + sin(phi_max) * s, 1);
    i_cell = fix(x(parametrization) / 1 * n);
    j_cell = fix(y(parametrization) / 1 * m);
    % (x_i,y_j) are cell centers of cells C_ij, where (x,y) lie
    x_i = i_cell * dx;
    y_j = j_cell * dy;
    linear_indices = sub2ind(size(density_xy), j_cell+1, i_cell+1);
    density_along_parametrized_curve = density_xy(linear_indices) + projected_derivative_wrt_x(linear_indices) .* (x(parametrization)-x_i) + projected_derivative_wrt_y(linear_indices) .* (y(parametrization)-y_j);
    
    figure;
    hold on;
    plot(parametrization, density_along_parametrized_curve, '.-');
    box on;
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 30,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    
%     figure(); hold on;
%     plot(x(parametrization), y(parametrization), '--w', 'LineWidth', 2);

end

function [] = PlotTwoDimMomentumField(density_all, n, m, l, dx, dy, dphi)

    % marginal density of spatial coordinates in UxU
    density_xy = squeeze(sum(density_all, 1) * dphi).';
    [~,max_idx] = max(density_xy(:));
    [max_row, max_col] = ind2sub(size(density_xy), max_idx);
    % planes with constant \varphi in UxUxT
    orientation_planes = repmat((0 : dphi : (l - 1) * dphi).', [1, n, m]);
    % \int_T cos(\varphi)f(r,\varphi,t) d\varphi
    momentum_x = squeeze(sum(cos(orientation_planes) .* density_all, 1) * dphi).';
    momentum_x = CenterTwoDimDensity(momentum_x, max_row, max_col);
    % \int_T sin(\varphi)f(r,\varphi,t) d\varphi
    momentum_y = squeeze(sum(sin(orientation_planes) .* density_all, 1) * dphi).';
    momentum_y = CenterTwoDimDensity(momentum_y, max_row, max_col);
    
    momentum_x = PadWithPeriodicBoundaries(momentum_x, n, m);
    momentum_y = PadWithPeriodicBoundaries(momentum_y, n, m);    
    [mesh_x, mesh_y] = meshgrid(0 : dx : (n - 0) * dx, 0 : dy : (m - 0) * dy);
    submesh_x = 1 : n/10 : n + 1;
    submesh_y = 1 : m/10 : m + 1;    
    
    auto_scale_factor = 1.5;
    figure;
    h_q = quiver(mesh_x(submesh_y, submesh_x), mesh_y(submesh_y, submesh_x), momentum_x(submesh_y, submesh_x), momentum_y(submesh_y, submesh_x), 'AutoScaleFactor', auto_scale_factor);
    figure(); hold on;
    U = h_q.UData;
    V = h_q.VData;
    X = h_q.XData;
    Y = h_q.YData;
    for j = 1 : size(X, 1)
        for i = 1 : size(X, 2)
            line_length = norm([U(j,i) V(j,i)]) * auto_scale_factor * 0.03; % 0.02
            head_length = norm([U(j,i) V(j,i)]) * auto_scale_factor * 8; % 3
            head_width = norm([U(j,i) V(j,i)]) * auto_scale_factor * 6; % 2
            ah = annotation('arrow',...
                'headStyle', 'cback1', 'HeadLength', head_length, 'HeadWidth', head_width, 'Color', 'White');
            set(ah,'parent', gca);
            set(ah,'position',[X(j,i) Y(j,i) line_length*U(j,i) line_length*V(j,i)]);
        end % i
    end % j
    box on; grid on;
    axis square tight;
%     xlabel('x'); ylabel('y');
%         set(gca, 'YTick', [], 'XTick', []);
%     set(gca,...
%         'Units', 'normalized',...
%         'FontUnits', 'points',...
%         'FontWeight', 'normal',...
%         'FontSize', 30,...
%         'FontName', 'Helvetica',...
%         'linew', 1,...
%         'XTickLabel', {'0', '1'},...
%         'YTickLabel', {'0', '1'},...
%         'XTick', [1, n],...
%         'YTick', [1, m]);
%         set(gcf, 'Color', 'None');

end