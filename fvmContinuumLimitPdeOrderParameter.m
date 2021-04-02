function [] = fvmContinuumLimitPdeOrderParameter
global n m l dt_recip t0 t1 dphi c_phi diff_vol

    clc;
    
    n = 40;
    m = 40;
    l = 256;
    dt_recip = 1;
    t0 = 0 * dt_recip;
    t1 = 1 * dt_recip;
    x_size = 1;
    y_size = 1;
    phi_size = 2 * pi;
    dx = x_size / n;
    dy = y_size / m;
    dphi = phi_size / l;
    c_phi = 2 * pi / phi_size;
    diff_vol = dx * dy * dphi;
    
%     folder = '/Users/nikita/Documents/Projects/spc2/spc2FiniteVolumeMethods/';
    folder = '/Volumes/Kruk/spc2/spc2FiniteVolumeMethods/Continuation/';
    file_name = strcat('dt_0.005_sigma_1_rho_0.3_alpha_1.45_Dphi_0.0075_', num2str(n), '_', num2str(m), '_', num2str(l), '_1000.txt');
    
%     hold on;
    figure;
%     [~, ~] = FromBinary(folder, file_name);
    [~] = FromText(folder, file_name);
    
    ylim([0 1]);
    box on;
    grid on;
    xlabel('t', 'Interpreter', 'Tex');
    ylabel('R', 'Interpreter', 'Tex');
    hold on;
    
    set(gca,...
        'Units', 'Normalized',...
        'FontUnits', 'Points',...
        'FontWeight', 'Normal',...
        'FontSize', 24,...
        'FontName', 'Helvetica',...
        'LineW', 1);
    
end

function [OrderParameter, OrderParameterNormalization] = FromBinary(folder, file_name)
global n m l dt_recip t0 t1 dphi c_phi diff_vol

    file_id = fopen(strcat(folder, file_name));
    fread(file_id, [1 + n * m * l, t0], 'float'); % skip t0 steps

    OrderParameter = zeros(t1 / dt_recip, 1);
    OrderParameterNormalization = zeros(t1 / dt_recip, 1);
    for t = t0 : dt_recip : t1 - 1
        ContinuumLimitPdf = fread(file_id, [1 + n * m * l, 1], 'float');
%         ContinuumLimitPdf = transpose(ContinuumLimitPdf);
        fprintf('t=%f\n', ContinuumLimitPdf(1, 1))
%         Density = reshape(ContinuumLimitPdf(2 : n * m * l + 1, 1), n, m, l);
        Density = reshape(ContinuumLimitPdf(2 : n * m * l + 1, 1), l, n, m);
        
        OrderParameterNormalization(t / dt_recip + 1, 1) = sum(sum(sum(Density))) * diff_vol;
        for p = 1 : 1 : l
            OrderParameter(t / dt_recip + 1, 1) = OrderParameter(t / dt_recip + 1, 1) + sum(sum((cos(c_phi * (p - 1) * dphi) + 1i * sin(c_phi * (p - 1) * dphi)) .* Density(p, :, :)));
        end
        OrderParameter(t / dt_recip + 1, 1) = OrderParameter(t / dt_recip + 1, 1) * diff_vol;
        
        fread(file_id, [1 + n * m * l, dt_recip - 1], 'float');
    end % for t
    fclose(file_id);
    
    plot(t0 / dt_recip : t1 / dt_recip - 1,...
        abs(OrderParameter(t0 / dt_recip + 1 : t1 / dt_recip, 1))...
            ./ OrderParameterNormalization(t0 / dt_recip + 1 : t1 / dt_recip, 1),...
        '-', 'LineWidth', 2.5);%, 'Color', [1 0.75 0] );
    
end

function [OrderParameter] = FromText(folder, file_name)

	OrderParameters = dlmread(strcat(folder, file_name));
    OrderParameter = OrderParameters(:, 2);
%     min(OrderParameters(2:end,1)-OrderParameters(1:end-1,1))
%     max(OrderParameters(2:end,1)-OrderParameters(1:end-1,1))
    
    plot(OrderParameters(:, 1), OrderParameters(:, 2), '-', 'LineWidth', 2.5);
%     hold on;
%     plot(OrderParameters(2:end-1,1), (OrderParameters(3:end, 2) - 2 * OrderParameters(2:end-1, 2) + OrderParameters(1:end-2, 2))/(OrderParameters(2,1)-OrderParameters(1,1)), '-', 'LineWidth', 2.5);

%     t_max = size(OrderParameters, 1);
%     mu = sum(OrderParameters(t_max - 100 : t_max, 2)) / length(OrderParameters(t_max - 100 : t_max, 2));
%     variance = sum((OrderParameters(t_max - 100 : t_max, 2) - mu) .^ 2) / (length(OrderParameters(t_max - 100 : t_max, 2)) - 1);
%     sprintf('%g, %g', mu, sqrt(variance))

end