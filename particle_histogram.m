clc;
    
N = 10000;
S = 3;
starting_time = 0;
dt_recip = 1;

folder = '/Users/nikita/Documents/Projects/spc2/spc2OdeIntegration/';
% folder = '/Volumes/Kruk/spc2/spc2OdeIntegration/continued/';
file_name = 'v0_1_sigma_1_rho_0.01_alpha_0_Dphi_0.01_N_10000_0_0.bin';
file_id = fopen(strcat(folder, file_name));
size_of_float = 4;
status = fseek(file_id,...
    (1 + S * N) * starting_time * dt_recip * size_of_float, 'bof'); % skip t_0 steps
assert(~status);
% fread(file_id, [1 + S * N, 1 * dt_recip], 'float'); % skip t_1 seconds
SystemState = fread(file_id, [1 + S * N, 1], 'float');
SystemState(1, 1)

figure%(5); hold on;
wrapped_phase_distribution = mod(SystemState(4 : S : end, 1), 2 * pi);
phase_mean = mean(wrapped_phase_distribution, 1);
wrapped_phase_distribution = mod(wrapped_phase_distribution - phase_mean + pi, 2 * pi);
histogram(wrapped_phase_distribution, 10, 'Normalization', 'Pdf');
xlim([0, 2 * pi]);
% title('Phase distribution');
% angle(mean(exp(1i*wrapped_phase_distribution)))
set(gca,...
    'Units', 'normalized',...
    'FontUnits', 'points',...
    'FontWeight', 'normal',...
    'FontSize', 20,...
    'FontName', 'Helvetica',...
    'linew', 1);

figure;
histogram2(SystemState(2 : S : end, 1), SystemState(3 : S : end, 1),...
    'FaceColor', 'Flat', 'Normalization', 'Pdf');
title('Position distribution');

% figure;
% % x = SystemState(2 : S : end, 1);
% % y = SystemState(3 : S : end, 1);
% vx = cos(wrapped_phase_distribution);
% vy = sin(wrapped_phase_distribution);
% histogram2(vx, vy, 25,...
%     'FaceColor', 'Flat', 'Normalization', 'Pdf');
% title('Momentum distribution');