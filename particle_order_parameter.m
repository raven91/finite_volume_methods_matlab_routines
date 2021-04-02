folder = '/Users/nikita/Documents/Projects/spc2/spc2OdeIntegration/';
% folder = '/Volumes/Kruk/spc2/spc2OdeIntegration/continued/';
file_name = strcat('summary_statistics_v0_1_sigma_1_rho_0.01_alpha_0_Dphi_0.01_N_10000_0_0.txt');
order_parameters = dlmread(strcat(folder, file_name));
figure;
plot(order_parameters(:, 1), order_parameters(:, 2), '-', 'LineWidth', 2.5);
hold on; 
grid on; box on;
ylim([0, 1]);