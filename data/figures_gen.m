% Load data for TRPV1
shell_values = 1:5;
load('rel_errors_TRPV1_shells_ell5.mat')
rel_errors_TRPV1 = rel_errors;

% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 4.5, 3.5]); % 1-column format

% Plot TRPV1
errorbar(shell_values, mean(rel_errors_TRPV1, 2), std(rel_errors_TRPV1, 0, 2), ...
    'o-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.2 0.7]);

hold on;

% Load data for 80S
load('rel_errors_80S_shells_ell5.mat')
rel_errors_80S = rel_errors;

% Plot 80S with different color/marker
errorbar(shell_values, mean(rel_errors_80S, 2), std(rel_errors_80S, 0, 2), ...
    's-', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);

% Labels and formatting
xlabel('Number of Shells', 'FontSize', 11);
ylabel('Average Relative Error', 'FontSize', 11);
title('Recovery Error vs. Number of Shells', 'FontSize', 11);
legend({'TRPV1', '80S'}, 'Location', 'northeast', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1);
box on;

% Set paper size for tight export
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperSize', [4.5 3.5]);

% Export as vector PDF (tight bounding box)
print(fig, 'freq_marching_error_vs_shells.pdf', '-dpdf', '-painters');
