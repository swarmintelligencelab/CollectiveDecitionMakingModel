clc; clear; close all;

% --- Load data ---
data = load('fish_synchronization_error_map.mat');

error_map = data.error_map(1:30,:);
n_range = data.n_range(1:30);
sigma_range = data.sigma_range;

[N_grid, S_grid] = meshgrid(sigma_range, n_range);

% --- Parameters ---
beta = data.beta;

% --- Compute sigma(n) for overlay ---
n_fine = linspace(min(n_range), max(n_range), 1000);
sigma_overlay = 4 * beta ./ (4 - 2*cos(2*pi./n_fine) - 2*cos(4*pi./n_fine));

%% --- 2D Heatmap View with Overlay ---
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6, 5]);

imagesc(sigma_range, n_range, error_map);
axis xy;
colormap(turbo);
colorbar_handle = colorbar;

hold on;
plot(sigma_overlay, n_fine, 'r-', 'LineWidth', 2);  % (sigma(n), n)

% --- Axis and colorbar styling ---
xlabel('$\sigma$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
ylabel('$n$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');

ylabel(colorbar_handle, '$\|e\|_2$', ...
    'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');

set(gca, 'FontSize', 25, ...
         'XColor', 'k', 'YColor', 'k', ...
         'TickLabelInterpreter', 'latex');
set(colorbar_handle, 'FontSize', 25, 'Color', 'k', ...
                     'TickLabelInterpreter', 'latex');

% --- Export to PDF ---
exportgraphics(fig1, 'sync_error_map_2D_with_sigma_curve.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
