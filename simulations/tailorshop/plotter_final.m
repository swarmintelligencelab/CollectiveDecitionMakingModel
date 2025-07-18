clear; close all; clc;

% List of file pairs
file_pairs = {
    'pairwise_basins_1_Ironer_Cotton_Boy.mat',       'pairwise_basins_2_Ironer_Cotton_Boy.mat';
    'pairwise_basins_1_Button_Machiner_Cotton_Boy.mat','pairwise_basins_2_Button_Machiner_Cotton_Boy.mat';
    'pairwise_basins_1_Tailor_Cotton_Boy.mat',       'pairwise_basins_2_Tailor_Cotton_Boy.mat';
    'pairwise_basins_1_Tailor_Ironer.mat',           'pairwise_basins_2_Tailor_Ironer.mat';
    'pairwise_basins_1_Tailor_Button_Machiner.mat',  'pairwise_basins_2_Tailor_Button_Machiner.mat';
    'pairwise_basins_1_Button_Machiner_Ironer.mat',  'pairwise_basins_2_Button_Machiner_Ironer.mat';
};

% Output folder
output_dir = 'original_basins';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:6
    for j = 1:2
        % Load data
        data = load(file_pairs{i, j});

        fig = figure('Visible', 'off', 'Color', 'w', ...
             'Units', 'inches', 'Position', [1, 1, 4, 4]);  % 4x4 inches
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.13 0.11 0.775 0.815]);

% Plot on fixed-size axes
contourf(ax, data.x_vals, data.v_vals, data.basin_idx, 2, 'LineColor', 'none');
colormap([1 0 0; 0 0 1]); % Red = +1, Blue = -1
axis(ax, 'square');
xlabel(ax, '$\bar{x}$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel(ax, '$\bar{v}$', 'Interpreter', 'latex', 'FontSize', 24);
set(ax, 'FontSize', 24);
xticks(ax, [-2 2]);
yticks(ax, [-2 2]);
set(ax, 'XColor', 'k', 'YColor', 'k');           % Axis tick color
ax.XRuler.TickLabelInterpreter = 'latex';       % Optional: if using LaTeX
ax.YRuler.TickLabelInterpreter = 'latex';

if strcmp(file_pairs{i, j}, 'pairwise_basins_2_Button_Machiner_Cotton_Boy.mat')
cb = colorbar(ax);
cb.Ticks = [-1 1];
cb.TickLabels = {'', ''};

end

name = sprintf('%s_%s_data_%d', data.group1, data.group2, j);
        name = strrep(name, ' ', '_');
        filepath = fullfile(output_dir, [name, '.pdf']);

% Save entire figure with controlled size
exportgraphics(fig, filepath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
