% Nu% Number of nodes
N = 10;

% Define edges (i → i+1 and i → i+3 mod N)
s = [];
t = [];
for i = 1:N
    s = [s, i, i]; % From node i
    t = [t, mod(i, N)+1, mod(i+2, N)+1]; % To i+1 and i+3
end

% Create directed graph
G = digraph(s, t);

% Custom node positions (you can keep this layout or use circle layout)
x = [0, 1, 1, 0, -1.5, 2.5, 2.5, -1.5, -2.5, 3.5];
y = [1, 1, 0, 0, 1, 1, 0, 0, 1, 1];

% Plot
figure;
h = plot(G, 'XData', x, 'YData', y, ...
         'NodeColor', 'r', 'NodeLabel', {}, ...
         'EdgeColor', 'k', 'ArrowSize', 15, ...
         'MarkerSize', 10, 'LineWidth', 1.2);

% Add white labels for nodes
hold on;
for i = 1:numnodes(G)
    text(x(i), y(i), num2str(i), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 10);
end

% Improve layout
padding = 0.5;
xlim([min(x)-padding, max(x)+padding]);
ylim([min(y)-padding, max(y)+padding]);

axis equal;
axis off;

% Export figure
exportgraphics(gcf, 'Figure_2_b.pdf', ...
    'BackgroundColor', 'none', 'ContentType', 'vector');
