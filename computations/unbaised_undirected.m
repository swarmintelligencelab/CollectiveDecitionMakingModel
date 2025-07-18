close all
clear all
clc

% Parameters
alpha = 2;
sigma_range = linspace(0.1, 5, 100);
beta_range = linspace(0.1, 5,100);
N = 10;
t_span = [0, 50];

% --- ODE Solver Options ---
RelTol = 1e-6;  % Relative tolerance
AbsTol = 1e-8;  % Absolute tolerance
ode_opts = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

% Graph
A_full = ones(N) - eye(N);
A = A_full;

L = diag(sum(A, 2)) - A;
l = eig(L);
lambda_max = max(l(2));

error_map = zeros(length(beta_range), length(sigma_range));

% Loop over beta and sigma
for b_idx = 1:length(beta_range)
    beta = beta_range(b_idx);
    for s_idx = 1:length(sigma_range)
        sigma = sigma_range(s_idx);
        
        X0 = rand(2*N, 1) - 0.5;
        dynamics = @(t, X) system_dynamics(t, X, N, alpha, beta, sigma, L);
        [t, X] = ode45(dynamics, t_span, X0);%, ode_opts);

        x_nodes = X(end, 1:2:end);
        x_avg = mean(x_nodes);
        error = mean(abs(x_nodes - x_avg));
        
        error_map(b_idx, s_idx) = error;
    end
end

% --- Heatmap Plot ---
figure;
imagesc(sigma_range, beta_range, error_map);
set(gca, 'YDir', 'normal');

ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLabelInterpreter = 'latex';
set(gca, 'FontSize', 20);

% c = colorbar;
% c = colorbar;
% 
% % Colorbar label
% ylabel(c, '$||e||_2$', 'Interpreter','latex', 'FontSize', 25, 'Color', 'k');
% 
% % Colorbar tick marks and text
% c.Color = 'k';                     % Tick color
% c.TickLabelInterpreter = 'latex'; % Consistent formatting
% c.FontSize = 20;

xlabel('$\sigma$', 'Interpreter','latex', 'FontSize', 25, 'Color', 'k');
ylabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
xticks([0,1,2,3,4,5])

hold on;
plot(sigma_range, 1/4 * sigma_range * lambda_max, 'r', 'LineWidth', 2);
hold off;

set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardcopy', 'off');

exportgraphics(gcf, 'Figure_2_c.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');

% --- Trajectory Plot ---
figure;
hold on

x = X(:, 1:2:end);
if any(~isfinite(x(:)))
    warning('NaN or Inf detected in trajectory data.');
end
x(~isfinite(x)) = NaN; % Clean non-finite values

for i = 1:size(x,2)
    plot(t, x(:, i), 'LineWidth', 1.5);
end

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
ylabel('$x_i(t)$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');

ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLabelInterpreter = 'latex';
set(gca, 'FontSize', 20);
set(gca, 'LooseInset', get(gca, 'TightInset'));
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardcopy', 'off');

exportgraphics(gcf, 'x_trajectory_full.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');


% --- System Dynamics Function ---
function dX = system_dynamics(~, X, N, alpha, beta, sigma, L)
    x = X(1:2:end);
    v = X(2:2:end);
    dx = v;
    dv = -alpha * v - 4 * beta * (x.^3 - x) - sigma * (L * x);
    dX = zeros(2*N, 1);
    dX(1:2:end) = dx;
    dX(2:2:end) = dv;
end


