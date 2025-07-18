clc; clear; close all;

%% Parameters
x_op = 1;
nu = 1;
sigma = 10;
N = 100;
p = 4*log(N)/N;

A_full = ones(N) - eye(N);
A = A_full;

L = diag(sum(A, 2)) - A;
lambda_max = max(eig(L));

beta_upper = (1 / (4 * x_op^2)) * ( sqrt( sigma * lambda_max * (sigma * lambda_max - 2*nu) ) - 1 );
beta = 1;
alpha_upper = 8/(3*sqrt(3))*beta*x_op^3;

alpha_values = [0, 0.5*alpha_upper, 0.8*alpha_upper, -0.8*alpha_upper];
figure_labels = {'3a', '3b', '3c', '3d'};

for idx = 1:length(alpha_values)
    alpha = alpha_values(idx);
    theta = acos(- (3*sqrt(3)*alpha) / (8*beta*x_op^3));
    x_star = 2*x_op/sqrt(3) * cos((theta + 2*pi*(0:2))/3);
    min_eq = min(x_star);
    max_eq = max(x_star);
    x_start_sorted = sort(abs(x_star));
    x_c = x_op/sqrt(3);
    r = abs(min(abs(x_start_sorted(2)) - x_c));
    
    [x_vals, v_vals] = meshgrid(linspace(-3,3,150), linspace(-3,3,150));
    basin_idx = zeros(size(x_vals));
    dx = zeros(size(x_vals));
    dv = zeros(size(v_vals));
    
    for i = 1:numel(x_vals)
        xi = x_vals(i);
        vi = v_vals(i);
        dx(i) = vi;
        dv(i) = -nu*vi - 4*beta*xi^3 + 4*beta*x_op^2*xi - alpha; 
        [~, X] = ode45(@(t,X) [X(2); -nu*X(2) - 4*beta*X(1)^3 + 4*beta*x_op^2*X(1) - alpha], [0 30], [xi; vi]);
        x_final = X(end,1);
        v_final = X(end,2);
        [~, min_idx] = min(abs(x_final - x_star));

        if abs(v_final) < 1e-3
            if norm(x_final - min_eq) < 1e-3
                basin_idx(i) = min_eq;
            elseif norm(x_final - max_eq) < 1e-3
                basin_idx(i) = max_eq;
            else
                basin_idx(i) = 0;
            end
        end
    end
    
    % Create a new figure for each alpha
    figure;
    contourf(x_vals, v_vals, basin_idx, 50, 'LineColor', 'none');
    colormap([1 0 0; 0 0 1]);
    
    % Only add colorbar on the final plot
    if idx == length(alpha_values)
        cb = colorbar;
        cb.Ticks = [];
    end
    
    xlabel('$\overline{x}$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
    ylabel('$\overline{v}$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
    
    % Title per alpha
    switch idx
        case 1
            title('$\alpha = 0$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
        case 2
            title('$\alpha = 0.5 \alpha_\mathrm{c}$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
        case 3
            title('$\alpha = 0.8 \alpha_\mathrm{c}$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
        case 4
            title('$\alpha = -0.8 \alpha_\mathrm{c}$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');
    end

    % Draw stability circles
    hold on;
    for k = 1:2
        x_center = [min(x_star), max(x_star)];
        r = abs(min(abs(x_center(k)) - x_c));
        plot(x_center(k) + r*cos(linspace(0, 2*pi, 100)), r*sin(linspace(0, 2*pi, 100)), ...
             'g-', 'LineWidth', 1.5);
    end
    hold off;

    set(gca, 'FontSize', 25, 'XColor', 'k', 'YColor', 'k');
    set(gcf, 'Color', 'none');
    set(gca, 'Color', 'none');

    exportgraphics(gcf, ['Figure_' figure_labels{idx} '.pdf'], ...
                   'BackgroundColor', 'none', 'ContentType', 'vector');
end

