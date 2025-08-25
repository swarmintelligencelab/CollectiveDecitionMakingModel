close all
clear all
clc

% Parameters
alpha = -1.06;
sigma = 1;   % fixed
n_range = linspace(10, 100, 80);
nu = 1;
beta = 1.0744;
t_span = [0, 50];
num_trials = 100;

% Fractions of agents initialized at (-1, 0)
frac_range = linspace(0, 1, 80);

% ODE Solver Options
RelTol = 1e-6;
AbsTol = 1e-8;
ode_opts = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

% Preallocation
conv_frac_map = zeros(length(n_range), length(frac_range));

% Tolerance for convergence check
x_tol = 0.05;   % tolerance for being near -1
v_tol = 0.05;   % tolerance for velocities being near 0

% Enable parallel pool
if isempty(gcp('nocreate'))
     parpool(6);
end

disp('Starting simulation...');
for n_idx = 1:length(n_range)
    N = round(n_range(n_idx));
    fprintf('N = %d/%d\n', N, round(n_range(end)));

    % Generate ring network and Laplacian
    A = generate_ring_network(N, min(4, N-1));
    L = diag(sum(A,2)) - A;

    for f_idx = 1:length(frac_range)
        frac = frac_range(f_idx);
        num_fixed = round(frac * N);

        conv_count = 0; % count number of successful convergences

        % --- Parallel loop over trials ---
        parfor trial = 1:num_trials
            % Initialize X0
            X0 = rand(2*N, 1)*2 - 1;  % random init

            % Force a fraction to (-1, 0)
            chosen = randperm(N, num_fixed);
            for k = chosen
                X0(2*k-1) = -0.8415;   % x = -1
                X0(2*k)   = 0;    % v = 0
            end

            dynamics = @(t, X) system_dynamics(t, X, N, nu, beta, sigma, L, alpha);
            [~, X] = ode45(dynamics, t_span, X0, ode_opts);

            % Final states
            x_nodes = X(end, 1:2:end);
            v_nodes = X(end, 2:2:end);

            % Check convergence: all agents near (-1,0)
            if all(abs(x_nodes - -0.8415) < x_tol) && all(abs(v_nodes) < v_tol)
                conv_success = 1;
            else
                conv_success = 0;
            end

            % Accumulate
            conv_count = conv_count + conv_success;
        end

        % Fraction of successful trials
        conv_frac_map(n_idx, f_idx) = conv_count / num_trials;
    end
end

save('fish_sync_convergence_fraction.mat', ...
    'conv_frac_map', 'n_range', 'frac_range', ...
    'alpha', 'beta', 't_span', 'sigma');

% --- Plot ---
figure;
set(gcf,'Color','w');   % <-- figure background white
[X, Y] = meshgrid(frac_range, n_range);
surf(X, Y, conv_frac_map, 'EdgeColor', 'none');

% Labels
xlabel('$p$', 'Interpreter','latex','FontSize',25,'Color','k');
ylabel('$n$', 'Interpreter','latex','FontSize',25,'Color','k');
ylim([10 max(n_range)]);

% Set axes appearance
set(gca,'FontSize',25,'XColor','k','YColor','k','ZColor','k','Color','w');

% Colormap and colorbar
colormap jet;
cb = colorbar;
ylabel(cb, '$q$', 'Interpreter','latex','FontSize',25,'Color','k');
set(cb, 'Color', 'k'); % ticks/labels black

% View from top
view(2);


% --- System Dynamics Function ---
function dX = system_dynamics(~, X, N, nu, beta, sigma, L, alpha)
    x = X(1:2:end);
    v = X(2:2:end);
    dx = v;
    dv = -nu * v - 4 * beta * (x.^3 - x) - sigma * (L * x) - alpha;
    dX = zeros(2*N, 1);
    dX(1:2:end) = dx;
    dX(2:2:end) = dv;
end

% --- Ring Network Generator ---
function A = generate_ring_network(n, k)
    A = zeros(n);
    half_k = floor(k / 2);
    for i = 1:n
        for j = 1:half_k
            A(i, mod(i-1 + j, n) + 1) = 1;
            A(i, mod(i-1 - j, n) + 1) = 1;
        end
    end
end
