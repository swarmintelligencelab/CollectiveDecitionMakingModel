close all
clear all
clc

% Parameters
alpha = 0.5;
sigma_range = linspace(0.1, 5, 391);
n_range = linspace(10, 300, 291);

beta = 0.01;
t_span = [0, 50];
num_trials = 2000;

% --- ODE Solver Options ---
RelTol = 1e-6;  % Relative tolerance
AbsTol = 1e-8;  % Absolute tolerance
ode_opts = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

% Preallocation
error_map = zeros(length(n_range), length(sigma_range));
trial_error_map = zeros(length(n_range), length(sigma_range), num_trials);
lambda_2_list = zeros(size(n_range));

% Loop over N and sigma
for n_idx = 1:length(n_range)
    N = round(n_range(n_idx));
    fprintf('N = %d/%d\n', N, round(n_range(end)));

    % Generate ring network and Laplacian
    A = generate_ring_network(N, min(4, N-1));
    L = diag(sum(A,2)) - A;

    % Store lambda_2
    lambda_vals = sort(eig(L));
    lambda_2_list(n_idx) = lambda_vals(2);

    for s_idx = 1:length(sigma_range)
        sigma = sigma_range(s_idx);
        trial_errors = zeros(1, num_trials);

        for trial = 1:num_trials
            X0 = rand(2*N, 1) - 0.5;  % Random initial condition
            dynamics = @(t, X) system_dynamics(t, X, N, alpha, beta, sigma, L);
            [~, X] = ode45(dynamics, t_span, X0, ode_opts);
            x_nodes = X(end, 1:2:end);
            x_avg = mean(x_nodes);
            trial_errors(trial) = mean(abs(x_nodes - x_avg));
        end

        % Store mean and all trial errors
        error_map(n_idx, s_idx) = mean(trial_errors);
        trial_error_map(n_idx, s_idx, :) = trial_errors;
    end
end

% --- Save results ---
save('fish_synchronization_error_map.mat', ...
    'error_map', 'trial_error_map', 'sigma_range', 'n_range', ...
    'alpha', 'beta', 't_span', 'lambda_2_list');

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
