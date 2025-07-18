clc; clear; close all;

alpha = 3; % Fixed alpha
sigma_range = linspace(0.1, 5, 100); % Range of sigma
beta_range = linspace(0.1, 5, 100);  % Range of beta
t_span = [0, 100]; % Time span for simulation

% --- ODE Solver Options ---
RelTol = 1e-6;  % Relative tolerance
AbsTol = 1e-8;  % Absolute tolerance
ode_opts = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

% --- Parameters ---
N = 10;
A = zeros(N);

% Each node connects to i+1 and i+3 (mod N)
for i = 1:N
    A(mod(i, N)+1, i) = 1;         % i → i+1
    A(mod(i+2, N)+1, i) = 1;       % i → i+3
end

% Verify balance and strong connectivity
assert(all(sum(A,1)' == sum(A,2)), 'Graph is not balanced');
G = digraph(A);

% Step 2: Add long-range connections to increase connectivity
for i = 1:N
    A(mod(i+2-1, N)+1, i) = 1;     % i → i+2
    A(mod(i-3, N)+1, i) = 1;       % i → i−3
end

% Step 3: Check balance
assert(all(sum(A,1)' == sum(A,2)), 'Graph is not balanced');

% Step 4: Check strong connectivity
G = digraph(A);
% assert(isconnected(G, 'strong'), 'Graph is not strongly connected');

% Optional: visualize the graph
figure;
plot(G, 'Layout', 'circle', 'NodeLabel', 1:N);
title('Strongly Connected & Balanced Digraph with Higher Connectivity');

% --- Algebraic connectivity ---
% Compute algebraic connectivity
out_deg = sum(A, 2);
L = diag(out_deg) - A;

a = algebraic_connectivity(L);

% --- Step 7: Initialize error storage ---
error_map = zeros(length(beta_range), length(sigma_range));

% --- Step 8: Loop over beta and sigma ---
for b_idx = 1:length(beta_range)
    beta = beta_range(b_idx);
    for s_idx = 1:length(sigma_range)
        sigma = sigma_range(s_idx);
        
        % Initial conditions for x and v (randomized)
        X0 = rand(2*N, 1) - 0.5;
        
        % Define the system dynamics
        dynamics = @(t, X) system_dynamics(t, X, N, alpha, beta, sigma, L);
        
        % Solve the system using ODE45 with specified tolerances
        [t, X] = ode45(dynamics, t_span, X0);%, ode_opts);
        
        % Extract final states of nodes
        x_nodes = X(end, 1:2:end); % Positions (x_i) at final time
        x_avg = mean(x_nodes);    % Average of the positions
        error = mean(abs(x_nodes - x_avg)); % Error metric
        
        % Store error in the map
        error_map(b_idx, s_idx) = error;
    end
end

% --- Step 9: Plot heat map of consensus error ---
figure;
imagesc(sigma_range, beta_range, error_map);
set(gca, 'YDir', 'normal');

ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLabelInterpreter = 'latex';
set(gca, 'FontSize', 20);

c = colorbar;
c = colorbar;

% Colorbar label
ylabel(c, '$||e||_2$', 'Interpreter','latex', 'FontSize', 25, 'Color', 'k');

% Colorbar tick marks and text
c.Color = 'k';                     % Tick color
c.TickLabelInterpreter = 'latex'; % Consistent formatting
c.FontSize = 20;

xlabel('$\sigma$', 'Interpreter','latex', 'FontSize', 25, 'Color', 'k');
ylabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 25, 'Color', 'k');

hold on;
plot(sigma_range, 1/4 * sigma_range * a, 'r', 'LineWidth', 2);
hold off;

set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardcopy', 'off');

exportgraphics(gcf, 'Figure_2_d.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');

% --- Step 10: Plot x_i(t) trajectories ---
figure;
hold on;
x = X(:, 1:2:end);
for i = 1:size(x,2)
    plot(t, x(:, i), 'LineWidth', 1.5)
end
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 25); 
ylabel('$x_i(t)$', 'Interpreter', 'latex', 'FontSize', 25); 
grid off;
set(gca, 'LooseInset', get(gca, 'TightInset'));
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca, 'FontSize', 20);
exportgraphics(gcf, ['x_trajectory_full_directed.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');

% --- System dynamics function ---
function dX = system_dynamics(~, X, N, alpha, beta, sigma, L)
    x = X(1:2:end); % Extract x_i
    v = X(2:2:end); % Extract v_i
    dx = v; % dx_i = v_i
    dv = -alpha * v - 4 * beta * (x.^3 - x) - sigma * (L * x); % dv_i
    dX = zeros(2*N, 1);
    dX(1:2:end) = dx;
    dX(2:2:end) = dv;
end


function a = algebraic_connectivity(L)
    % Check size
    N = size(L, 1);

    % Compute left eigenvector beta of L corresponding to eigenvalue 0
    [V, D] = eig(L');
    zero_eig_idx = find(abs(diag(D)) < 1e-10, 1);
    beta = real(V(:, zero_eig_idx));
    beta = beta / sum(beta);  % normalize such that beta' * 1 = 1

    % Construct weighted Laplacian: L_hat = 0.5 * (diag(beta)*L + L'*diag(beta))
    B = diag(beta);
    L_hat = 0.5 * (B * L + L' * B);

    % Define constraint x' * beta = 0 using null space
    Q = null(beta');  % N x (N-1) orthonormal basis s.t. Q'*beta = 0

    % Project generalized eigenproblem to subspace orthogonal to beta
    A = Q' * L_hat * Q;
    M = Q' * B * Q;

    % Solve the generalized eigenvalue problem: A z = lambda * M z
    eigvals = eig(A, M);
    eigvals = real(eigvals);  % Remove small imaginary parts due to round-off

    % Filter out nonpositive or near-zero eigenvalues (numerical threshold)
    tol = 1e-10;
    eigvals = eigvals(eigvals > tol);

    if isempty(eigvals)
        error('All computed eigenvalues are non-positive or numerically zero.');
    end

    % Algebraic connectivity is the smallest strictly positive eigenvalue
    a = min(eigvals);
end


