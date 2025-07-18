clear all
close all


%% Parameters
x_op = 1;
nu = 1;
sigma = 3;
N = 39; % Number of nodes
beta = 1;
t_span = [0, 50];
grid_size = 270;
[x_vals, v_vals] = meshgrid(linspace(-3, 3, grid_size), linspace(-3, 3, grid_size));

groups = cell(1,39);
groups{19} = {'Head Tailor'};
groups{16} = {'Cutter'};
groups([25,26]) = {'Button Machiner'};
groups([29,33,39]) = {'Ironer'};
groups([30:32,34:38]) = {'Cotton Boy'};
groups([1:3,5:7,9,11:14,21,24,4,10,17,18,8,15,20,22:23,27:28]) = {'Tailor'};


% Define the group pairs
group_pairs = {
    {'Tailor', 'Button Machiner'};
    {'Tailor', 'Ironer'};
    {'Tailor', 'Cotton Boy'};
    {'Button Machiner', 'Ironer'};
    {'Button Machiner', 'Cotton Boy'};
    {'Ironer', 'Cotton Boy'};
    };

% --- Start plotting ---
figure;
tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

for pair_id = 1:length(group_pairs)
    group1 = group_pairs{pair_id}{1};
    group2 = group_pairs{pair_id}{2};

    % Find nodes in selected groups
    group1_idx = find(strcmp(groups, group1));
    group2_idx = find(strcmp(groups, group2));
    fixed_group_idx = unique([group1_idx, group2_idx]);

    % Head Tailor and Cutter fixed
    X_fixed = zeros(2*N, 1);
    X_fixed(N*2-1) = -1;
    X_fixed(N*2) = 0;
    % X_fixed(19*2-1) = -1; X_fixed(19*2) = 0; % Head Tailor
    % X_fixed(fixed_group_idx*2 - 1) = 1;     % group1 + group2: x = 1
    % X_fixed(fixed_group_idx*2)     = 0;     % group1 + group2: v = 0

    % Load graph (pick kapfts1.dat)
    A = load('kapfts1.dat');
    D = diag(sum(A, 2));
    L = D - A;

    basin_idx = zeros(grid_size, grid_size);

    % Loop over grid
    for i = 1:grid_size
        for j = 1:grid_size
            X0 = X_fixed;

            for k = 1:N
                if ~any(k == fixed_group_idx)
                    continue;
                end
                X0(2*k - 1) = x_vals(i,j);
                X0(2*k)     = v_vals(i,j);
            end

            dynamics = @(t, X) system_dynamics(t, X, N, nu, beta, sigma, L);
            [~, X_sol] = ode45(dynamics, t_span, X0);

            x_final = mean(X_sol(end, 1:2:end));
            v_final = mean(X_sol(end, 2:2:end));

            if abs(v_final) < 1e-3
                if abs(x_final - (-1)) < 1e-3
                    basin_idx(i,j) = -1; % blue
                elseif abs(x_final - 1) < 1e-3
                    basin_idx(i,j) = 1;  % red
                else
                    basin_idx(i,j) = 0;  % unknown
                end
            end
        end
    end

    % --- Plot subplot ---
    nexttile;
    contourf(x_vals, v_vals, basin_idx, 50, 'LineColor', 'none');
    colormap([1 0 0; 1 1 1; 0 0 1]); % Red, White, Blue
    title(sprintf('%s + %s', group1, group2), 'Interpreter', 'none');
    axis square;
    xlabel('$\bar{x}$','Interpreter','latex');
    ylabel('$\bar{v}$','Interpreter','latex');
    set(gca, 'FontSize', 14);

    filename = sprintf('pairwise_basins_1_%s_%s.mat', ...
        strrep(group1, ' ', '_'), ...
        strrep(group2, ' ', '_'));
    save(filename, 'basin_idx', 'x_vals', 'v_vals', 'group1', 'group2');
end

% Save figure
exportgraphics(gcf, 'pairwise_group_basins.pdf', 'ContentType', 'vector');


%%
clear all
close all


%% Parameters
x_op = 1;
nu = 1;
sigma = 3;
N = 39; % Number of nodes
beta = 1;
t_span = [0, 50];
grid_size = 270;
[x_vals, v_vals] = meshgrid(linspace(-3, 3, grid_size), linspace(-3, 3, grid_size));

groups = cell(1,39);
groups{19} = {'Head Tailor'};
groups{16} = {'Cutter'};
groups([25,26]) = {'Button Machiner'};
groups([29,33,39]) = {'Ironer'};
groups([30:32,34:38]) = {'Cotton Boy'};
groups([1:3,5:7,9,11:14,21,24,4,10,17,18,8,15,20,22:23,27:28]) = {'Tailor'};


% Define the group pairs
group_pairs = {
    {'Tailor', 'Button Machiner'};
    {'Tailor', 'Ironer'};
    {'Tailor', 'Cotton Boy'};
    {'Button Machiner', 'Ironer'};
    {'Button Machiner', 'Cotton Boy'};
    {'Ironer', 'Cotton Boy'};
    };

% --- Start plotting ---
figure;
tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

for pair_id = 1:length(group_pairs)
    group1 = group_pairs{pair_id}{1};
    group2 = group_pairs{pair_id}{2};

    % Find nodes in selected groups
    group1_idx = find(strcmp(groups, group1));
    group2_idx = find(strcmp(groups, group2));
    fixed_group_idx = unique([group1_idx, group2_idx]);

    % Head Tailor and Cutter fixed
    X_fixed = zeros(2*N, 1);
    X_fixed(N*2-1) = -1;

    % Load graph (pick kapfts1.dat)
    A = load('kapfts2.dat');
    D = diag(sum(A, 2));
    L = D - A;

    basin_idx = zeros(grid_size, grid_size);

    % Loop over grid
    for i = 1:grid_size
        for j = 1:grid_size
            X0 = X_fixed;

            for k = 1:N
                if ~any(k == fixed_group_idx)
                    continue;
                end
                X0(2*k - 1) = x_vals(i,j);
                X0(2*k)     = v_vals(i,j);
            end

            dynamics = @(t, X) system_dynamics(t, X, N, nu, beta, sigma, L);
            [~, X_sol] = ode45(dynamics, t_span, X0);

            x_final = mean(X_sol(end, 1:2:end));
            v_final = mean(X_sol(end, 2:2:end));

            if abs(v_final) < 1e-3
                if abs(x_final - (-1)) < 1e-3
                    basin_idx(i,j) = -1; % blue
                elseif abs(x_final - 1) < 1e-3
                    basin_idx(i,j) = 1;  % red
                else
                    basin_idx(i,j) = 0;  % unknown
                end
            end
        end
    end

    % --- Plot subplot ---
    nexttile;
    contourf(x_vals, v_vals, basin_idx, 50, 'LineColor', 'none');
    colormap([1 0 0; 1 1 1; 0 0 1]); % Red, White, Blue
    title(sprintf('%s + %s', group1, group2), 'Interpreter', 'none');
    axis square;
    xlabel('$\bar{x}$','Interpreter','latex');
    ylabel('$\bar{v}$','Interpreter','latex');
    set(gca, 'FontSize', 14);

    filename = sprintf('pairwise_basins_2_%s_%s.mat', ...
        strrep(group1, ' ', '_'), ...
        strrep(group2, ' ', '_'));
    % save(filename, 'basin_idx', 'x_vals', 'v_vals', 'group1', 'group2');
end




function dX = system_dynamics(~, X, N, nu, beta, sigma, L)
x = X(1:2:end);
v = X(2:2:end);
dx = v;
dv = -nu * v - 4 * beta * (x.^3 - x) - sigma * (L * x);
dX = zeros(2*N, 1);
dX(1:2:end) = dx;
dX(2:2:end) = dv;
end
