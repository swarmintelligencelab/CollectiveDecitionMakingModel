%% heatmap_c_sigma_Omega1_with_theory.m
% x-axis: c (sublevel size)
% y-axis: sigma (coupling strength)
% color : p_escape(c,sigma) = proportion of trials that ever leave the "well"
% starting from Omega_1(c) (k=1).
%
% Outside theorem-valid region:
%  - If P0(sigma) is not PD, we PROJECT it to SPD for sampling (so Omega is still compact).
%  - If r_k(sigma) is not positive, we fall back to separatrix radius |x_k^* - x_mid|.
%
% Model (networked, your sign convention):
%   xdot = v
%   vdot = -nu v - 4 beta (x.^3 - xop^2 x) - alpha*1 - sigma L x

clear; clc; close all; rng(1);

%% ------------ Parameters (edit) ------------
n     = 10;
beta  = 1.0;
xop   = 1.0;
alpha = 0.0;

k     = 1;        % focus on Omega_1 around x_1^*
tEnd  = 80;

Nc    = 200;       % # c points
Ns    = 200;       % # sigma points
Nmc   = 50;       % MC per (c,sigma)

c_vec     = linspace(0.01, 15.0, Nc);
sigma_vec = linspace(0.5, 5.0, Ns);

topology  = 'alltoall';   % 'alltoall','star','path','tree','nearest'
odeOpts = odeset('RelTol',1e-6,'AbsTol',1e-8);

epsSPD  = 1e-8;     % eigenvalue floor for SPD projection
use_r99 = true;     % use r=0.99*r_k in c_th
%% ------------------------------------------

%% Graph / Laplacian (fixed while sigma varies)
L = make_laplacian(n, topology);
A = [0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 0 1;
     1 0 0 0 0 0 0 0 0 0];
N = n;
for i = 1:N
    neighbors = mod((i + [1 2]) - 1, N) + 1;  % next 2 nodes (wrap around)
    A(i, neighbors) = 1;
end

Dout = diag(sum(A,2)); % sum across columns of each row = out-degree
L = Dout - A;

lam = sort(eig(L+L'));
lambda2 = lam(2);
lambdan = lam(end);

fprintf('Topology: %s | lambda2=%.4f, lambdan=%.4f\n', topology, lambda2, lambdan);

sigma_crit = (4*beta*xop^2)/lambda2;  % cond: sigma*lambda2 > 4 beta xop^2


alpha_c = (8*beta*xop^3)/(3*sqrt(3));

if abs(alpha) >= alpha_c
    warning('Not bistable: |alpha|=%.4g >= alpha_c=%.4g. Setting entire map to NaN.', abs(alpha), alpha_c);
    % if you already allocated your arrays:
    Pescape(:) = NaN;
    c_th_curve(:) = NaN;
    % you can still plot (it will be empty/NaN) and return:
    % return;
end


%% Intrinsic equilibria (stable) and separatrix (independent of sigma)
[xneg, xpos, xmid] = intrinsic_equilibria(alpha, beta, xop);

% Handle tie (alpha ~ 0 gives |xneg|=|xpos|)
if abs(abs(xneg)-abs(xpos)) < 1e-10
    xStars = [xneg, xpos];   % deterministic: x1* = negative, x2* = positive
else
    xStars = sort_by_abs([xneg, xpos]);  % ordered by |.| increasing
end

xk = xStars(k);
ak = abs(xk);

fprintf('alpha=%.3f: x1*=%.6f, x2*=%.6f | using k=%d, xk=%.6f\n', alpha, xStars(1), xStars(2), k, xk);

if ~isfinite(xmid)
    warning('Could not find separatrix root x_mid; fallback tube will be disabled.');
end

r_sep = abs(xk - xmid);  % intrinsic separatrix distance (always meaningful if bistable)

%% Preallocate
Pescape    = nan(Ns, Nc);     % escape probability heatmap
c_th_curve = nan(Ns, 1);      % theoretical c_th(sigma) (only where valid)

if isempty(gcp('nocreate'))
    parpool(6);   % uses default number of workers
end
sig_max = sqrt(max(abs(eig(L*L'))));
%% Main sweep
parfor is = 1:Ns
    Progress = is/Ns*100
    sigma = sigma_vec(is);
    ell2 = 4*beta*xop^2 - sigma*lambda2;

    nu = 1.1*sigma*sig_max/2/abs(ell2)*(1+abs(ell2));


    % Theorem-built P0(sigma)
    [P0_raw, lambda_x_raw] = build_P0(beta, xop, sigma, nu, lambda2, lambdan, sig_max);

    % --- SPD check for 2x2 matrix: P0_raw = [p1 p2; p2 p3]
    P0s = (P0_raw + P0_raw')/2;     % symmetrize (just in case)
    p1 = P0s(1,1); p2 = P0s(1,2); p3 = P0s(2,2);

    isPD = (p1 > 0) && (p3 > 0) && ((p1*p3 - p2^2) > 0);

    if isPD
        P0_use = P0s;
    else
        P0_use = eye(2);           % <-- your choice: identity fallback
    end

    % Tube radius for escape test:
    rk_raw = compute_rk_from_theorem(ak, beta, xop, sigma, nu, lambda2, lambdan);

    if isfinite(rk_raw) && rk_raw > 0
        r_tube = rk_raw;
    else
        r_tube = r_sep;   % fallback: separatrix-based "well" radius
    end

    % Theoretical boundary c_th(sigma) only where theorem makes sense
    if (sigma*lambda2 > 4*beta*xop^2) && isfinite(lambda_x_raw) && lambda_x_raw > 0 && isfinite(rk_raw) && rk_raw > 0
        r_use = rk_raw;
        if use_r99, r_use = 0.99*r_use; end
        c_th_curve(is) = lambda_x_raw * (r_use^2);
    else
        c_th_curve(is) = NaN;
    end

    % Pack params for ODE
    % Pack params for ODE (LOCAL to worker)
params = struct( ...
    'n', n, ...
    'beta', beta, ...
    'xop', xop, ...
    'sigma', sigma, ...
    'nu', nu, ...
    'alpha', alpha, ...
    'L', L );


    % Sweep c
    for ic = 1:Nc
        c = c_vec(ic);
        Pescape(is, ic) = escape_rate_MC(c, xk, r_tube, n, P0_use, params, tEnd, odeOpts, Nmc);
    end
end

%% Plot heatmap + overlays
figure('Color','w'); hold on;
imagesc(c_vec, sigma_vec, Pescape);
set(gca,'YDir','normal');
xlabel('c (sublevel size of \Omega_1)');
ylabel('\sigma (coupling strength)');
title(sprintf('Escape probability from \\Omega_1 (k=1), n=%d, %s, \\alpha=%.3f', n, topology, alpha));
cb = colorbar; cb.Label.String = 'p_{escape}(c,\sigma)';
grid on;

% Overlay theorem boundaries
plot(c_th_curve, sigma_vec, 'w', 'LineWidth', 2);
hold on;
plot([min(c_vec) max(c_vec)], [sigma_crit sigma_crit], 'r', 'LineWidth', 2);
text(max(c_vec), sigma_crit, '  \sigma_{crit}', 'Color','r', ...
     'VerticalAlignment','bottom','HorizontalAlignment','left');


%% ===================== Local functions =====================

function rate = escape_rate_MC(c, xStar, rTube, n, P0, params, tEnd, odeOpts, Nmc)
    if ~isfinite(rTube) || rTube <= 0
        % If the tube radius is nonsensical, count everything as escape
        rate = 1.0;
        return;
    end

    leaveCount = 0;

    for m = 1:Nmc
        X0 = sample_X0_in_Omega(n, xStar, P0, c);

        try
            [~,X] = ode45(@(t,x) dynamics_ode(t,x,params), [0 tEnd], X0, odeOpts);
            xTraj = X(:,1:n);
            maxDev = max(abs(xTraj - xStar), [], 2);
            if any(maxDev > rTube)
                leaveCount = leaveCount + 1;
            end
        catch
            % If the solver fails (stiff blow-up etc.), treat as escape
            leaveCount = leaveCount + 1;
        end
    end

    rate = leaveCount / Nmc;
end

function dX = dynamics_ode(~, X, p)
    n = p.n;
    x = X(1:n);
    v = X(n+1:end);

    xdot = v;
    vdot = -p.nu*v ...
           - 4*p.beta.*(x.^3 - p.xop^2*x) ...
           - p.alpha ...
           - p.sigma*(p.L*x);

    dX = [xdot; vdot];
end

function X0 = sample_X0_in_Omega(n, xStar, P0, c)
    % Sample approximately uniformly in { sum_i delta_i' P0 delta_i <= c }
    if c <= 0
        X0 = [xStar*ones(n,1); zeros(n,1)];
        return;
    end

    U = chol(P0,'lower');     % P0 = U*U'
    dim = 2*n;

    z = randn(dim,1); z = z/norm(z);
    rad = (rand()^(1/dim))*sqrt(c);
    y = rad*z;

    K = kron(eye(n), U);
    delta = K \ y;

    dx = delta(1:2:end);
    dv = delta(2:2:end);

    X0 = [xStar + dx; dv];
end

function xs = sort_by_abs(x)
    [~,idx] = sort(abs(x),'ascend');
    xs = x(idx);
end

function [xneg, xpos, xmid] = intrinsic_equilibria(alpha, beta, xop)
    % 0 = -4 beta (x^3 - xop^2 x) - alpha
    % -> 4 beta x^3 - 4 beta xop^2 x + alpha = 0
    coeff = [4*beta, 0, -4*beta*xop^2, +alpha];
    rts = roots(coeff);
    rts = rts(abs(imag(rts))<1e-8);
    rts = sort(real(rts));

    % stable if |x| > xop/sqrt(3)
    stab = rts(abs(rts) > xop/sqrt(3) + 1e-10);

    if numel(stab) < 2
        xneg = NaN; xpos = NaN; xmid = NaN;
        warning('Not bistable for this alpha (need |alpha| < alpha_c).');
        return;
    end

    xneg = stab(1);
    xpos = stab(end);

    if numel(rts) == 3
        xmid = rts(2);
    else
        xmid = NaN;
    end
end

function [P0, lambda_x] = build_P0(beta, xop, sigma, nu, lambda2, lambdan, sig_max)
    ell2 = 4*beta*xop^2 - sigma*lambda2;
    elln = 4*beta*xop^2 - sigma*lambdan;

    p3 = 1.0;
    q = p3^2*abs(ell2)*(2*nu-sigma*sig_max)^2;
    d = p3^2*sigma*sig_max*(2*nu-sigma*sig_max);
    p2 = p3*(q+d)/(2*q)*(2*nu-sigma*sig_max);
    p1 = nu*p2 - p3*(ell2 + elln)/2;

    P0 = [p1 p2; p2 p3];
    lambda_x = p1 - (p2^2)/p3;
end

function rk = compute_rk_from_theorem(ak, beta, xop, sigma, nu, lambda2, lambdan)
    dnu = xop^2 + (sigma^2*(lambda2 + lambdan)^2)/(16*beta*nu^2);
    disc = dnu - (3/4)*ak^2;
    rk = (3/2)*ak - sqrt(max(disc,0));
end


function [L,A] = make_laplacian(n, topology, varargin)
    topology = lower(string(topology));
    k = 1;
    treeMode = "binary";

    for ii = 1:2:numel(varargin)
        name = lower(string(varargin{ii}));
        val  = varargin{ii+1};
        switch name
            case "k"
                k = val;
            case "mode"
                treeMode = lower(string(val));
            otherwise
                error('Unknown option "%s".', name);
        end
    end

    A = zeros(n);

    switch topology
        case {"alltoall","complete","fullyconnected"}
            A = ones(n) - eye(n);

        case "star"
            for i = 2:n, A(1,i)=1; A(i,1)=1; end

        case {"path","line"}
            for i = 1:n-1, A(i,i+1)=1; A(i+1,i)=1; end

        case {"nearest","nearestneighbors","knearest","ring"}
            k = max(1, min(k, floor((n-1)/2)));
            for i = 1:n
                for d = 1:k
                    j1 = mod(i-1 + d, n) + 1;
                    j2 = mod(i-1 - d, n) + 1;
                    A(i,j1)=1; A(j1,i)=1;
                    A(i,j2)=1; A(j2,i)=1;
                end
            end

        case {"tree","binarytree"}
            switch treeMode
                case "binary"
                    for i = 1:n
                        c1 = 2*i; c2 = 2*i + 1;
                        if c1 <= n, A(i,c1)=1; A(c1,i)=1; end
                        if c2 <= n, A(i,c2)=1; A(c2,i)=1; end
                    end
                case "random"
                    for i = 2:n
                        j = randi(i-1);
                        A(i,j)=1; A(j,i)=1;
                    end
                otherwise
                    error('Unknown tree mode "%s".', treeMode);
            end

        otherwise
            error('Unknown topology "%s".', topology);
    end

    deg = sum(A,2);
    L = diag(deg) - A;
end