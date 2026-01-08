clear; clc; close all; rng(1);

%% ------------ Parameters (edit) ------------
n     = 10;
beta  = 1.0;
xop   = 1.0;
sigma = 0.40;
nu    = 2.0;

alpha = -0.2;  
tEnd  = 80;

Nmc   = 10000;       % Monte Carlo per c
Nc    = 100;       % number of c points
c_vec = linspace(0, 8, Nc);

% choose topology
topology = 'alltoall';      % 'alltoall','star','path','tree','nearest'
% [optional] for nearest: kNN = 2;  L = make_laplacian(n,'nearest','k',kNN);
%% ------------------------------------------

%% Graph / Laplacian
L = make_laplacian(n, topology);   % first output is L (works with your call style)

lam = sort(eig(L));
lambda2 = lam(2);
lambdan = lam(end);

fprintf('Topology: %s | lambda2=%.4f, lambdan=%.4f\n', topology, lambda2, lambdan);

%% Build P0 from theorem formulas
ell2 = 4*beta*xop^2 - sigma*lambda2;
elln = 4*beta*xop^2 - sigma*lambdan;

p3 = 1.0;
p2 = (nu*p3)/2;
p1 = nu*p2 - p3*(ell2 + elln)/2;

P0 = [p1 p2; p2 p3];
lambda_x = p1 - (p2^2)/p3;

if min(eig(P0)) <= 0
    warning('P0 not PD with these parameters. Adjust (nu,sigma,beta,xop).');
end

%% Intrinsic equilibria (stable) for THIS alpha (your sign: -alpha)
[xneg, xpos, ~] = intrinsic_equilibria(alpha, beta, xop);
xStars = sort_by_abs([xneg, xpos]);   % x1*, x2* ordered by |.| increasing
x1 = xStars(1);
x2 = xStars(2);

fprintf('alpha=%.3f: x1*=%.6f, x2*=%.6f (ordered by |.|)\n', alpha, x1, x2);

%% Theoretical radii r1,r2 and thresholds c_th1,c_th2 (optional overlays)
dnu = xop^2 + (sigma^2*(lambda2 + lambdan)^2)/(16*beta*nu^2);

rk_vec   = zeros(2,1);
c_th_vec = zeros(2,1);

for kk = 1:2
    ak = abs(xStars(kk));
    disc = dnu - (3/4)*ak^2;

    % r_k = (3/2)a_k - sqrt(dnu - (3/4)a_k^2)
    % clamp for numerics (to avoid complex if disc<0)
    rk = (3/2)*ak - sqrt(max(disc,0));
    rk_vec(kk) = rk;

    r = 0.99*rk;                 % strict tube
    c_th_vec(kk) = lambda_x*r^2;  % c_th = lambda_x r^2
end

r1 = rk_vec(1);  r2 = rk_vec(2);
c_th1 = c_th_vec(1); c_th2 = c_th_vec(2);
fprintf('r1=%.4f, c_th1=%.4f | r2=%.4f, c_th2=%.4f\n', r1, c_th1, r2, c_th2);

%% Monte Carlo sweep for both wells
odeOpts = odeset('RelTol',1e-6,'AbsTol',1e-8);

params.n=n; params.beta=beta; params.xop=xop;
params.sigma=sigma; params.nu=nu; params.alpha=alpha; params.L=L;

leave_rate_1 = nan(size(c_vec));
leave_rate_2 = nan(size(c_vec));

% bands (optional)
lo1 = nan(size(c_vec)); hi1 = nan(size(c_vec));
lo2 = nan(size(c_vec)); hi2 = nan(size(c_vec));
parpool(6)
parfor ic = 1:numel(c_vec)
    c = c_vec(ic);

    [p1hat, lo, hi] = leave_rate_for_well(c, x1, r1, n, P0, params, tEnd, odeOpts, Nmc);
    leave_rate_1(ic) = p1hat; lo1(ic) = lo; hi1(ic) = hi;

    [p2hat, lo, hi] = leave_rate_for_well(c, x2, r2, n, P0, params, tEnd, odeOpts, Nmc);
    leave_rate_2(ic) = p2hat; lo2(ic) = lo; hi2(ic) = hi;
end

%% Plot
figure('Color','w'); hold on; grid on;


hB = plot(c_vec, leave_rate_1, 'b-', 'LineWidth', 2);
hR = plot(c_vec, leave_rate_2, 'r-', 'LineWidth', 2);

ylim([0 1.05]);
xlabel('c (sublevel size of \Omega_k)');
ylabel('Proportion of trials escaping \Omega_k');

% theoretical thresholds (if positive/finite)
if isfinite(c_th1) && c_th1>0, xline(c_th1, 'b--', 'LineWidth', 1.8, 'Label', 'c_{th,1}'); end
if isfinite(c_th2) && c_th2>0, xline(c_th2, 'r--', 'LineWidth', 1.8, 'Label', 'c_{th,2}'); end

% if showBands
%     legend([hB hR], {'\Omega_1 around x_1^*','\Omega_2 around x_2^*'}, 'Location','best');
% else
%     legend([hB hR], {'\Omega_1 around x_1^*','\Omega_2 around x_2^*'}, 'Location','best');
% end

title(sprintf('Escape vs c (n=%d, %s, \\alpha=%.3f)', n, topology, alpha));
delete(gcp('nocreate'))
%% ===================== Local functions =====================

function [rate, lo, hi] = leave_rate_for_well(c, xStar, rk, n, P0, params, tEnd, odeOpts, Nmc)
% Returns:
%   rate = (# trials that ever violate tube) / Nmc
%   lo,hi = ~95% normal-approx CI for binomial proportion (clamped to [0,1])

    if ~isfinite(rk) || rk <= 0
        rate = NaN; lo = NaN; hi = NaN;
        return;
    end

    leaveCount = 0;

    for m = 1:Nmc
        X0 = sample_X0_in_Omega(n, xStar, P0, c);

        [~,X] = ode45(@(t,x) dynamics_ode(t,x,params), [0 tEnd], X0, odeOpts);
        xTraj = X(:,1:n);

        maxDev = max(abs(xTraj - xStar), [], 2);
        if any(maxDev > rk)
            leaveCount = leaveCount + 1;
        end
    end

    rate = leaveCount / Nmc;

    % simple 95% band (binomial SE)
    z  = 1.96;
    se = sqrt(max(rate*(1-rate),0) / Nmc);
    lo = max(0, rate - z*se);
    hi = min(1, rate + z*se);
end

function dX = dynamics_ode(~, X, p)
    n = p.n;
    x = X(1:n);
    v = X(n+1:end);

    xdot = v;

    % your sign convention: "... - alpha"
    vdot = -p.nu*v ...
           - 4*p.beta.*(x.^3 - p.xop^2*x) ...
           - p.alpha ...
           - p.sigma*(p.L*x);

    dX = [xdot; vdot];
end

function X0 = sample_X0_in_Omega(n, xStar, P0, c)
% Sample approximately uniformly in { sum_i delta_i' P0 delta_i <= c }
% where delta_i = [x_i - xStar; v_i].

    if c <= 0
        X0 = [xStar*ones(n,1); zeros(n,1)];
        return;
    end

    U = chol(P0,'lower');     % P0 = U*U'
    dim = 2*n;

    % uniform-ish in 2n ball radius sqrt(c)
    z = randn(dim,1); z = z/norm(z);
    rad = (rand()^(1/dim))*sqrt(c);
    y = rad*z;

    % delta = (I kron U)^{-1} y
    K = kron(eye(n), U);
    delta = K \ y;

    dx = delta(1:2:end);
    dv = delta(2:2:end);

    x0 = xStar + dx;
    v0 = dv;

    X0 = [x0; v0];
end

function xs = sort_by_abs(x)
    [~,idx] = sort((x),'ascend');
    xs = x(idx);
end

function [xneg, xpos, xmid] = intrinsic_equilibria(alpha, beta, xop)
% Equilibria of intrinsic dynamics with YOUR sign:
%   0 = -4 beta (x^3 - xop^2 x) - alpha
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

function [L,A] = make_laplacian(n, topology, varargin)
%MAKE_LAPLACIAN  Build Laplacian for common undirected networks.
% Usage:
%   L = make_laplacian(n,'alltoall')
%   L = make_laplacian(n,'nearest','k',2)
%   [L,A] = make_laplacian(...)

    topology = lower(string(topology));
    k = 1;                  % for nearest-neighbors
    treeMode = "binary";    % for tree

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
            for i = 2:n
                A(1,i) = 1; A(i,1) = 1;
            end

        case {"path","line"}
            for i = 1:n-1
                A(i,i+1) = 1; A(i+1,i) = 1;
            end

        case {"nearest","nearestneighbors","knearest","ring"}
            k = max(1, min(k, floor((n-1)/2)));
            for i = 1:n
                for d = 1:k
                    j1 = mod(i-1 + d, n) + 1;
                    j2 = mod(i-1 - d, n) + 1;
                    A(i,j1) = 1; A(j1,i) = 1;
                    A(i,j2) = 1; A(j2,i) = 1;
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
                    error('Unknown tree mode "%s". Use "binary" or "random".', treeMode);
            end

        otherwise
            error('Unknown topology "%s".', topology);
    end

    deg = sum(A,2);
    L = diag(deg) - A;
end