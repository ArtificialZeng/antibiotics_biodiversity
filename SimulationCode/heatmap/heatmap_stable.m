% Heatmap of presence/absence of interior attractor (fixed point or limit cycle)
% Community: {PD, SP, RS}, KP=40, KD=4, g1=1, costs=0.
% Scans g2 (x-axis) and g3 (y-axis) in [0.01, 10] (log-spaced).

clear; clc;

%% Parameters
KP   = 40;
KD   = 4;
g1   = 1;

% Grids (log-spaced)
Ngrid   = 100;                       % adjust resolution if desired
g2vals  = logspace(-1, 1, Ngrid);   % 0.01 -> 10
g3vals  = logspace(-1, 1, Ngrid);
% g2vals  = exp(linspace(0.01, log(10), Ngrid));   % 0.01 -> 10
% g3vals  = exp(linspace(0.01, log(10), Ngrid));

% Simulation settings
num_inits   = 10;      % random initial conditions per (g2,g3)
T_total     = 3000;    % total iterations
T_burn      = 1500;    % burn-in
L_window    = 200;     % analyze last L_window points
eps_vertex  = 1e-3;    % threshold to call a vertex attractor
tol_fp      = 1e-6;    % fixed-point stability tolerance (std threshold)
tol_cycle   = 5e-4;    % cycle detection tolerance
max_period  = 6;       % search cycles up to period 6

%% Matrices (match the R code)
S    = [0 0 0; 1 0 0; 0 1 0];   % susceptibility
Pmat = [1 0 0; 0 1 0; 0 0 0];   % producers
Dmat = [0 1 0; 0 0 0; 0 0 0];   % degraders
% Rmat = [0 0 0; 0 0 0; 1 0 0];  % intrinsic resistance (unused when costs=0)
Rmat = [0 0 0; 0 0 0; 0 0 0];   % intrinsic resistance (unused when costs=0)

% S    = [0 0 1; 1 0 0; 0 1 0];   % susceptibility
% Pmat = [1 0 0; 0 1 0; 0 0 1];   % producers
% Dmat = [0 1 0; 0 0 1; 1 0 0];   % degraders
% Rmat = [0 0 0; 0 0 0; 0 0 0];   % intrinsic resistance (unused when costs=0)


% Pack env so functions don’t depend on script scope
env = struct('S',S,'Pmat',Pmat,'Dmat',Dmat,'Rmat',Rmat,'KP',KP,'KD',KD);

%% Scan grid
classMap = zeros(Ngrid, Ngrid); % 0 none, 1 limit cycle, 2 fixed point

for iy = 1:Ngrid
    g3 = g3vals(iy);
    for ix = 1:Ngrid
        g2 = g2vals(ix);
        g = [g1, g2, g3];

        found_fixed = false;
        found_cycle = false;

        for trial = 1:num_inits
            X = rand_simplex(1); % 1x3
            % iterate
            traj = zeros(T_total, 3);
            traj(1,:) = X;
            for t = 2:T_total
                traj(t,:) = step_map(traj(t-1,:), g, env);
            end

            tail = traj((T_total - T_burn + 1):end, :);
            wnd  = tail(max(1,size(tail,1)-L_window+1):end,:);

            if is_vertex_attractor(tail, eps_vertex)
                % keep trying other inits
            else
                if is_interior_fixed_point(wnd, tol_fp)
                    found_fixed = true; break;
                else
                    [iscycle, ~] = is_small_limit_cycle(wnd, tol_cycle, max_period);
                    if iscycle
                        found_cycle = true; % keep scanning for fixed point just in case
                    end
                end
            end
        end

        if found_fixed
            classMap(iy, ix) = 2;
        elseif found_cycle
            classMap(iy, ix) = 1;
        else
            classMap(iy, ix) = 0;
        end
    end
    fprintf('Row %d/%d done.\n', iy, Ngrid);
end

%% Plot
figure('Color','w');
imagesc(log10(g2vals), log10(g3vals), classMap); axis xy;
colormap([0.95 0.95 0.95;  % 0: none  (light grey)
          0.30 0.65 1.00;  % 1: cycle (blue)
          0.10 0.70 0.10]);% 2: fixed (green)
caxis([0 2]);
colorbar('Ticks',[0,1,2],'TickLabels',{'none','cycle','fixed point'});
xlabel('log_{10}(g_2)'); ylabel('log_{10}(g_3)');
% xlabel('g_2'); ylabel('g_3');
title('Interior attractor map for {PD, SP, RS} (KP=40, KD=4, g1=1)');
xticks(log10([0.01 0.1 1 10])); xticklabels({'0.01','0.1','1','10'});
yticks(log10([0.01 0.1 1 10])); yticklabels({'0.01','0.1','1','10'});
% xticks(log([0.01 0.1 1 10])); xticklabels({'0.01','0.1','1','10'});
% yticks(log([0.01 0.1 1 10])); yticklabels({'0.01','0.1','1','10'});
grid on; set(gca,'Layer','top');

%% ---- Local functions (stateless; take inputs explicitly) ----
function Xnext = step_map(X, g, env)
    % costs are zero, so r = g
    r = g(:)'; % row 1x3

    % Community-level effects
    Aj = X * env.Dmat;  % 1x3 (amount of degraders acting on each antibiotic j)
    Bj = X * env.Pmat;  % 1x3 (amount of producers of each antibiotic j)

    % Susceptibility & survival
    pij = zeros(3,3);
    for i = 1:3
        for j = 1:3
            pij(i,j) = env.S(i,j) * exp(-env.KD * Aj(j)) * (1 - exp(-env.KP * Bj(j)));
        end
    end
    surv = prod(1 - pij, 2)';     % 1x3
    fi = r .* surv;               % 1x3
    fi(fi < 0) = 0;

    denom = sum(X .* fi);
    if ~isfinite(denom) || denom <= 1e-12
        Xnext = X;
    else
        Xnext = (X .* fi) / denom;
    end

    % numerical safeguard to keep on simplex
    Xnext = max(Xnext, 0);
    Xnext = Xnext / sum(Xnext);
end

function X = rand_simplex(n)
    % draw n random points from 3-simplex (Dirichlet(1,1,1))
    E = exprnd(1, n, 3);
    X = E ./ sum(E, 2);
end

function tf = is_vertex_attractor(traj, eps_vertex)
    % traj: Tx3 (post burn-in). Vertex if some coord ~1 for all last steps.
    xmax = max(traj, [], 2);
    tf = all(xmax > (1 - eps_vertex));
end

function tf = is_interior_fixed_point(traj, tol_fp)
    % small variance in last window AND all components bounded away from 0
    if size(traj,1) < 20
        tf = false; return;
    end
    if any(mean(traj(end-19:end,:), 1) < 1e-4) % quick boundary check
        tf = false; return;
    end
    s = std(traj, 0, 1);
    tf = max(s) < tol_fp;
end

function [tf, pbest] = is_small_limit_cycle(traj, tol_cycle, max_period)
    % look for small period cycles up to max_period via mean recurrence
    if size(traj,1) < (max_period+10)
        tf = false; pbest = NaN; return;
    end
    % reject near-vertex behavior
    if any(mean(traj(end-49:end,:), 1) < 1e-4)
        tf = false; pbest = NaN; return;
    end
    tf = false; pbest = NaN;
    for p = 2:max_period
        d = vecnorm(traj((p+1):end,:)-traj(1:end-p,:), 2, 2);
        if mean(d) < tol_cycle
            tf = true; pbest = p; return;
        end
    end
end


% g2        g3  
% ______    ______
% 4.5349    6.8926
% 6.8926    5.4623