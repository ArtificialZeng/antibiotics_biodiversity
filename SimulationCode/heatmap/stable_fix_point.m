% Find & test stability of the interior fixed point for given (g2,g3)
% Community {PD,SP,RS}, KP=40, KD=4, g1=1, costs=0 (same matrices as before)
clear; clc;

%% Parameters
KP = 40; KD = 4; g1 = 1;

% Matrices (match the R code you shared)
S    = [0 0 0; 1 0 0; 0 1 0];
Pmat = [1 0 0; 0 1 0; 0 0 0];
Dmat = [0 1 0; 0 0 0; 0 0 0];
Rmat = [0 0 0; 0 0 0; 1 0 0]; %#ok<NASGU> (unused when costs=0)
env = struct('S',S,'Pmat',Pmat,'Dmat',Dmat,'KP',KP,'KD',KD);

% Example pairs
pairs = [ 4.5349 6.8926
          6.8926 5.4623 ];

for k = 1:size(pairs,1)
    g2 = pairs(k,1); g3 = pairs(k,2);
    [ok, Xstar, reason] = analytic_fp(g1,g2,g3,KP,KD);
    fprintf('\n=== g2=%.4f, g3=%.4f ===\n', g2, g3);
    if ~ok
        fprintf('No feasible interior fixed point: %s\n', reason);
        continue;
    end
    fprintf('Candidate interior FP: x* = [%.6f, %.6f, %.6f]\n', Xstar);

    % Stability via 2D Jacobian on the simplex (x3 = 1 - x1 - x2)
    J = jac2D(@(x1,x2) step12([x1,x2,1-x1-x2], [g1,g2,g3], env), Xstar(1), Xstar(2));
    ev = eig(J);
    sr = max(abs(ev));
    fprintf('Eigenvalues (2D): %.6f%+.6fi,  %.6f%+.6fi  | spectral radius = %.6f\n', ...
        real(ev(1)), imag(ev(1)), real(ev(2)), imag(ev(2)), sr);
    if sr < 1
        fprintf('=> STABLE interior fixed point ✅\n');
    else
        fprintf('=> UNSTABLE interior fixed point ❌\n');
    end
end

%% ---------- Functions ----------
function [ok, Xstar, reason] = analytic_fp(g1,g2,g3,KP,KD)
% Analytic interior FP using f1=f2=f3 with costs=0 and the given community
% f1 = g1
% f2 = g2 * exp(-KP*x1)
% f3 = g3 * (1 - exp(-KD*x1)*(1 - exp(-KP*x2)))

    ok = false; Xstar = [NaN,NaN,NaN]; reason = '';

    if g2 <= g1
        reason = 'Need g2>g1 for x1*>0 (from f1=f2).'; return;
    end
    x1 = (1/KP)*log(g2/g1);              % from g1 = g2*e^{-KP x1}
    a  = exp(-KD*x1);                    % a in (0,1)
    num = (g1/g3 + a - 1);
    if ~(g3 > g1)
        reason = 'Need g3>g1 for f1=f3 feasibility.'; return;
    end
    if num <= 0 || num >= a
        reason = 'b not in (0,1): violates 0 < (g1/g3 + a - 1) < a.'; return;
    end
    b  = num / a;                         % b = e^{-KP x2}
    x2 = -(1/KP)*log(b);
    x3 = 1 - x1 - x2;

    if x1<=0 || x2<=0 || x3<=0
        reason = 'One or more coordinates not strictly positive.'; return;
    end

    Xstar = [x1,x2,x3];
    ok = true;
end

function Xnext = step_map(X, g, env)
    % costs=0 => r=g
    r = g(:)';                        % 1x3
    Aj = X * env.Dmat;                % 1x3
    Bj = X * env.Pmat;                % 1x3

    pij = zeros(3,3);
    for i = 1:3
        for j = 1:3
            pij(i,j) = env.S(i,j) * exp(-env.KD * Aj(j)) * (1 - exp(-env.KP * Bj(j)));
        end
    end
    surv = prod(1 - pij, 2)';         % 1x3
    fi = r .* surv; fi(fi<0) = 0;     % 1x3

    denom = sum(X .* fi);
    if ~isfinite(denom) || denom<=1e-12
        Xnext = X;
    else
        Xnext = (X .* fi) / denom;
    end
    Xnext = max(Xnext,0); Xnext = Xnext / sum(Xnext); % keep on simplex
end

function x12next = step12(X, g, env)
    Xn = step_map(X, g, env);
    x12next = Xn(1:2).';   % return as 2x1 column
end

function J = jac2D(F, x1, x2)
    % Numerical 2x2 Jacobian of F = [x1'; x2'] at (x1,x2)
    h = 1e-7;
    f0  = F(x1, x2); %#ok<NASGU>
    f1p = F(x1+h, x2);
    f1m = F(x1-h, x2);
    f2p = F(x1, x2+h);
    f2m = F(x1, x2-h);

    dFdx1 = (f1p - f1m) / (2*h);  % 2x1
    dFdx2 = (f2p - f2m) / (2*h);  % 2x1

    J = [dFdx1, dFdx2];           % 2x2 matrix
end
