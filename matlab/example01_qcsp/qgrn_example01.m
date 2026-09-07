% QGRN_EXAMPLE01  MATLAB walkthrough of the QuantumGRN model.
%
%   A line-by-line counterpart of test/01_example.py, written to be read:
%   every quantum step goes through the MATLAB Support Package for Quantum
%   Computing (quantumCircuit / ryGate / cryGate / simulate / querystates)
%   rather than hand-rolled linear algebra, and the gradient is obtained with
%   the parameter-shift rule instead of hand-derived derivative circuits.
%
%   Run from this folder:
%       >> qgrn_example01
%
%   Requires: MATLAB Support Package for Quantum Computing.
%   Validate the conventions first with:  check_conventions
%
%   NOTE: this is the paper's order-SENSITIVE model. Genes are mapped to qubits
%   by activation ratio and the c-Ry layers do not commute across controls, so
%   the fitted network depends on that ordering. How much: over 41 orderings
%   the max edge spread is 0.0162 with no sign flips, but two edges straddle
%   the pruning threshold -- IRF4-AICDA appears in 24 of 41 orderings and
%   REL-RELA in 4 of 41. See README.md, "This model is order-sensitive".

clear; clc;
assert(exist('quantumCircuit', 'class') == 8 || exist('quantumCircuit', 'file') > 0, ...
    'MATLAB Support Package for Quantum Computing not found.');

dataFile = fullfile('..', '..', 'dataset', 'expr_matrix_pearsonresidual_7.txt');

%% 1. Data -> observed distribution ------------------------------------------
% Binarize each gene at 0, order genes by how often they are active, then count
% how often each of the 2^n activation patterns occurs. That histogram is the
% distribution the circuit is asked to reproduce.
fprintf('1. reading %s\n', dataFile);
[genes, pObs, act, ncells] = qgrn_read_data(dataFile, 0);
n = numel(genes);
fprintf('   %d genes x %d cells\n', n, ncells);
fprintf('   gene order (by activation): %s\n', strjoin(genes, ' '));
fprintf('   activation ratios         : %s\n', num2str(act', '%.3f '));

%% 2. Parameter initialisation ------------------------------------------------
% Diagonal angles encode each gene's activation ratio exactly:
%   Ry(theta)|0> = cos(theta/2)|0> + sin(theta/2)|1>, so P(active) = sin^2(theta/2)
%   and theta = 2*asin(sqrt(activation)) reproduces the marginal on qubit k.
% Off-diagonal (regulation) angles start at zero: no interactions yet.
theta0 = zeros(n);
theta0(1:n+1:end) = 2 * asin(sqrt(act));

pInit = qgrn_pout(theta0);
fprintf('\n2. initial circuit: %d gates, KL = %.6e\n', ...
        numel(qgrn_circuit(theta0).Gates), ...
        qgrn_loss(theta0, pInit, pObs, ncells, 1));

%% 3. Train --------------------------------------------------------------------
% Plain gradient descent. Every gradient entry is built from circuit
% simulations at shifted angles -- see qgrn_dpraw for why controlled
% rotations need a four-term shift rule rather than the usual two-term one.
fprintf('\n3. training (parameter-shift gradients)\n');
[theta, hist] = qgrn_train(theta0, pObs, ncells, ...
                           epochs = 100, lr = 0.1, trainEncoder = false);

%% 4. Results ------------------------------------------------------------------
pOut = qgrn_pout(theta);
fprintf('\n4. final KL = %.6e   sum sq. error = %.6e\n', hist.loss(end), hist.err(end));

disp(' ');
disp('   fitted theta (row = control, col = target):');
disp(array2table(round(theta, 4), 'VariableNames', genes, 'RowNames', genes));

%% 5. Figures ------------------------------------------------------------------
fprintf('5. writing figures\n');

fig = figure('Color', 'w', 'Position', [100 100 900 380]);
tiledlayout(1, 2, 'TileSpacing', 'compact');

nexttile;
keep = find(pObs > 0.01 | pOut > 0.01);
bar([pObs(keep) pOut(keep)]);
legend({'p^{obs}', 'p^{out}'}, 'Location', 'northeast'); box off;
set(gca, 'XTick', 1:numel(keep), 'XTickLabel', ...
    compose('%s', dec2bin(keep-1, n)), 'XTickLabelRotation', 90);
ylabel('probability'); title('observed vs fitted');

nexttile;
semilogy(hist.loss, 'LineWidth', 1.5); hold on;
semilogy(hist.err, 'LineWidth', 1.5); box off;
xlabel('epoch'); legend({'KL divergence', 'sq. error'});
title('training');
exportgraphics(fig, '01_comparison_matlab.pdf', 'ContentType', 'vector');
fprintf('  wrote 01_comparison_matlab.pdf\n');

qgrn_draw_network(genes, theta, 1, '01_qgrn_network_matlab.pdf');

%% 6. Export -------------------------------------------------------------------
% The GRN itself: one row per gene pair with its fitted angle.
[A, B] = ndgrid(1:n, 1:n);
m = A < B;
T = table(genes(A(m))', genes(B(m))', theta(m), ...
          'VariableNames', {'gene1', 'gene2', 'theta'});
T = sortrows(T, 'theta', 'descend', 'ComparisonMethod', 'abs');
writetable(T, '01_qgrn_edges_matlab.csv');
fprintf('  wrote 01_qgrn_edges_matlab.csv\n');
disp(' '); disp('   strongest edges:'); disp(head(T, 8));

%% 7. Cross-check against the Python package -----------------------------------
% ref_python_theta_final.csv was produced by test/01_example.py under the same
% settings (100 epochs, lr 0.1, encoder frozen). The two implementations should
% land on the same network.
if isfile('ref_python_theta_final.csv')
    Tpy = readmatrix('ref_python_theta_final.csv');
    d = max(abs(theta(:) - Tpy(:)));
    fprintf('\n7. max|theta_matlab - theta_python| = %.3e\n', d);
    if d < 1e-8
        fprintf('   -> MATLAB and Python agree to machine precision.\n');
    else
        fprintf('   -> divergence larger than expected; check settings.\n');
    end
end
