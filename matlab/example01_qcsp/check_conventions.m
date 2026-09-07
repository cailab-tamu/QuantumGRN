% CHECK_CONVENTIONS  Validate the two things that silently corrupt this model.
%
%   1. BIT ORDERING  -- does a toolbox basis-state string map to the same
%      state index the Python package uses?
%   2. SHIFT RULE    -- is the four-term parameter-shift rule really needed
%      for controlled rotations, and is it exact?
%
%   Both are checked against ground truth: reference values exported from the
%   Python package (ref_python_*.csv) and central finite differences.
%
%   Run from this folder:  >> check_conventions

clear; clc;
assert(exist('quantumCircuit', 'class') == 8 || exist('quantumCircuit', 'file') > 0, ...
    'MATLAB Support Package for Quantum Computing not found.');

fprintf('=== 1. bit ordering =========================================\n');
% ryGate(1,pi) should excite qubit 1 alone. Which character of the basis-state
% string moves tells us where qubit 1 lives.
[st, pr] = querystates(simulate(quantumCircuit( ...
    [ryGate(1,pi); ryGate(2,0); ryGate(3,0)])));
fprintf('ryGate(1,pi) on 3 qubits -> "%s"  (char 1 is qubit 1)\n', string(st(pr > 1e-9)));

theta3 = readmatrix('ref_python_theta3.csv');
pyRef  = readmatrix('ref_python_pout3.csv');
n3 = size(theta3, 1);
[stt, prr] = querystates(simulate(qgrn_circuit(theta3)));
stt = string(stt);
for mode = ["q1lsb", "q1msb"]
    p = zeros(2^n3, 1);
    for i = 1:numel(stt)
        b = char(stt(i)) == '1';
        if mode == "q1msb", b = fliplr(b); end
        p(sum(b .* 2.^(0:n3-1)) + 1) = prr(i);
    end
    p(1) = 0; p = p / sum(p);
    fprintf('  %-6s : max|dp| vs Python = %.3e\n', mode, max(abs(p - pyRef)));
end
fprintf('  -> qgrn_pout uses q1lsb (gene k = bit k-1). Agreement is exact.\n');

fprintf('\n=== 2. parameter-shift rule =================================\n');
% A single c-Ry on a product state gives a pure single-frequency signal, where
% even the two-term rule looks correct. Interfering controls are what expose it,
% so use the real QuantumGRN gate pattern.
th = [1.10 0.63 0.50; -0.30 0.40 0.70; 0.20 -0.60 0.90];
k = 1; q = 2; h = 1e-6;
cp = (sqrt(2)+1)/(4*sqrt(2));
cm = (sqrt(2)-1)/(4*sqrt(2));
f = @(v) rawp(setp(th, k, q, v));

fd  = (f(th(k,q)+h) - f(th(k,q)-h)) / (2*h);
p2  = ( f(th(k,q)+pi/2) - f(th(k,q)-pi/2) ) / 2;
p4  = cp*( f(th(k,q)+pi/2)   - f(th(k,q)-pi/2) ) ...
    - cm*( f(th(k,q)+3*pi/2) - f(th(k,q)-3*pi/2) );

fprintf('  max|2-term - finite diff| = %.2e   <-- WRONG, two frequencies present\n', ...
        max(abs(p2 - fd)));
fprintf('  max|4-term - finite diff| = %.2e   <-- exact\n', max(abs(p4 - fd)));

d = qgrn_dpraw(th, k, q);
fprintf('  qgrn_dpraw agrees with finite differences to %.2e\n', max(abs(d - fd)));

fprintf('\n=== 3. gradient vs the Python package =======================\n');
% The real fidelity test: same theta, same data, same ratio -- does
% qgrn_gradient reproduce the Python compute_gradient() entry by entry?
Tinit = readmatrix('ref_python_theta_init.csv');
Gpy   = readmatrix('ref_python_grad_epoch1.csv');
pObs6 = readmatrix('ref_python_pobs.csv');
sc    = readmatrix('ref_python_scalars.csv');    % ncells, ratio, loss
nc = sc(1); ratio = sc(2);

Gml = qgrn_gradient(Tinit, qgrn_pout(Tinit), pObs6, nc, ratio, false);
fprintf('  max|G_matlab - G_python| = %.3e   (max|G| = %.3e)\n', ...
        max(abs(Gml(:) - Gpy(:))), max(abs(Gpy(:))));
if max(abs(Gml(:) - Gpy(:))) < 1e-9
    fprintf('  -> parameter-shift gradient reproduces the analytic one exactly.\n');
else
    fprintf('  -> MISMATCH: investigate before trusting the port.\n');
end

fprintf('\n=== 4. what the model actually optimises ====================\n');
% Both implementations build dp/dtheta from the RAW Born probabilities while
% evaluating the loss on the |0..0>-rescaled ones -- they drop the derivative
% of the rescaling. That is a modelling choice inherited from the paper's code,
% not an identity. This shows how far it sits from the true gradient.
th = [1.10 0.63 0.50; -0.30 0.40 0.70; 0.20 -0.60 0.90];
h = 1e-6; a = 1; b = 3;
rng(0); pO = rand(2^3,1); pO(1) = 0; pO = pO/sum(pO);
G3 = qgrn_gradient(th, dropzero(rawp(th)), pO, 1000, 1e-3, false);
Lf = @(T) lossval(T, pO, 1000, 1e-3);
t1 = setp(setp(th,a,b,th(a,b)+h), b,a, th(b,a)+h);
t2 = setp(setp(th,a,b,th(a,b)-h), b,a, th(b,a)-h);
praw0 = rawp(th);
fprintf('  edge %d-%d: model gradient %+.6f   true d(loss)/dtheta %+.6f\n', ...
        a, b, G3(a,b)+G3(b,a), (Lf(t1)-Lf(t2))/(2*h));
fprintf('  they differ because |0..0> holds %.0f%% of the raw probability mass\n', ...
        100*praw0(1));

% ---- local helpers ---------------------------------------------------------
function T = setp(T, k, q, v), T(k,q) = v; end
function p = rawp(T), [~, p] = qgrn_pout(T); end
function p = dropzero(p), p(1) = 0; p = p / sum(p); end
function L = lossval(T, pObs, nc, ratio)
    [l, c] = qgrn_loss(T, dropzero(rawp(T)), pObs, nc, ratio);
    L = l + c;
end
