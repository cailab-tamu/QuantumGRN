% RUN_COMPARISON  MATLAB/Octave port of run_comparison.py
%   Model A (original qscGRN, order-sensitive) vs Model B (Ising Born
%   machine, order-free) on the 6-gene LCL data. See ../Provenance.md.
%   Run from test/order_free_rzz/matlab/.  Runtime ~ a few minutes.
clear; clc;
datafile = fullfile('..','..','..','dataset','expr_matrix_pearsonresidual_7.txt');
[X, genes] = read_expr(datafile);
Xb = double(X > 0); m = size(Xb,2);
act = mean(Xb, 2);
[~, ord] = sort(act, 'descend');                  % paper ordering (by activation ratio)
names = genes(ord); Xb = Xb(ord,:); act = act(ord);
n = numel(names); Z = zz_signs(n);
fprintf('activation order: %s\n', strjoin(names', ' '));
pobs = build_pobs(Xb);
theta0 = diag(2*asin(sqrt(act)));                 % Algorithm 2 init

% Fig. 3c of the paper, in the same gene order (IRF4 REL PAX5 RELA PRDM1 AICDA)
paper = [ 0    -.070  0     0    .166 -.018;
         -.070  0    .078   0   -.129  .084;
          0    .078  0     .043 -.101  .067;
          0     0    .043   0     0     0  ;
         .166 -.129 -.101   0     0   -.078;
         -.018 .084  .067   0   -.078   0 ];
iu = triu(true(n),1); nz = paper(iu) ~= 0;
concord = @(th) sum(sign(th(iu)) == sign(paper(iu)) & nz);

% ---- Model A, paper order
fprintf('\n[A] original qscGRN, paper order\n');
[thA, LA] = qgrn_fit(@pout_qscgrn, theta0, pobs, m, struct('train_diag',false));
fprintf('[A] KL = %.2e   sign concordance with Fig.3c: %d/%d\n', LA, concord(thA), sum(nz));
disp(print_theta(thA, names));

% ---- Model A, random orderings (mapped back to paper order)
if exist('rng'), rng(1); else, rand('seed',1); end   % MATLAB / Octave
nperm = 4; stack = zeros(n,n,nperm);
for r = 1:nperm
    perm = randperm(n); inv(perm) = 1:n;
    [th, L] = qgrn_fit(@pout_qscgrn, diag(2*asin(sqrt(act(perm)))), build_pobs(Xb(perm,:)), m, ...
                       struct('train_diag',false,'verbose',false));
    stack(:,:,r) = th(inv,inv);
    fprintf('[A] perm %d (%s)  KL = %.2e\n', r, strjoin(names(perm)',' '), L);
end
rng_ = max(stack,[],3) - min(stack,[],3);
fprintf('[A] max range of any theta_kp across %d orderings: %.4f\n', nperm, max(rng_(iu)));

% ---- Model B, Ising Born machine, diagonals trained
mixes = [pi/6 pi/4 pi/3]; thB = cell(size(mixes));
for i = 1:numel(mixes)
    mix = mixes(i);
    fprintf('\n[B] Ising Born machine, mixer = %.3f, diagonals trained\n', mix);
    [thB{i}, LB] = qgrn_fit(@(t) pout_ising_bm(t, mix, Z), theta0, pobs, m, struct('train_diag',true));
    fprintf('[B] KL = %.2e   sign concordance with Fig.3c: %d/%d\n', LB, concord(thB{i}), sum(nz));
    disp(print_theta(thB{i}, names));
end
save('results_matlab.mat', 'thA', 'stack', 'thB', 'mixes', 'names', 'paper');
