function [theta, Lkl, hist] = qgrn_fit(model, theta0, pobs, m, opts)
%QGRN_FIT  Gradient descent on the qscGRN loss with symmetrised gradient (Eq. 8).
%   model  : function handle theta -> pout
%   opts   : .iters (1200) .lr (0.05) .train_diag (false) .h (1e-6) .verbose
%   Gradient by central finite differences on the symmetric parameterisation:
%   perturbing theta_kp and theta_pk together gives d_kp + d_pk, which is
%   twice the (grad + grad')/2 entry used in the paper.
    if nargin < 5, opts = struct(); end
    d = struct('iters',1200,'lr',0.05,'train_diag',false,'h',1e-6,'verbose',true,'stop',-1);
    f = fieldnames(d); for i = 1:numel(f), if ~isfield(opts,f{i}), opts.(f{i}) = d.(f{i}); end, end
    n = size(theta0,1); theta = theta0;
    pobs_s = (pobs*m + 1) / (m + numel(pobs));
    [ku, pu] = find(triu(ones(n),1)); pairs = [ku pu];
    hist = zeros(opts.iters,1);
    for s = 1:opts.iters
        [~, Lkl, Lc] = qgrn_loss(theta, model(theta), pobs_s, m, 0);
        lam = Lkl / Lc;                            % dynamic lambda (paper)
        Lfun = @(th) qgrn_loss(th, model(th), pobs_s, m, lam);
        grad = zeros(n);
        for i = 1:size(pairs,1)
            k = pairs(i,1); p = pairs(i,2);
            tp = theta; tp(k,p) = tp(k,p)+opts.h; tp(p,k) = tp(p,k)+opts.h;
            tm = theta; tm(k,p) = tm(k,p)-opts.h; tm(p,k) = tm(p,k)-opts.h;
            g = (Lfun(tp) - Lfun(tm)) / (2*opts.h) / 2;
            grad(k,p) = g; grad(p,k) = g;
        end
        if opts.train_diag
            for k = 1:n
                tp = theta; tp(k,k) = tp(k,k)+opts.h;
                tm = theta; tm(k,k) = tm(k,k)-opts.h;
                grad(k,k) = (Lfun(tp) - Lfun(tm)) / (2*opts.h);
            end
        end
        theta = theta - opts.lr * grad;
        hist(s) = Lkl;
        if opts.verbose && mod(s,200) == 0, fprintf('  iter %4d  KL = %.3e\n', s, Lkl); end
        if opts.stop > 0 && Lkl < opts.stop, hist = hist(1:s); break; end
    end
    [~, Lkl] = qgrn_loss(theta, model(theta), pobs_s, m, 0);
end
