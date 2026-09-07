function [theta, hist] = qgrn_train(theta0, pObs, ncells, opts)
%QGRN_TRAIN  Gradient descent on the QuantumGRN angles.
%
%   [THETA, HIST] = QGRN_TRAIN(THETA0, POBS, NCELLS, OPTS)
%
%   OPTS fields (defaults match the Python model class):
%     epochs        100      iterations
%     lr            0.1      learning rate
%     trainEncoder  false    if false the encoder angles stay at their
%                            activation-ratio initialisation
%     lossThreshold 1e-16    relative-plateau tolerance for early stopping
%     verbose       true
%
%   HIST has fields loss, cons, err and theta (epochs x n x n).
%
%   The stopping test mirrors the Python: over a trailing window of 10 epochs,
%   delta = |loss(e-10) - mean(loss(e-9:e))| / mean(...), and training stops
%   once delta drops below lossThreshold. With the 1e-16 default this
%   effectively means "run all epochs unless the loss stops changing at all".

    arguments
        theta0  (:,:) double
        pObs    (:,1) double
        ncells  (1,1) double
        opts.epochs        (1,1) double = 100
        opts.lr            (1,1) double = 0.1
        opts.trainEncoder  (1,1) logical = false
        opts.lossThreshold (1,1) double = 1e-16
        opts.verbose       (1,1) logical = true
    end

    theta = theta0;
    n = size(theta, 1);
    E = opts.epochs;
    window = 10;

    hist.loss  = zeros(E,1);
    hist.cons  = zeros(E,1);
    hist.err   = zeros(E,1);
    hist.theta = zeros(E, n, n);

    t0 = tic;
    for e = 1:E
        pOut = qgrn_pout(theta);

        [loss, cons1, err] = qgrn_loss(theta, pOut, pObs, ncells, 1);
        ratio = loss / cons1;          % keep the two terms comparable in size
        cons  = ratio * cons1;         % the constraint is linear in ratio

        hist.loss(e)      = loss;
        hist.cons(e)      = cons;
        hist.err(e)       = err;
        hist.theta(e,:,:) = theta;

        if e > window + 1
            iloss = mean(hist.loss(e-window+1:e));
            floss = hist.loss(e-window);
            delta = abs(floss - iloss) / iloss;
        else
            delta = opts.lossThreshold;
        end

        if opts.verbose && (mod(e, 10) == 0 || e == 1)
            fprintf('  epoch %3d/%d   KL = %.6e   err = %.3e\n', e, E, loss, err);
        end

        stop = (opts.lossThreshold > delta) || (e == E);
        if stop
            hist.loss  = hist.loss(1:e);
            hist.cons  = hist.cons(1:e);
            hist.err   = hist.err(1:e);
            hist.theta = hist.theta(1:e,:,:);
            if opts.verbose
                why = 'number of epochs reached';
                if opts.lossThreshold > delta, why = 'loss threshold reached'; end
                fprintf('  optimisation stopped after %d epochs (%s), %.1f s\n', ...
                        e, why, toc(t0));
            end
            break
        end

        G = qgrn_gradient(theta, pOut, pObs, ncells, ratio, opts.trainEncoder);
        theta = theta - opts.lr * G;
    end
end
