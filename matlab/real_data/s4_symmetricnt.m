%required input variables
% ---------------------------------------------------------
%f0=sum(y>0)./size(X,2);   % per gene initial activate freq.
%pt=n./size(X,2);          % target cell state freq.    
% ---------------------------------------------------------

%function s4_symmetricnt(f0,pt)

n=length(f0);
ncells=size(X,2);          % cells, for the paper's Laplace smoothing (alpha = 1)

layer1=[];
for k=1:n, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end

np=nchoosek(n,2);
% Start every coupling angle at 0, as the published method does (qscgrn's
% theta_init default is method="zeros"). This is not cosmetic. From a random
% start in +-pi/2, with bounds at +-pi and no L_cons barrier pulling the
% angles back, fmincon settles in a different basin with |theta| pinned at
% pi: paper-style KL 3.4e-01 and only 4/10 Fig. 3c signs. From zeros it
% recovers the published network -- KL 4.45e-03, 10/10 signs, max|theta|
% 0.166 (published: 4.2504e-03 and 0.166).
%
% For random multi-start instead, restore i_randpmpi(np,1) and keep the
% lowest-KL fit over several restarts -- but narrow that range, or tighten lb
% and ub to +-pi/2, or the large-angle solutions will keep coming back.
%
% NOTE these values only seed layer2 below; since the i_fullcirc_symm
% indexing fix, every coupling gate is rewritten on each objective call, so
% x0 acts purely as the optimiser's starting point.
x0=zeros(np,1);

% Coupling layer, one c-Ry per ordered pair. The layout here is "pair-major":
% both directions of genes (1,2), then both of (1,3), and so on. The two
% directions of a pair share one angle x0(c), which is the paper's
% symmetrisation theta_kp = theta_pk.
%
% (For reference only: Algorithm 4 of the paper emits these gates grouped by
% CONTROL instead -- every gate with control 1, then every gate with control
% 2. Because c-Ry gates with different controls do not commute, that is a
% slightly different circuit; it was measured at KL 4.2702e-03 against
% 4.2976e-03 for the layout below, but it also moved two near-threshold edges
% and is not used here.)
a=nchoosek(1:n,2);
layer2=[];
c=1;
for k=1:size(a,1)
    layer2=[layer2; cryGate(a(k,1),a(k,2),x0(c))];
    layer2=[layer2; cryGate(a(k,2),a(k,1),x0(c))];
    c=c+1;
end
C = quantumCircuit([layer1; layer2]);



methodid=1;


% --- the paper's L_cons barrier and its dynamic lambda -------------------
% Published loss:  L = L_KL + lambda*L_cons,  with
%     L_cons = sum 1/(theta^4 - (pi/2)^4)^2
% and lambda = L_KL/L_cons recomputed each iteration so the two terms stay the
% same order of magnitude.
%
% TRAP: lambda must be FROZEN while the solver runs. Recompute it from the
% same theta inside the objective and the barrier cancels identically,
%     lambda(th)*L_cons(th) = [L_KL(th)/L_cons(th)]*L_cons(th) = L_KL(th)
% so L collapses to 2*L_KL -- same argmin, constraint silently absent. The
% paper recomputes lambda per iteration but holds it constant within one, so:
% freeze it, solve, recompute at the new theta, repeat. It settles in a few
% rounds.
%
% NOTE L_cons is a singular SPIKE at |theta| = pi/2, not a wall: past that
% point theta^4-(pi/2)^4 grows again and the term decays. So it confines the
% angles only if the solver never steps across, which is why lb/ub are kept
% strictly inside pi/2 rather than at the old +-pi.
BOUND = (pi/2)*(1 - 1e-6);
nLambdaRounds = 4;

for round = 1:nLambdaRounds
    po_r = i_fullcirc_symm(x0,C);
    lambda = i_kldiverg(pt,po_r,true,ncells) / i_lcons(x0);   % frozen for this solve
    objFcn = @(x) i_obj(x,pt,C,lambda,ncells);

    switch methodid
        case 99
            % surrogateopt https://www.mathworks.com/help/matlab/math/ground-state-protein-folding-using-variational-quantum-eigensolver-vqe.html
            options = optimoptions("surrogateopt",...
                "MaxFunctionEvaluations",10, ...
                "PlotFcn","optimplotfval",...
                "InitialPoints",pi*ones(np,1));
        
            lb = repmat(-BOUND,np,1);
            ub = repmat(BOUND,np,1);
            [xa, fval] = surrogateopt(objFcn,lb,ub,[],[],[],[],[],options);
            % angles,minEnergy
        case 1
            options = optimoptions('fmincon','Display','none', ...
                'PlotFcns','optimplotfval','OutputFcn', @myoutput);
            lb=-BOUND*ones(np,1);
            ub=BOUND*ones(np,1);
            [xa,fval] = fmincon(objFcn,x0,[],[],[],[], ...
                lb,ub,[],options);
        case 2
            options = optimset('Display','none', ...
                'OutputFcn', @myoutput);        
            [xa,fval] = fminsearch(objFcn,x0,options);
        case 3
            options = optimset('Display','none','PlotFcns','optimplotfval', ...
                'TolX',1e-16, 'FunValCheck','on', 'MaxFunEvals',2e3, ...
                'MaxIter',2e3, 'TolFun',2e-6,'OutputFcn', @myoutput);
            [xa,fval] = fminunc(objFcn,x0,options);
    end
    fprintf('lambda round %d/%d: lambda = %.4e, fval = %.6f, max|theta| = %.4f\n', ...
            round, nLambdaRounds, lambda, fval, max(abs(xa)));
    x0 = xa;                                   % warm start the next round
end



[poa,f1a,statev]=i_fullcirc_symm(xa,C);



subplot(2,2,4)
bar([pt poa])
set(gca,'XTick',1:length(statev));
set(gca,'XTickLabel',statev);
ylabel('# of cells');
xlabel('Expression pattern');
legend({'Target','Observed'})
title(sprintf('kl=%f',i_kldiverg(pt,poa,true,ncells)));

At=zeros(n);
c=1;
a=nchoosek(1:n,2);
for k=1:size(a,1)
    At(a(k,1),a(k,2))=xa(c);
    At(a(k,2),a(k,1))=xa(c);
    c=c+1;
end


% At now holds the fitted angles theta_kp exactly as the optimiser returned
% them, which is what the paper reports (Fig. 3c). Removed a preceding
%   At=At.*abs(abs(At)-pi);
% which multiplied every edge by ~pi for small theta and by ~0 as
% |theta|->pi, so it rescaled all magnitudes nonlinearly.

figure;
subplot(2,2,1)
    % Published pruning rule: keep |theta| > pi/180 (~0.01745). That is the
    % absolute cutoff the shipped package applies in
    % qscgrn/visualization/qsc_grn.py:46 with its default threshold=1.
    % (Whether the paper text intends pi/180 or pi/360 is unsettled; the
    % package value is the one that is verifiable in code.)
    %
    % This replaced ten.e_filtadjc(At,0.65,false), which is a different
    % criterion: it renormalises At by the mean of its top 10 entries and
    % then keeps whatever exceeds the 65th percentile of |nonzeros|. Being a
    % quantile, it retains a roughly fixed FRACTION of edges however strong
    % or weak the fit is, rather than testing an absolute effect size.
    At2 = At;
    At2(abs(At2) <= pi/180) = 0;
    % graph, not digraph: theta_kp = theta_pk holds by construction here
    % (both directions of a pair share one angle), so At2 is symmetric and
    % digraph drew every edge twice, as a pair of opposing arrows implying a
    % directionality the symmetrised model cannot express.
    plot(graph(At2,genes,'omitselfloops'));   
    %imagesc(abs(A)+abs(A'))
    %imagesc(A); title('Ground truth')
subplot(2,2,2)
    %imagesc(abs(At)+abs(At'))
    imagesc(At); title('qGRN')
subplot(2,2,3)
    imagesc(sc_grn(X)); title('PCR')          % sc_pcnet -> sc_grn (default 'pcrnet')
subplot(2,2,4)
    bar([f0' f1a'])
    %set(gca,'XTick',1:n);
    %set(gca,'XTickLabel',statev);
    %ylabel('Gene');
    xlabel('Activation frequency');
    legend({'Target','Observed'})
    title(sprintf('kl=%f',i_kldiverg(f0,f1a)));


function [y]=i_obj(x,pt,C,lambda,ncells)
    % Published objective: the cell-state distribution term ONLY. The paper
    % fits KL(p_out || p_obs) and lets the c-Ry layers move the single-gene
    % marginals wherever they must -- the encoder pins them at
    % 2*asin(sqrt(f0)) only at initialisation, and Algorithm 2 freezes it
    % there rather than defending it.
    %
    % A second term used to be added here,
    %     y2 = i_kldiverg(f0, f1, false);
    % penalising exactly that drift. It is not in the paper and it fights the
    % first term: measured on this data, paper-style KL 4.4536e-03 with it
    % against 4.2990e-03 without (published 4.2504e-03). Signs were 10/10
    % either way, so it cost accuracy rather than correctness.
    %
    % f0 is gone from the signature: every solver now receives objFcn, a
    % closure carrying pt, C and the frozen lambda, so nothing is passed
    % positionally any more. The marginal drift is still reported in
    % subplot(2,2,4) -- now as a diagnostic of the fit rather than as
    % something the fit was told to minimise.
    po=i_fullcirc_symm(x,C);          % one output: skips the f1 marginals
    y=i_kldiverg(pt,po,true,ncells) + lambda*i_lcons(x);
end


function Lc = i_lcons(x)
%I_LCONS Paper's constraint term, Eq. 7: sum 1/(theta^4 - (pi/2)^4)^2.
%   Diverges as |theta| -> pi/2, which is what keeps the coupling angles
%   inside that box.
%
%   Summed over the free coupling angles only. The encoder angles are frozen
%   (Algorithm 2) so including them would add a constant, and counting each
%   pair once rather than once per gate direction only rescales Lc -- which
%   the L_KL/L_cons normalisation of lambda absorbs either way.
    Lc = sum(1./(x(:).^4 - (pi/2)^4).^2);
end


function stop = myoutput(x, optimValues, state,~,~,~)
stop = false;
fprintf('%.2f ',x);
fprintf('\n');
 %hold on;
 %plot(x(1),'.');
 %drawnow
end

