function [po,f1,statev]=i_fullcirc_symm(theta0,C)
    np=length(theta0);
    n=(1+sqrt(1+8*np))/2;
    assert(nchoosek(n,2)==np)

    % C must be the encoder followed by all n(n-1) coupling gates.
    assert(numel(C.Gates) == n + 2*np, ...
        'i_fullcirc_symm: expected %d gates (%d encoder + %d coupling), got %d.', ...
        n + 2*np, n, 2*np, numel(C.Gates));

    % Which shared parameter belongs to each unordered pair.
    a = nchoosek(1:n,2);
    pairIdx = zeros(n);
    for k = 1:np
        pairIdx(a(k,1),a(k,2)) = k;
        pairIdx(a(k,2),a(k,1)) = k;
    end

    % Rewrite every coupling gate in place, reading its OWN control and
    % target to decide which angle it takes. So the caller may lay the
    % coupling layer out control-major (Algorithm 4), pair-major, or any
    % other way, and each gate still gets the parameter for its own pair.
    %
    % This deliberately assumes nothing about slot order. The original code
    % did -- it wrote both directions of pair k to the same slot C.Gates(n+c)
    % with c advancing once per PAIR, so the second overwrote the first and
    % only 15 of the 30 coupling gates were ever updated; the other 15 kept
    % whatever angles C was constructed with. Nothing errored. A slot-index
    % scheme also has to change whenever the layout does, and the gate-count
    % assertion above cannot catch a layout mismatch because both layouts
    % have the same number of gates.
    for slot = 1:2*np
        g = C.Gates(n+slot);
        ctrl = g.ControlQubits(1);
        targ = g.TargetQubits(1);
        C.Gates(n+slot) = cryGate(ctrl, targ, theta0(pairIdx(ctrl,targ)));
    end
    S = simulate(C);
    [statev,po] = querystates(S);
    if nargout>1        
        f1=zeros(1,n);
        for k=1:n
            f1(k)=probability(S,k,"1");
        end
    end
end


% https://www.mathworks.com/matlabcentral/answers/1987423-change-of-parameters-of-gates-in-quantumcircuit-does-not-take-effect?s_tid=mlc_ans_email_view#answer_1265263