function psi = apply_c1q(psi, U, c, t, n)
%APPLY_C1Q  Controlled-U: apply U to target t on the subspace where control c = 1.
    idx = 0:2^n-1;
    sel = bitget(idx, c) == 1;                    % states with control bit set
    sub = psi(sel);                               % (n-1)-qubit sub-vector, same ordering
    tt = t - (t > c);                             % target index after removing control axis
    sub = apply_1q(sub, U, tt, n-1);
    psi(sel) = sub;
end
