function p = pout_qscgrn(theta)
%POUT_QSCGRN  Model A: original qscGRN circuit (Roman-Vicharra & Cai 2023).
%   Ry(theta_kk) encoder, then regulation layers L_1..L_n in order:
%   for control k, for target p ~= k: c-Ry(theta_kp). Order-sensitive.
    n = size(theta,1); G = gate_defs();
    psi = zeros(2^n,1); psi(1) = 1;
    for k = 1:n, psi = apply_1q(psi, G.Ry(theta(k,k)), k, n); end
    for k = 1:n
        for p = 1:n
            if p ~= k, psi = apply_c1q(psi, G.Ry(theta(k,p)), k, p, n); end
        end
    end
    p = abs(psi).^2; p(1) = 0; p = p / sum(p);
end
