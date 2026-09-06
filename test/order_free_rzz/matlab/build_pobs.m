function pobs = build_pobs(Xb)
%BUILD_POBS  Observed distribution over 2^n activation states.
%   Xb: n x m binary matrix, row k = qubit k. Qubit 1 is the least
%   significant bit of the state index (MATLAB column-major convention,
%   consistent with apply_1q / zz_signs). |0>^n is set to 0 and rescaled.
    n = size(Xb,1); m = size(Xb,2);
    lab = zeros(1, m);
    for k = 1:n, lab = lab + Xb(k,:) * 2^(k-1); end
    pobs = accumarray(lab(:)+1, 1, [2^n 1]) / m;
    pobs(1) = 0; pobs = pobs / sum(pobs);
end
