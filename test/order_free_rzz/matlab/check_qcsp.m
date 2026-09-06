% CHECK_QCSP  Cross-check the plain statevector models against the MATLAB
%   Quantum Computing Support Package at the fitted parameters, and export
%   the Model B circuit to OpenQASM. Requires MATLAB + the support package.
%   Run from test/order_free_rzz/matlab/ after run_comparison.m.
clear; clc;
assert(exist('quantumCircuit','class') == 8 || exist('quantumCircuit','file') > 0, ...
    'MATLAB Support Package for Quantum Computing not found');
R = load(fullfile('..','results','results_matlab.mat'));
n = numel(R.names); Z = zz_signs(n);

for bo = {'q1lsb', 'q1msb'}
    bitorder = bo{1};
    dA = max(abs(pout_qscgrn(R.thA) - pout_qscgrn_qcsp(R.thA, bitorder)));
    dB = zeros(1, numel(R.mixes));
    for i = 1:numel(R.mixes)
        dB(i) = max(abs(pout_ising_bm(R.thB{i}, R.mixes(i), Z) - ...
                        pout_ising_bm_qcsp(R.thB{i}, R.mixes(i), bitorder)));
    end
    fprintf('bitorder %s : max|dp|  Model A = %.2e   Model B = %s\n', ...
            bitorder, dA, mat2str(dB, 3));
    if dA < 1e-10 && all(dB < 1e-10)
        fprintf('  -> plain statevector code and support package AGREE (convention: %s)\n', bitorder);
        ok = true;
    end
end
if ~exist('ok','var')
    warning(['No bit ordering gives agreement. Likely cause: rzzGate/ryGate angle ' ...
             'convention differs from exp(-i*theta/2*G). Compare a 2-qubit case by hand.']);
end

% Export the fitted Model B (mixer pi/6) circuit for hardware / Qiskit
[~, C] = pout_ising_bm_qcsp(R.thB{1}, R.mixes(1));
qasm = generateQASM(C);
fid = fopen(fullfile('..','results','ising_bm_mix_pi6.qasm'), 'w'); fwrite(fid, qasm); fclose(fid);
fprintf('Exported ../results/ising_bm_mix_pi6.qasm  (%d gates, %d qubits)\n', numel(C.Gates), C.NumQubits);
